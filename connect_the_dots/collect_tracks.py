import trackpy
from itertools import product
import pandas as pd
from collections import Iterable
import copy
import h5py
import pandas as pd
import re

from operator import mul
from functools import reduce

import scipy.stats as sp
from scipy.stats import mannwhitneyu

import matplotlib.pyplot as plt
import numpy as np

import bioformats
import javabridge # for interfacing with java (required for bioformats)
javabridge.start_vm(class_path=bioformats.JARS)

###############################################################################

###############################################################################
def get_CZI_voxel(czi_file,channel,frame,z,y,x,windowZYXdims,img_info=None):
    """
    Parameters
    ----------
    czi_file : string
        File path to czi file containing z-stack for a single frame and channel
        of the 3D image time-series.

    channel : int
        Integer specifying the channel for which to obtain zstack.

    frame : int
        Integer specifying the frame for which to obtain zstack.

    z, y ,x : float
        Coordinates of the center (best guess) of the dot positions in the z-stack.

    windowZYXdims : tuple of ints
        Lengths of the "voxel" around `(z, y, x)` coordinates for which to crop out data.
        `windowZYXdims` dimensions correspond to 1/2 the total window length around
        the chosen coordinates for each dimension.

    img_info : list of 5-tuples of ints, optional
        A list of tuples contains size-related metadata for each of the
        files in `in_filenames`. The following image values:
        (sizeX,sizeY,sizeZ,sizeT,num_channels).

    Output
    ------
    zstack : ndarray of floats
        A cutout of the larger 3D volume with maximum dimensions of `windowZYXdims`. The size is smaller than `windowZYXdims` if the
        `z`, `y` , `x` are too close to the boundary.

    """
#     # start java virtual machine (used by bioformats to read czi files)
#     javabridge.start_vm(class_path=bioformats.JARS)
    
    # get image metadata
    if img_info is None:
        metadata = xml2dict(bioformats.get_omexml_metadata(czi_file))
        sizeT = metadata['OME']['Image']['Pixels']['SizeT']
        sizeX = metadata['OME']['Image']['Pixels']['SizeX']
        sizeY = metadata['OME']['Image']['Pixels']['SizeY']
        sizeZ = metadata['OME']['Image']['Pixels']['SizeZ']
        num_channels = len(metadata['OME']['Image']['Pixels']['Channel'])
    else:
        (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info

    minX = int(np.max([x-windowZYXdims[2],0]))
    minY = int(np.max([y-windowZYXdims[1],0]))
    minZ = int(np.max([z-windowZYXdims[0],0]))
    maxX = int(np.min([x+windowZYXdims[2],sizeX]))
    maxY = int(np.min([y+windowZYXdims[1],sizeY]))
    maxZ = int(np.min([z+windowZYXdims[0],sizeZ]))
    
    zstack = np.zeros((maxZ-minZ,maxY-minY,maxX-minX))
       
    # get zstack
    with bioformats.ImageReader(czi_file) as reader:
        for z in np.arange(minZ,maxZ):
            zstack[z-minZ,:,:] = reader.read(t=frame,z=z,c=channel)[minY:maxY,minX:maxX]
            
#     # kill java virtual machine
#     javabridge.kill_vm()
 
    return zstack, (minZ, minY, minX)


def getHDF5ZstackVoxel(hdf5_file,channel,frame,z,y,x,windowZYXdims):
    """
    Parameters
    ----------
    hdf5_file : string
        File path to HDF5 file containing z-stack for a single frame of the
        3D image time-series.The HDF5 file must have zstacks stored in
        datasets with the following key naming convention: "frame_{}_channel_{}"

    channel : int
        Integer specifying the channel for which to obtain zstack.

    frame : int
        Integer specifying the frame for which to obtain zstack.

    z, y ,x : float
        Coordinates of the center (best guess) of the dot positions in the z-stack.

    windowZYXdims : tuple of ints
        Lengths of the "voxel" around `(z, y, x)` coordinates for which to crop out data.
        `windowZYXdims` dimensions correspond to 1/2 the total window length around
        the chosen coordinates for each dimension.

    Output
    ------
    zstack : ndarray of floats
        A cutout of the larger 3D volume with maximum dimensions of `windowZYXdims`. The size is smaller than `windowZYXdims` if the
        `z`, `y` , `x` are too close to the boundary.

    (minZ, minY, minX) : tuple of ints
        Specifies the absolute minimum coorinates of the `zstack` bounding box.

    """
    zstack = 0
    with h5py.File(hdf5_file,'r') as h:
        zstack = h['frame_{}_channel_{}'.format(frame,channel)]
        dims = zstack.shape
        minX = int(np.max([x-windowZYXdims[2],0]))
        minY = int(np.max([y-windowZYXdims[1],0]))
        minZ = int(np.max([z-windowZYXdims[0],0]))
        maxX = int(np.min([x+windowZYXdims[2],dims[2]]))
        maxY = int(np.min([y+windowZYXdims[1],dims[1]]))
        maxZ = int(np.min([z+windowZYXdims[0],dims[0]]))
        zstack = h['frame_{}_channel_{}'.format(frame,channel)][minZ:maxZ, minY:maxY, minX:maxX]
        #print(dims)
    return zstack, (minZ, minY, minX)


def create_trajectories(table_filename, data_filename, \
            save_fig_name_before_filtering, save_fig_name_after_filtering, \
            dot_trajectory_table, \
            search_range=(5,20,20), \
            memory=3, min_track_length = 100,  \
            pval_min = 1e-5,verbose=False):
    """
    Creates dot trajectories based on pre-identified positions.
    Dot positions are linked into trajectories (within and across channels).
    Frames with missing dots data are "restored" and dot positions are refined.


    Parameters
    ----------
    table_filename : str
        Path to the file containing the object feature predictions csv file.

    data_filename : str
        Path to the file containing the HDF5-formatted image series data.

    save_fig_name_before_filtering : str
        Path/name of the file for which to output figure of dot tracks.

    save_fig_name_after_filtering : str
        Path/name of the file for which to output figure of dot tracks.

    search_range : tuple of ints
        Maximum search range in pixels for the "next" dot for creating particle
        trajectories. Dimensions are (Z, Y, X).

    memory : int
        The number of frames for which particles can disappear from the image
        segmentation before a particle trajectory is terminated.

    min_track_length : int
        The minimum number of frames required for processing dot trajectories.
        This value helps also serve as a filter for false-positive detections.

    pval_min : float
        Minimum p-value from which to filter out step size distributions for
        which the distributions across channel0 and channel1 are dissimilar.
        The metrics used are both the mannwhitneyu and the f_test.

    verbose : T/F
        Displays text describing function progress.

    Output
    ------
    Two figures corresponding to: `save_fig_name_before_filtering` and
    `save_fig_name_after_filtering`.

    DataFrame of

    """
    colour_scheme = ['r','b','g','m','c']
    windowZYXdims = search_range

    # load csv files with segmentation positions ouputs
    df_particles = pd.read_csv(table_filename)
    df_particles['timestep'] = df_particles['frame']

    # get special positions
    c = df_particles['channel'].values
    t = df_particles['timestep'].values
    x = df_particles['Center of the object_0'].values
    y = df_particles['Center of the object_1'].values
    z = df_particles['Center of the object_2'].values

    # make new, reduced data frame
    df_short = pd.DataFrame({'channel': c, 'frame': t,'x': x, 'y': y, 'z': z})

    # do for each channel
    channels = set(df_short['channel'].values)


    # link trajectories across channels:


    # artificially shift the frames and create trajectories
    linked_df = df_short.copy()
    linked_df['frame'] = linked_df['frame']*len(channels)+linked_df['channel']

    linked_df = trackpy.link_df(linked_df, \
                search_range=tuple([x*2 for x in search_range]), \
                memory=memory*2)

    # reverse shift the frame numbers
    linked_df['frame'] = (linked_df['frame']-linked_df['channel'])/len(channels)

    # get unique IDs for sets of linked dots
    particle_set = set(linked_df.particle)

    # average out the linked dot positions in each channel
    # this is used as the "initial guess" for the "refined" dot positions
    avg_linked_df = linked_df.groupby(['particle','frame']).mean()
    avg_linked_df.reset_index(inplace=True)

    refined_df_list = []

    plt.figure(figsize=(10,10))
    for ch in list(channels)[::-1]:
        if verbose == True:
            print("Doing channel {}".format(ch))

        # 1) go through each frame in tracks_df,
        # 2) crop out z-stack around predicted dot positions
        # 3) refine the dot positions
        # 4) fill in dot positions for missing particles

        for p in particle_set:

            if len(linked_df[linked_df.particle == p]) < min_track_length:
                continue

            # get average positions
            X = avg_linked_df[avg_linked_df.particle == p].x.values
            Y = avg_linked_df[avg_linked_df.particle == p].y.values
            Z = avg_linked_df[avg_linked_df.particle == p].z.values
            T = avg_linked_df[avg_linked_df.particle == p].frame.values
            T = np.array(T, dtype=int)

            if len(T) == 0:
                continue

            # get index of sorted frames
            std = np.argsort(T)

            # plot the trajectories "before"
            plt.plot(X,Y,'k:')
            plt.title("Before filtering")

            # refined positions
            traj_X = []
            traj_Y = []
            traj_Z = []
            traj_T = []
            near_edge = []

            count = 0
            for frame in np.arange(T[std[0]],T[std[-1]]):

                # get dot position estimate
                if frame == T[std[count]]:
                    count += 1
                    # use known current location
                    z, y, x = (Z[std[count]], Y[std[count]], X[std[count]])
                else:
                    # use last known location
                    z, y, x = (Z[std[count]], Y[std[count]], X[std[count]])

                # extract voxel surrounding (estimated) dot location
                zstack, (z_low, y_low, x_low) = getHDF5ZstackVoxel(data_filename,channel=ch,frame=frame,\
                                            x=x,y=y,z=z,windowZYXdims=windowZYXdims)

                # classify dots as being "single" or "double"

                # classify dots as being near the edge otherwise
                is_near_edge = (reduce(mul,zstack.shape)/reduce(mul,windowZYXdims) < 0.7)
                near_edge.append(is_near_edge)

                # refine the dot localization estimate
                z_rel, y_rel, x_rel = np.unravel_index(np.argmax(zstack),zstack.shape)


                z_refined = z_low + z_rel
                y_refined = y_low + y_rel
                x_refined = x_low + x_rel

                traj_X.append(x_refined)
                traj_Y.append(y_refined)
                traj_Z.append(z_refined)
                traj_T.append(frame)

                # classify dots as being "uni-" or "multi-" modal

            tmp_df = pd.DataFrame.from_dict({'particle':[p]*len(traj_T),\
                                             'channel':[ch]*len(traj_T), 'frame':traj_T, \
                                             'x':traj_X, 'y':traj_Y, 'z':traj_Z, \
                                            'near_edge':near_edge})
            refined_df_list.append(tmp_df)

            # plot the trajectories "after"
            plt.plot(traj_X,traj_Y,colour_scheme[np.mod(ch,len(colour_scheme))])
    plt.savefig(save_fig_name_before_filtering)
    if verbose == True:
        print("Done creating rough trajectories")

    # filter good from bad dot trajectories based on step size stats.
    # use mean and variance statistics on the stepsize distributions
    # to see if the ditributions are similar enough in each channel
    refined_df = pd.concat(refined_df_list)
    plt.figure(figsize=(10,10))
    for p in particle_set:
        step_size_dist = []
        X_list = []
        Y_list = []
        for ch in list(channels)[::-1]:
            X = refined_df[(refined_df.particle == p) & \
                            (refined_df.channel == ch)].x.values
            Y = refined_df[(refined_df.particle == p) & \
                            (refined_df.channel == ch)].y.values
            Z = refined_df[(refined_df.particle == p) & \
                            (refined_df.channel == ch)].z.values
            X_list.append(X)
            Y_list.append(Y)
            step_size_dist.append(np.sqrt(np.diff(X)**2 + np.diff(Y)**2 + np.diff(Y)**2))
        try:
            ustat, upval = sp.mannwhitneyu(step_size_dist[0],step_size_dist[1])
            fstat, fpval = sp.f_oneway(step_size_dist[0],step_size_dist[1])

            if upval > pval_min and fpval > pval_min:
                for ch in list(channels)[::-1]:
                    plt.plot(X_list[ch],Y_list[ch],\
                            colour_scheme[np.mod(ch,len(colour_scheme))])
        except:
            1
    plt.savefig(save_fig_name_after_filtering)

    refined_df.to_csv(dot_trajectory_table,index=False)

        # classify dots as being near the edge, on the edge, or otherwise

        # return "fuller" data frame of tracks

        # make plot of trajectories with minimum length
        # make histogram of gaps per trajectory (before and after filtering)
        
def create_trajectories_from_czi(table_filename, data_filename, \
            save_fig_name_before_filtering, save_fig_name_after_filtering, \
            dot_trajectory_table, \
            search_range=(5,20,20), \
            memory=3, min_track_length = 100,  \
            pval_min = 1e-5,verbose=False, img_info=None):
    """
    Creates dot trajectories based on pre-identified positions.
    Dot positions are linked into trajectories (within and across channels).
    Frames with missing dots data are "restored" and dot positions are refined.


    Parameters
    ----------
    table_filename : str
        Path to the file containing the object feature predictions csv file.

    data_filename : str
        Path to the file containing the HDF5-formatted image series data.

    save_fig_name_before_filtering : str
        Path/name of the file for which to output figure of dot tracks.

    save_fig_name_after_filtering : str
        Path/name of the file for which to output figure of dot tracks.

    search_range : tuple of ints
        Maximum search range in pixels for the "next" dot for creating particle
        trajectories. Dimensions are (Z, Y, X).

    memory : int
        The number of frames for which particles can disappear from the image
        segmentation before a particle trajectory is terminated.

    min_track_length : int
        The minimum number of frames required for processing dot trajectories.
        This value helps also serve as a filter for false-positive detections.

    pval_min : float
        Minimum p-value from which to filter out step size distributions for
        which the distributions across channel0 and channel1 are dissimilar.
        The metrics used are both the mannwhitneyu and the f_test.

    verbose : T/F
        Displays text describing function progress.
        
        
    img_info : list of 5-tuples of ints, optional
        A list of tuples contains size-related metadata for each of the raw data
        files in `data_filename`. The following image values:
        (sizeX,sizeY,sizeZ,sizeT,num_channels).

    Output
    ------
    Two figures corresponding to: `save_fig_name_before_filtering` and
    `save_fig_name_after_filtering`.

    DataFrame of

    """
    colour_scheme = ['r','b','g','m','c']
    windowZYXdims = search_range

    # load csv files with segmentation positions ouputs
    df_particles = pd.read_csv(table_filename)
    df_particles['timestep'] = df_particles['frame']

    # get special positions
    c = df_particles['channel'].values
    t = df_particles['timestep'].values
    x = df_particles['Center of the object_0'].values
    y = df_particles['Center of the object_1'].values
    z = df_particles['Center of the object_2'].values

    # make new, reduced data frame
    df_short = pd.DataFrame({'channel': c, 'frame': t,'x': x, 'y': y, 'z': z})

    # do for each channel
    channels = set(df_short['channel'].values)

    # obtain image metadata
    if img_info is None:
        img_info, _ = get_CZI_metadata(data_filename)

    # artificially shift the frames and create trajectories
    linked_df = df_short.copy()
    linked_df['frame'] = linked_df['frame']*len(channels)+linked_df['channel']

    linked_df = trackpy.link_df(linked_df, \
                search_range=tuple([x*2 for x in search_range]), \
                memory=memory*2)

    # reverse shift the frame numbers
    linked_df['frame'] = (linked_df['frame']-linked_df['channel'])/len(channels)

    # get unique IDs for sets of linked dots
    particle_set = set(linked_df.particle)

    # average out the linked dot positions in each channel
    # this is used as the "initial guess" for the "refined" dot positions
    avg_linked_df = linked_df.groupby(['particle','frame']).mean()
    avg_linked_df.reset_index(inplace=True)

    refined_df_list = []

    plt.figure(figsize=(10,10))
    for ch in list(channels)[::-1]:
        if verbose == True:
            print("Doing channel {}".format(ch))

        # 1) go through each frame in tracks_df,
        # 2) crop out z-stack around predicted dot positions
        # 3) refine the dot positions
        # 4) fill in dot positions for missing particles

        for p in particle_set:

            if len(linked_df[linked_df.particle == p]) < min_track_length:
                continue

            # get average positions
            X = avg_linked_df[avg_linked_df.particle == p].x.values
            Y = avg_linked_df[avg_linked_df.particle == p].y.values
            Z = avg_linked_df[avg_linked_df.particle == p].z.values
            T = avg_linked_df[avg_linked_df.particle == p].frame.values
            T = np.array(T, dtype=int)

            if len(T) == 0:
                continue

            # get index of sorted frames
            std = np.argsort(T)

            # plot the trajectories "before"
            plt.plot(X,Y,'k:')
            plt.title("Before filtering")

            # refined positions
            traj_X = []
            traj_Y = []
            traj_Z = []
            traj_T = []
            near_edge = []

            count = 0
            for frame in np.arange(T[std[0]],T[std[-1]]):

                # get dot position estimate
                if frame == T[std[count]]:
                    count += 1
                    # use known current location
                    z, y, x = (Z[std[count]], Y[std[count]], X[std[count]])
                else:
                    # use last known location
                    z, y, x = (Z[std[count]], Y[std[count]], X[std[count]])

                # extract voxel surrounding (estimated) dot location
                zstack, (z_low, y_low, x_low) = get_CZI_voxel(data_filename,channel=ch,frame=frame,\
                                            x=x,y=y,z=z,windowZYXdims=windowZYXdims,img_info=img_info)

                # classify dots as being "single" or "double"

                # classify dots as being near the edge otherwise
                is_near_edge = (reduce(mul,zstack.shape)/reduce(mul,windowZYXdims) < 0.7)
                near_edge.append(is_near_edge)

                # refine the dot localization estimate
                z_rel, y_rel, x_rel = np.unravel_index(np.argmax(zstack),zstack.shape)


                z_refined = z_low + z_rel
                y_refined = y_low + y_rel
                x_refined = x_low + x_rel

                traj_X.append(x_refined)
                traj_Y.append(y_refined)
                traj_Z.append(z_refined)
                traj_T.append(frame)

                # classify dots as being "uni-" or "multi-" modal

            tmp_df = pd.DataFrame.from_dict({'particle':[p]*len(traj_T),\
                                             'channel':[ch]*len(traj_T), 'frame':traj_T, \
                                             'x':traj_X, 'y':traj_Y, 'z':traj_Z, \
                                            'near_edge':near_edge})
            refined_df_list.append(tmp_df)

            # plot the trajectories "after"
            plt.plot(traj_X,traj_Y,colour_scheme[np.mod(ch,len(colour_scheme))])
    plt.savefig(save_fig_name_before_filtering)
    if verbose == True:
        print("Done creating rough trajectories")

    # filter good from bad dot trajectories based on step size stats.
    # use mean and variance statistics on the stepsize distributions
    # to see if the ditributions are similar enough in each channel
    refined_df = pd.concat(refined_df_list)
    plt.figure(figsize=(10,10))
    for p in particle_set:
        step_size_dist = []
        X_list = []
        Y_list = []
        for ch in list(channels)[::-1]:
            X = refined_df[(refined_df.particle == p) & \
                            (refined_df.channel == ch)].x.values
            Y = refined_df[(refined_df.particle == p) & \
                            (refined_df.channel == ch)].y.values
            Z = refined_df[(refined_df.particle == p) & \
                            (refined_df.channel == ch)].z.values
            X_list.append(X)
            Y_list.append(Y)
            step_size_dist.append(np.sqrt(np.diff(X)**2 + np.diff(Y)**2 + np.diff(Y)**2))
        try:
            ustat, upval = sp.mannwhitneyu(step_size_dist[0],step_size_dist[1])
            fstat, fpval = sp.f_oneway(step_size_dist[0],step_size_dist[1])

            if upval > pval_min and fpval > pval_min:
                for ch in list(channels)[::-1]:
                    plt.plot(X_list[ch],Y_list[ch],\
                            colour_scheme[np.mod(ch,len(colour_scheme))])
        except:
            1
    plt.savefig(save_fig_name_after_filtering)

    refined_df.to_csv(dot_trajectory_table,index=False)

        # classify dots as being near the edge, on the edge, or otherwise

        # return "fuller" data frame of tracks

        # make plot of trajectories with minimum length
        # make histogram of gaps per trajectory (before and after filtering)
