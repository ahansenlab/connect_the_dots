## Input/output functions. 
# Written by Hugo B. Brandao (hbrandao@g.harvard.edu)
# (c) 2021, Hugo B. Brandao

## Import required modules
import numpy as np # for manipulating arrays
import os # for making/deleting directories
import bioformats # for reading image series
import javabridge # for interfacing with java (required for bioformats)
from tifffile import xml2dict # for parsing the metadata from bioformats
import pickle # for saving python objects and other data
from pathlib import Path

# for plotting images
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use("Agg")
from matplotlib.patches import Rectangle
from matplotlib import animation
import pandas as pd

import gc # trying to get garbage collection to work with matplotlib
import time

javabridge.start_vm(class_path=bioformats.JARS) # start java virtual machine


def search_for(start_folder, name_includes=None, name_excludes=None):
    """
    Search directories for files containing or excluding specific strings.
    
    Parameters
    ----------
    start_folder : str
        Parent directory for recursive filename search.
    
    name_includes : list of strings
        List of strings to include in the search data. Search does not do case matching. 

    name_excludes : list of strings
        List of strings to exclude from the search. Search does not do case matching. 
        
    Returns
    -------
    Two lists: 1) list of file names, 2) list of the absolute file path.
    
    """
    filenames_list = []
    filepath_list = []
    for path in Path(start_folder).rglob('*.*'):

        parent_folder = path.parent.name
        if parent_folder in name_excludes:
            continue   

        if (all([name.lower() in path.name.lower() for name in name_includes])==False) or \
            any([name.lower() in path.name.lower() for name in name_excludes])==True:
            continue
            
        filenames_list.append(path.name)
        filepath_list.append(str(path.parent))
    return filenames_list, filepath_list


def get_CZI_metadata(filename,filepath=None,verbose=False):
    """
    Obtains the metadata from a CZI image series.
    
    Parameters
    ----------
    filename : str
               Name of the file from which to retrieve the z-stack.
            
    filepath : str, optional
               Path to the file.
           
    verbose : {T,F}, optional
              If true, prints (sizeX,sizeY,sizeZ,sizeT,num_channels) to standard output
            
    Returns
    -------
    (sizeX,sizeY,sizeZ,sizeT,num_channels) : tuple of ints
                Information on the length of the sizes of the `X`, `Y`, `Z` (spatial) and `T` 
                (temporal) dimensions of the image series and the number of channels, `num_channels`.
                In case of failutre to load, returns a 5-tuple of values 0. 

    metadata : dict, or None
               Dictionary containing the full metadata formatted in the Bioformats OME style. 
               If loading is unsuccessful, `None` is returned.
    """
    
    if not filepath is None:
        czi_image = os.path.join(filepath,filename)
    else:
        czi_image = filename

    if not os.path.exists(czi_image):
        return (0,0,0,0,0), None
        
    metadata = xml2dict(bioformats.get_omexml_metadata(czi_image))
    sizeT = metadata['OME']['Image']['Pixels']['SizeT']
    sizeX = metadata['OME']['Image']['Pixels']['SizeX']
    sizeY = metadata['OME']['Image']['Pixels']['SizeY']
    sizeZ = metadata['OME']['Image']['Pixels']['SizeZ']
    num_channels = len(metadata['OME']['Image']['Pixels']['Channel'])
    
    if verbose:
        print(sizeX,sizeY,sizeZ,sizeT,num_channels)
    
    return (sizeX,sizeY,sizeZ,sizeT,num_channels), metadata


def get_CZI_zstack(filename,frame,channel,filepath=None,img_info=None):
    """
    Obtains a single z-stack from a 3D imaging time-series for a specified time and channel.
    
    Parameters
    ----------
    filename : str
               Name of the file from which to retrieve the z-stack.

    frame : int
               The temporal slice of the image series from which to retrieve the z-stack.
            
    channel : int
               The channel from which to retrieve the z-stack.
            
    filepath : str, optional
               Path to the file.

    img_info : tuple of ints, optional 
               5-tuple containing lengths of the `X`, `Y`, `Z` (spatial), `T` (temporal) dimensions
               of the image series, and the number of channels, `num_channels`.
               E.g. (sizeX,sizeY,sizeZ,sizeT,num_channels). See output of get_CZI_metadata(). 
               Pass these pre-computed values for increased speed in batch processing.
                
    Returns
    -------
    zstack : numpy.ndarray, or None
             Z-stack of the image series specified by the desired `frame`; contains 3 spatial 
             dimensions. If loading is unsuccessful, `None` is returned.
    """
    # prepare file name, check that file exists 
    if not (filepath is None):
        czi_image = os.path.join(filepath,filename)
    else:
        czi_image = filename
    if not os.path.exists(czi_image):
        return None
    
    # retrieve image dimensions, and number of channels
    if img_info is None:
        (sizeX,sizeY,sizeZ,sizeT,num_channels), _ = get_CZI_metadata(filename,filepath=filepath)
    else:
        assert len(img_info) == 5
        (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info
        
    
    # make sure frame and channel are in bounds
    assert frame < sizeT
    assert channel < num_channels
        
    #initialize array and load z-stack
    zstack = np.zeros((sizeZ, sizeY,sizeX))
    with bioformats.ImageReader(czi_image) as reader:
        for z in range(sizeZ):
            zstack[z,:,:] = reader.read(t=frame,z=z,c=channel)

    return zstack

def get_CZI_zstack_timeseries(filename,frames,channel,filepath=None,img_info=None):
    """
    Obtains a time-series of z-stacks for a specified channel from .CZI data.
    
    Parameters
    ----------
    filename : str
               Name of the file from which to retrieve the z-stack.

    frames : range or list of ints
               The temporal slices of the image series from which to retrieve the z-stack.
            
    channel : int
               The channel from which to retrieve the z-stack.
            
    filepath : str, optional
               Path to the file.

    img_info : tuple of ints, optional 
               5-tuple containing lengths of the `X`, `Y`, `Z` (spatial), `T` (temporal) dimensions
               of the image series, and the number of channels, `num_channels`.
               E.g. (sizeX,sizeY,sizeZ,sizeT,num_channels). See output of get_CZI_metadata(). 
               Pass these pre-computed values for increased speed in batch processing.
                
    Returns
    -------
    zstack : numpy.ndarray, or None
             Z-stack of the image time-series; contains 4 spatial dimensions. 
             If loading is unsuccessful, `None` is returned. Dimensions are (Z,Y,X,T).
             
    """
    # prepare file name, check that file exists 
    if not (filepath is None):
        czi_image = os.path.join(filepath,filename)
    else:
        czi_image = filename
    if not os.path.exists(czi_image):
        return None
    
    # retrieve image dimensions, and number of channels
    if img_info is None:
        (sizeX,sizeY,sizeZ,sizeT,num_channels), _ = get_CZI_metadata(filename,filepath=filepath)
    else:
        assert len(img_info) == 5
        (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info
    
    # make sure frame and channel are in bounds
    for frame in frames:
        assert frame < sizeT
    assert channel < num_channels
        
    #initialize array and load z-stack
    num_slices = len(frames)
    zstack = np.zeros((sizeZ, sizeY,sizeX,num_slices))
    with bioformats.ImageReader(czi_image) as reader:
        for frame in frames:
            for z in range(sizeZ):
                zstack[z,:,:,frame] = reader.read(t=frame,z=z,c=channel)

    return zstack    


# get metadata from image
def make_movie(Z_dict,
               linked_df,
               metadata,               
               output_filename='test.mp4',
               desired_bar_length_um=1,
               percent_cutoff=99.99,
               disp_line_len=25,
               verbose=False,
               max_axis=0,
               text_dist=10,
               adaptive_text=False,
               line_alpha=0.5,
               adaptive_LUT=True,
               millisecond_per_frame=200):
    """
    Makes XY maximum intensity projection movies of the data and the identified trajectories.
    
    Parameters
    ----------
    Z_dict : dict of np.ndarrays
        Dictionary of the z-stack timeseries data
    
    linked_df : pandas.DataFrame
        DataFrame of dot locations for valid trajectories. If None is provided, no trajectories are overlayed.

    metadata : dict
        Dictionary containing the OME-formatted metadata for the timeseries data. See
        e.g. `connect_the_dots.io.get_CZI_metadata()`.

    output_filename : str
        Name of the output file. Should contain the extension .mp4      
        
    desired_bar_length_um : float
        Desired length of the scale bar in microns.
        
    percent_cutoff : float
        Intensity percentile cutoff (0-100) for the normalization of the time-series image. 
        Normalization is done based on the maximum intensity projection of the first image.
    
    disp_line_len : int
        Number of time-points to display for the trajectory line.
        
    verbose : boolean {T,F}
        Print progress for the creation of the movie to standard output.
        
    max_axis : int
        Maximum intensity projection axis to display.
        
    text_dist : float
        Number of pixels away from the dot localizations whereby to place the dot trajectory ID.
        
    adaptive_text : boolean {T,F}
        Adjusts the text position to be less overlapping with the trajectory (if True). Otherwise,
        the dot trajectory ID is at a fixed position.
    
    line_alpha : float
        Value for the transparency of the dot trajectory lines.
        
    adaptive_LUT : boolean {T,F}
        Linear adjustment the image threshold (if True). Otherwise, fix to the first image. 
    
    Results
    -------
    movie :
    
    
    """
    sizeZ, sizeY, sizeX, sizeT = Z_dict[0].shape
    
    bar_pixel_width_um = metadata['OME']['Image']['Pixels']['PhysicalSizeX']
    bar_width = desired_bar_length_um/bar_pixel_width_um # for length in pixels/um
    bar_height = 5 # in pixels

    # if linked_df is None, create a dummy list
    if linked_df is None:
        linked_df = pd.DataFrame({'x': {},
         'y': {},
         'z': {},
         'mean_intensity': {},
         'max_intensity': {},
         'dot_size_in_pixels': {},
         'frame': {},
         'channel': {},
         'particle': {}})
    # make initial image
    fig = plt.figure(figsize=(10,10))
    frame = 0
    channel = 0
    img = np.max(Z_dict[0][:,:,:,0],axis=max_axis)#get_zstack_max(os.path.join(input_filepath,input_file),0,0)
    img2 = np.max(Z_dict[1][:,:,:,0],axis=max_axis)#get_zstack_max(os.path.join(input_filepath,input_file),0,1)
    img_last = np.max(Z_dict[0][:,:,:,-1],axis=max_axis)#get_zstack_max(os.path.join(input_filepath,input_file),0,0)
    img2_last = np.max(Z_dict[1][:,:,:,-1],axis=max_axis)#get_zstack_max(os.path.join(input_filepath,input_file),0,1)      
 
    # normalization values (image cutoff value)
    ic = np.percentile(img,percent_cutoff)
    ic2 = np.percentile(img2,percent_cutoff)
    ic_last = np.percentile(img_last,percent_cutoff)
    ic2_last = np.percentile(img2_last,percent_cutoff)
    
    def get_normalization(frame):         
        """
        Returns the linearly interpolated normalization.
        """
        if adaptive_LUT == True:
            return ic - (ic-ic_last)*frame/sizeT , ic2 - (ic2-ic2_last)*frame/sizeT
        else:
            return ic, ic2
    
    ic, ic2 = get_normalization(0)
    img = img/ic
    img2 = img2/ic2
    img[img>1] = 1
    img2[img2>1] = 1
    img_array = np.zeros((img.shape[0],img.shape[1],3))
    img_array[:,:,0] = img2
    img_array[:,:,1] = img
    img_array[:,:,2] = img2
    im = plt.imshow(img_array, animated=True)
    txt = plt.text(20,30,f'Frame: {0}',color='w',)

    line1_dict = {}
    line0_dict = {}
    text0_dict = {}
    particle_set0 =  linked_df[(linked_df.channel==0)].particle.unique()
    particle_set1 =  linked_df[(linked_df.channel==1)].particle.unique()
    

    for p in particle_set0:
        line0_dict[p] = {}
        line0, = plt.plot([0],[0],color='cyan',alpha=line_alpha)
        line0_dict[p] = line0
        
        text0 = plt.text(0,0,f'{p}',color='white',) 
        text0_dict[p] = {}
        text0_dict[p] = text0
        text0_dict[p].set_visible(False)
    for p in particle_set1:
        line1_dict[p] = {}
        line1, = plt.plot([1],[1],color='orange',alpha=line_alpha)
        line1_dict[p] = line1


    rect = Rectangle((sizeX-bar_width-20,sizeY-20), bar_width, bar_height,color='w', angle=0.0)
    plt.gca().add_patch(rect)
    plt.text(sizeX-bar_width//2-20,sizeY-30,r'{} $\mu$m'.format(desired_bar_length_um),color='w',)

    fig.gca().get_xaxis().set_visible(False)
    fig.gca().get_yaxis().set_visible(False)
    
    if verbose ==True:
        print('starting...')
    
    # update the figure
    img_array = np.zeros((img.shape[0],img.shape[1],3))
    
    def updatefig_2colour(*args):
        img = np.max(Z_dict[0][:,:,:,args[0]],axis=max_axis)
        img2 = np.max(Z_dict[1][:,:,:,args[0]],axis=max_axis)
        frame = args[0]
        
        if adaptive_LUT==True:
            ic, ic2 = get_normalization(frame)
        else:
            ic, ic2 = get_normalization(0)

        img = img/ic
        img2 = img2/ic2
        img[img>1] = 1
        img2[img2>1] = 1
        
        img_array[:,:,0] = img2
        img_array[:,:,1] = img
        img_array[:,:,2] = img2
        im.set_array(img_array)
        txt.set_text(f'Frame: {frame}')

        for p in particle_set0:
            X = linked_df[(linked_df.particle==p) & (linked_df.channel==0) & linked_df.frame.between(frame-disp_line_len,frame, inclusive=True)].x.values
            Y = linked_df[(linked_df.particle==p) & (linked_df.channel==0) & linked_df.frame.between(frame-disp_line_len,frame, inclusive=True)].y.values
            line0 = line0_dict[p]
            line0.set_ydata(Y)
            line0.set_xdata(X)
            
            if len(X)>0:
                # position the number to not overlap with the trajectory
                # position the number such that it falls within the frame
                if adaptive_text==True:
                    mx = np.nanmedian(X)
                    my = np.nanmedian(Y)
                    dx = X[-1]-mx
                    dy = Y[-1]-my
                    norm = np.linalg.norm((dx,dy))
                    if norm == 0:
                        dx = text_dist/np.sqrt(2)
                        dy = text_dist/np.sqrt(2)
                    else:
                        dx = dx/norm*text_dist
                        dy = dy/norm*text_dist
                else:
                    dx = text_dist/np.sqrt(2)
                    dy = text_dist/np.sqrt(2)
                    
                    
                if X[-1]+dx <= sizeX-text_dist:
                    text0_dict[p].set_x(X[-1]+dx)
                else:
                    text0_dict[p].set_x(X[-1]-dx)
                    
                if Y[-1]+dy <= sizeY-text_dist:
                    text0_dict[p].set_y(Y[-1]+dy)
                else:
                    text0_dict[p].set_y(Y[-1]-dy)
                    
                text0_dict[p].set_visible(True)
            
        for p in particle_set1:
            X = linked_df[(linked_df.particle==p) & (linked_df.channel==1) & linked_df.frame.between(frame-disp_line_len,frame, inclusive=True)].x.values
            Y = linked_df[(linked_df.particle==p) & (linked_df.channel==1) & linked_df.frame.between(frame-disp_line_len,frame, inclusive=True)].y.values
            line1 = line1_dict[p]
            line1.set_ydata(Y)
            line1.set_xdata(X)
                        
        if verbose==True:
            print(args[0])
        return im,

    ani = animation.FuncAnimation(fig, updatefig_2colour,
                                  frames=np.arange(0,sizeT,1), 
                                  interval=millisecond_per_frame, blit=True)
    #plt.show() # this line was causing issues with snakemake
    ani.save(f'{output_filename}') 
    
    
    fig.clf()
    plt.close()
    gc.collect() 
    print("sleeping 3 seconds")
    time.sleep(3) # give some time
    
# internal helper function to help with plotting
def gridspec_inches(
    wcols,
    hrows,
    wspace=0.75,
    hspace=0.5,
    fig_kwargs={}):

    fig = plt.figure()
    fig_height_inches = (
        sum(hrows)
        )

    fig_width_inches = (
        sum(wcols)
        )

    fig=plt.figure(
        figsize=(fig_width_inches,fig_height_inches),
        subplotpars=matplotlib.figure.SubplotParams(
        left=0,
        right=1,
        bottom=0,
        top=1,
        wspace =0,
        hspace = 0.0),
        **fig_kwargs)
    fig.set_size_inches(fig_width_inches,fig_height_inches,forward=True)

    gs = matplotlib.gridspec.GridSpec(
        len(hrows),
        len(wcols),
        left=0,
        right=1,
        top=1,
        bottom=0,
        wspace=wspace,
        hspace=hspace,
        width_ratios=wcols,
        height_ratios=hrows
        )
    return fig, gs