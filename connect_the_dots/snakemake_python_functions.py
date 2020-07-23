import h5py
import os
import bioformats
import numpy as np
import javabridge
import re
from tifffile import xml2dict
from itertools import product
import pandas as pd

################################################################################
# HELPER FUNCTIONS
################################################################################
def make_dynamic_file_name(in_files,dynamic_str):
    """
    Makes a file "dynamic".

    in_files : list of str
        List of files and paths
    dynamic_str : str
        String appended to the file name in each of `in_files` files.
        String should contain snakemake wild-cards for use with dynamic().
        e.g. dynamic_str='-frame_{fr}_channel_{ch}_Object Predictions.h5'
        appends the schema name of subfolders in `in_files` files.

    """

    # make dynamic file names from hdf5 subfolders
    dynamic_file_name = []
    for file in in_files:
        dynamic_file_name.append(file+dynamic_str)
    return dynamic_file_name


def filenames_by_extension(in_files_list, out_files_path=None,ext='.h5'
                           ,append_filename_to_path=False):
    """
    Creates new file names based on the file names in `in_files_list`.
    New file names will be placed in out_files_path with a specified
    extension.

    Parameters
    ----------
    in_files_list : list of str
        List of file names (and the full path).
    out_files_path : list of str, or str
        File path(s) that will associated with the output file names.
        By default, the file will be placed in the same folder as
        paths in each of `in_files_list` in a subfolder with the file name.
    ext : 'str'
        File name extension
    append_filename_to_path : {T,F}, optional
        Appends the file

    Returns
    -------
    list
        List containing the names of files specified by the above rules.
    """

    # check if out_files is a list of the correct length
    if isinstance(out_files_path, list) == True:
        assert len(files_list) == len(metadata_list)
    elif out_files_path is None:
        out_files_path = [os.path.dirname(x) for x in in_files_list]
    else:
        out_files_path = [out_files_path]*len(in_files_list)


    intermediates = []
    for file, filepath in zip(in_files_list,out_files_path):
        base = os.path.splitext(os.path.basename(file))[0]
        base_ext = base+ext

        # optionally append file a name to the out_files_path
        if append_filename_to_path == True:
            full_name = os.path.join(filepath,base,base_ext)
        else:
            full_name = os.path.join(filepath,base_ext)
        intermediates.append(full_name)
    return intermediates


def czi_dump_zstack_info(in_filenames,out_filename='file_info.csv',
                            out_filepath=None):
    """
    Retrieves partial metadata from CZI files and dumps to csv. Data goes into
    a subsubdirectory with the `in_filenames` file name in the directory
    specified by `out_file_path`.

    Parameters
    ----------
    in_filenames : set
        A set containing the input CZI files (and their path) for which to
        retrieve metadata z-stacks and save into hdf5 format. `in_filenames`
        should contain the full path (and name) of the file.

    out_filename : string, optional
        Name of the output files.

    out_filepath : string, optional
        Path to the desired output files. By default, it is the same
        path as in `in_filenames`

    Returns
    -------
    csv file
        File containing metadata on number of frames, channels and z-stack size.
    """

    # start java virtual machine (used by bioformats to read czi files)
    javabridge.start_vm(class_path=bioformats.JARS)

    # iterate over the czi_image files provided
    for f in in_filenames:
        czi_image = f

        # get image metadata
        metadata = xml2dict(bioformats.get_omexml_metadata(czi_image))
        sizeT = metadata['OME']['Image']['Pixels']['SizeT']
        sizeX = metadata['OME']['Image']['Pixels']['SizeX']
        sizeY = metadata['OME']['Image']['Pixels']['SizeY']
        sizeZ = metadata['OME']['Image']['Pixels']['SizeZ']
        num_channels = len(metadata['OME']['Image']['Pixels']['Channel'])

        # create partial metadata DataFrame
        df = pd.DataFrame.from_dict({'sizeT':[sizeT], 'sizeX':[sizeX],
                            'sizeY':[sizeY], 'num_channels': [num_channels]})

        if out_filepath is None:
            out_filepath = os.path.dirname(in_filenames)

        # create save folder if it does not exist
        base = os.path.splitext(os.path.basename(czi_image))[0]
        save_dir_name = os.path.join(out_filepath, base)
        if not os.path.exists(save_dir_name):
            os.mkdir(save_dir_name)

        # save csv file with the metadata
        df.to_csv(os.path.join(save_dir_name, out_filename))

    # kill java virtual machine
    javabridge.kill_vm()
    return

def filenames_from_frames_channels_range(files_list, filepath, first_frame,
                                        last_frame, num_channels,ext='.h5'):
    """
    Creates a list of hdf5 files provided various input informations.

    Parameters
    ----------
    files_list : list of str
        List of file names.
    filepath : str
        Absolute path to the files.
    first_frame : int
        First frame for which to start creating file names.
    last_frame : int
        Last frame for which to end creating file names.
    num_channels : int
        Number of channels for which to create file names
    ext : str, optional
        Specifies the file extention to append to the created files. By
        default  it is '.h5'

    Returns
    -------
    list
        List containing the names of files specified by the above rules.
    """
    intermediates = []
    frame_range = np.arange(first_frame,last_frame)
    channel_range = range(num_channels)
    for (i,f,c) in product(files_list,frame_range,channel_range):
        fname = '{}/{}/frame_{}_channel_{}{}'.format(filepath,i,f,c,ext)
        intermediates.append(fname)
    return intermediates

def filenames_from_metadata(files_list, metadata_list,filepath=None,ext='.h5'):
    """
    Creates a list of hdf5 files provided various input informations. Assumes
    that files_list and metadata_list are properly ordered relative to each
    other.

    Parameters
    ----------
    files_list : list of str
        List of file names.
    metadata_list : list of str
        File path to the metadata associated with a file in `files_list`.
        Metadata must be a csv file containing the fields `sizeT` and
        `num_channels`.

    Returns
    -------
    list
        List containing the names of files specified by the above rules.
    """
    assert len(files_list) == len(metadata_list)

    for file, metadata in zip(files_list,metadata_list):
        df = pd.read_csv(metadata)
        last_frame = df['sizeT'].iloc[0]
        num_channels = df['num_channels'].iloc[0]
        intermediates = []
        frame_range = range(last_frame)
        channel_range = range(num_channels)

        if filepath is None:
            filepath = os.path.dirname(file)

        for (f,c) in product(frame_range,channel_range):
            fname = '{}/{}/frame_{}_channel_{}{}'.format(filepath,file,f,c,ext)
            intermediates.append(fname)
        return intermediates

def czi_to_hdf5(in_filenames, out_filenames, img_info=None,channels_list=None):
    """
    Retrieves z-stacks from CZI files and saves z-stacks to hdf5 format,
    with one z-stack per "frame" folder.

    Parameters
    ----------
    in_filenames : set
        A set containing the input CZI files (and their path) for which to
        remove z-stacks and save into hdf5 format.

    out_filenames : set
        Name of the output files (including the file path). Files in `out_files`
        must specify the `channel` and `frame`, and are formatted as
        "/some/path/frame_{}_channel_{}.h5".

    img_info : list of 5-tuples of ints, optional
        A list of tuples contains size-related metadata for each of the
        files in `in_filenames`. The following image values:
        (sizeX,sizeY,sizeZ,sizeT,num_channels).

    channels_list : list of ints, optional
        A list of integers specifying which channels of the czi image series
        to convert to hdf5 format. By default `channels_list` is None, in which
        case all channels are saved.

    Returns
    -------
    hdf5 file
        z-stacks from the CZI files saved into hdf5 format into folder `zstack`.
    """

    # start java virtual machine (used by bioformats to read czi files)
    javabridge.start_vm(class_path=bioformats.JARS)

    # iterate over the czi_image files provided
    for czi_image, h5_image in zip(in_filenames,out_filenames):

        # get image metadata
        if img_info is None:
            metadata = xml2dict(bioformats.get_omexml_metadata(czi_image))
            sizeT = metadata['OME']['Image']['Pixels']['SizeT']
            sizeX = metadata['OME']['Image']['Pixels']['SizeX']
            sizeY = metadata['OME']['Image']['Pixels']['SizeY']
            sizeZ = metadata['OME']['Image']['Pixels']['SizeZ']
            num_channels = len(metadata['OME']['Image']['Pixels']['Channel'])
        else:
            (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info

        # initialize a z-stack
        zstack = np.zeros((sizeZ,sizeY,sizeX))

        # make sure output file directory exists
        if not os.path.exists(os.path.dirname(h5_image)):
            os.makedirs(os.path.dirname(h5_image))

        # create h5 files, with subfolders for each z-stack and channel
        h5f = h5py.File(h5_image, 'w')
        reader = bioformats.ImageReader(czi_image)

        # retrieve z-stacks from the czi file
        #if channels_list is None:
        #    num
        channels = range(num_channels)

        for frame, channel in product(range(sizeT),channels):
            for z in range(sizeZ):
                zstack[z,:,:] = reader.read(t=frame,z=z,c=channel)
            # save z-stack to output hdf5 file
            out_filename = 'frame_{}_channel_{}'.format(frame,channel)
            h5f.create_dataset(out_filename, data=zstack,compression='gzip')
        reader.close()
        h5f.close()

    # kill java virtual machine
    javabridge.kill_vm()

def czi_zstack_to_hdf5(in_filenames, out_files, img_info=None):
    """
    Retrieves z-stacks from CZI files and saves z-stacks to hdf5 format.

    Parameters
    ----------
    in_filenames : set
        A set containing the input CZI files (and their path) for which to
        remove z-stacks and save into hdf5 format.

    out_files : set
        Name of the output files (including the file path). Files in `out_files`
        must specify the `channel` and `frame`, and are formatted as
        "/some/path/frame_{}_channel_{}.h5".

    img_info : list of 5-tuples of ints, optional
        A list of tuples contains size-related metadata for each of the
        files in `in_filenames`. The following image values:
        (sizeX,sizeY,sizeZ,sizeT,num_channels).

    Returns
    -------
    hdf5 file
        z-stacks from the CZI files saved into hdf5 format into folder `zstack`.
    """

    # start java virtual machine (used by bioformats to read czi files)
    javabridge.start_vm(class_path=bioformats.JARS)

    # iterate over the czi_image files provided
    while in_filenames:
        czi_image = in_filenames.pop()

        # get image metadata
        if img_info is None:
            metadata = xml2dict(bioformats.get_omexml_metadata(czi_image))
            sizeT = metadata['OME']['Image']['Pixels']['SizeT']
            sizeX = metadata['OME']['Image']['Pixels']['SizeX']
            sizeY = metadata['OME']['Image']['Pixels']['SizeY']
            sizeZ = metadata['OME']['Image']['Pixels']['SizeZ']
            num_channels = len(metadata['OME']['Image']['Pixels']['Channel'])
        else:
            (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info

        zstack = np.zeros((sizeZ,sizeY,sizeX))

        # determine what zstacks to pull out based on `out_files` names
        basename = os.path.basename(czi_image)
        out_files_partial = [x for x in out_files if basename in x]

        for out_file in out_files_partial:
            frame_re = re.search('Frame_(\d+)', out_file, re.IGNORECASE)
            channel_re = re.search('channel_(\d+)', out_file, re.IGNORECASE)
            frame = frame_re.group(1)
            channel = channel_re.group(1)
            out_filepath = os.path.dirname(out_file)

            if not os.path.exists(out_filepath):
                os.makedirs(out_filepath)

            # get zstack
            with bioformats.ImageReader(czi_image) as reader:

                for z in range(sizeZ):
                    zstack[z,:,:] = reader.read(t=frame,z=z,c=channel)

            # save z-stack to output hdf5 file
            out_filename = 'frame_{}_channel_{}.h5'.format(frame,channel)
            with h5py.File(os.path.join(out_filepath,out_filename), 'w') as h5f:
                h5f.create_dataset('zstack', data=zstack)
                h5f.close()

    # kill java virtual machine
    javabridge.kill_vm()


def hdf5_zstack_to_max_projection(in_filenames,
                                out_file_suffix='_max_projection.h5',
                                axis=0):
    """
    Retrieves z-stacks from hdf5 files and saves a maximum projection image
    hdf5 format.

    Parameters
    ----------
    in_filenames : set
        A set containing the input hdf5 files for which to retrieve z-stacks
        and calculate maximum projection.

    out_file_suffix : set
        Specifies the extension and suffix of the output files.

    axis : int
        Defines the axes for which to do the maximum projection. By default,
        it is set to 0 (which is nominally the "Z" dimension).

    Returns
    -------
    hdf5 file
        Maximum projection of z-stacks from the input hdf5 files.
    """

    for file in in_filenames:

        # make output file name
        filepath = os.path.dirname(file)
        basename, extension = os.path.splitext(os.path.basename(file))
        output_file = os.path.join(filepath,basename+out_file_suffix)

        with h5py.File(file,'r') as h, h5py.File(output_file, 'w') as hh:
            # for each dataset in the hdf5 file
            for key in h.keys():
                zstack = np.array(h[key])
                hh.create_dataset(key, data=np.max(zstack,axis=axis),compression='gzip')
