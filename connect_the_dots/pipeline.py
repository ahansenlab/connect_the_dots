
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
        filepath_list.append(os.path.join(start_folder,path.parent.name))
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
    
    

def natural_sort(l):
    # found here: https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def wavelet_filter_zstack(zstack,filtered_levels=[0,4],wavelet_type='haar'):
    """
    Perform 3D (or 4D) discrete wavelet transform and coefficient filtering.
    
    Parameters
    ----------
    zstack : ndarray
        Numpy ndarray containing 3 (or 4) dimensions. Dimensions should all be even for best performance!        
    
    filtered_levels : list of ints
        Wavelet coefficient levels to be removed for the image reconstruction.
        
    Returns
    -------
    rec : ndarray 
        Filtered zstack
        
    num_levels : int
        Number of levels that could be accessed
    """

    ## Do Wavelet-decomposition of the image
    levels = filtered_levels 
    coeffs = pywt.wavedecn(zstack,wavelet_type)

    ## Reconstruct image, after filtering different wavelet coefficients 
    for level in levels:
        if level == 0:
            coeffs[0] = np.zeros_like(coeffs[0])
        else:
            try:
                coeffs[level] = {k: np.zeros_like(v) for k, v in coeffs[level].items()} 
            except:
                warnings.warn(f"Wavelet level {level} does not exist. There are only {len(levels)} levels." + \
                              "Check image input size. max_levels = log2(smallest_dimension_of_data/(filter_size-1))")
    rec = pywt.waverecn(coeffs, wavelet_type) 

    ## crop rec to match zstack 
    if (rec.shape != zstack.shape):
        if len(zstack.shape) == 3:
            rec = rec[0:zstack.shape[0],0:zstack.shape[1],0:zstack.shape[2]]
        elif len(zstack.shape) == 4:
            rec = rec[0:zstack.shape[0],0:zstack.shape[1],0:zstack.shape[2],0:zstack.shape[3]]
        elif len(zstack.shape) == 2:
            rec = rec[0:zstack.shape[0],0:zstack.shape[1]]
            
    return rec

def filter_zstack_DoG(zstack,dog_sigma1 = 1.5,dog_sigma2 = 15,absolute_value=True):
    """
    Applies Difference of Gaussian (DoG) filtering on a single z-stack. 

    Parameters
    ----------
    zstack : numpy.ndarray [sizeZ by sizeY by sizeX]
             Z-stack of the image series for a single channel (containing 3 spatial dimentions)
    
    dog_sigma1 : float, optional
                 Standard deviation of the first Gaussian distribution of the DoG filter. 
                 `dog_sigma1` should be close in size to the "dots" being tracked.
             
    dog_sigma2 : float, optional
                 Standard deviation of the second Gaussian distribution of the DoG filter. 
                 `dog_sigma2` should be ~10-times larger than `dog_sigma1`; it helps to smooth 
                 local sources noise and background of the image. 
    
    absolute_value : {T,F}, optional
                     Toggles on/off taking the absolute value of the DoG filter result. 
                     
    Returns
    -------
    filtered_zstack : numpy.ndarray
                      Absolute value of Difference of Gaussian filtered z-stack. 
    """
    filtered_zstack = gaussian_filter(zstack,dog_sigma1)- gaussian_filter(zstack,dog_sigma2)
    
    if absolute_value==True:
        filtered_zstack = np.abs(filtered_zstack)
    
    return filtered_zstack


def filter_zstack_LoG(img,log_sigma=2,num_std_threshold=0):
    """
    Returns the Laplacian of Gaussian (LoG) filtered image and the image of 
    edges detected by the LoG method applied on the image.
    
    Parameters
    ----------
    img : numpy.ndarray
            Image on which to perform the edge detection.
            
    log_sigma : float
            The standard deviation of the gaussian used in 'Laplacian of Gaussian'. It is
            akin to a smoothing parameter.
            
    num_std_threshold : float
            The number of standard deviations used above "zero" to classify a zero-crossing
            event of the 'Laplacian of Gaussian' as being an edge. Only set `num_std_threshold`
            greater than zero; also, only 
    
    Returns
    -------
    log_img
    
    edges_img
    
    """
    log_img = gaussian_laplace(img, log_sigma)
    threshold = np.absolute(log_img).std() * num_std_threshold
    edges_img = np.zeros(img.shape)
    w = edges_img.shape[1]
    h = edges_img.shape[0]
    for y in range(1, h - 1):
        for x in range(1, w - 1):
            region = log_img[y-1:y+2, x-1:x+2]
            val = log_img[y, x]
            max_val = region.max()
            min_val = region.min()
            if (val > 0):
                zerocross = True if min_val < 0 else False
            else:
                zerocross = True if max_val > 0 else False
            if ((max_val - min_val) > threshold) and zerocross:
                edges_img[y, x] = 1
    return log_img, edges_img 

    
def get_image_objects(zstack_to_mask,zstack,min_size_pixels=20,
                      max_size_pixels=10000,
                      percentile_threshold=99.9,
                      do_dilation=True,
                      return_mask=False):
    """
    Creates a binary mask on the zstack, and segments the binary image. 3D objects are identified
    and selected based on their pixel sizes.
    
    Parameters
    ----------
    zstack_to_mask : numpy.ndarray [sizeZ by sizeY by sizeX]
        zstack used to create the binary mask (possibly filtered)
        
    zstack : numpy.ndarray [sizeZ by sizeY by sizeX]
        zstack used to obtain the intensity values of the dot localization (raw data)
        
    min_size_pixels : int
        Minimum size of object to be extracted by the segmentation
        
    max_size_pixels : int
        Maximum size of object to be extracted by the segmentation
    
    percentile_threshold : float (0 to 100)
        Percentile threshold with which to binarize the zstack and create object masks.
        
    do_dilation : {T,F}
        Perform morphological filtering and size selection on image after binary dilation. Default == True.
        
    return_mask : {T,F}
        Returns the mask. By default, the mask is not returned.
        
    Returns
    -------
    A data frame containing the following column headers:    
    x, y, z, mean_intensity, max_intensity, dot_size_pixels
    
    If return_mask == True:
    blobs : np.ndarray
        The thresholded zstack
    blobs_filt : np.ndarray
        The thresholded zstack, after feature size selection
    
    """
    
    ## Apply threshold to binarize image 
    stack_threshold = percentile_threshold
    blobs = zstack_to_mask.copy()
    blobs[blobs < np.percentile(blobs,stack_threshold)] = 0
    blobs[blobs>0] = 1

    
    if do_dilation == True:
        blobs = dilation(blobs)
        
    
    ## Measure blobs/region properties 
    blobs_labels = measure.label(blobs, background=0) 
    blob_metrics = measure.regionprops(blobs_labels, zstack) ## note that we use the unfiltered zstack as input
    
    # print('labels created')
    # get object centroids i.e. (z,y,x) coordinates, volumes and intensities
    centroids = [tuple(np.array(x.weighted_centroid,dtype=float)) for x in blob_metrics]
    vol = [x.area for x in blob_metrics]
    max_intensity = [x.max_intensity for x in blob_metrics]
    mean_intensity = [x.mean_intensity for x in blob_metrics]

    ## Filter localized dots based on size
    good_dots = [vi for vi, v in enumerate(vol) if ((v >= min_size_pixels) and (v <= max_size_pixels))]

    ## Update the list of localized dots
    X = [centroids[i][2] for i in  good_dots]
    Y = [centroids[i][1] for i in  good_dots]
    Z = [centroids[i][0] for i in  good_dots]

    mean_intensity = [mean_intensity[i] for i in  good_dots]
    max_intensity = [max_intensity[i] for i in  good_dots]
    dot_size_pixels = [vol[i] for i in  good_dots]

    
    df = pd.DataFrame({'x':X, 'y':Y,'z':Z,'mean_intensity':mean_intensity,
                  'max_intensity':max_intensity,'dot_size_in_pixels':dot_size_pixels})
    
    if return_mask == False:
        return df
    
#     blobs_filt = np.zeros_like(blobs)
#     for v in good_dots:
#         blobs_filt[blobs_labels==(v+1)] = 1
    
    return df, blobs#, blobs_filt


def get_localizations_iterative(filtered_zstack_timeseries, 
                                zstack_timeseries,
                                frames, 
                                channel,
                                percentile_threshold=99.95,
                                max_dot_size = 10000,
                                min_dot_size = 20, 
                                search_range=(5,15,15),
                                min_track_length = 15,
                                max_iterations=5,
                                min_dot_size_increment=5,
                                percentile_threshold_increment=0.01,
                                current_iteration=0,
                                verbose=True):
    """
    Use trackpy to create trajectories from dot localizations.
    
    Parameters
    ----------
    filtered_zstack_timeseries : numpy.ndarray
        Post-filtering 4-D zstack of the image series data used to make a binary mask.
        
    zstack_timeseries : numpy.ndarray
        Raw 4-D zstack of the image series data.
    
    frames : range or list of int
        Frames to process.
        
    channel : int
        
    percentile_threshold : float 
        Percentile from 0 to 100 by which to threshold the `filtered_zstack_timeseries`.
        
    max_dot_size : int
        Maximum size of the dots in numbers of pixels
        
    min_dot_size : int
        Minimum size of the dots in numbers of pixels
    
    search_range : 3-tuple of ints
        Maximum frame-to-frame displacement (Z,Y,X) of dots
        
    memory : int
        Number of frames a dot may "disappear" and still be called the same trajectory.
        
    min_track_length : int
        Minimum length of trajectories kept
        
    max_iterations : int
        Maximum number of iterations of parameter updates for trajectory creation
        
    min_dot_size_increment : int
        Increment added to `min_dot_size` if trajectory creation fails.
        
    percentile_threshold_increment : float
        Increment added to `percentile_threshold` if trajectory creation fails by `dot_size_increment`.
        
    current_iteration : int
        Keeps track of the number of iterations. 
    
    Returns
    -------
    dict : dictionary of particle localizations
    dict : dictionary of particle trajectories
    """
    
    if verbose == True:
        print(f"Starting iteration: {current_iteration}")
                  
    loc_df = []
    for frame in frames:
        if verbose==True:
            print(f"Getting objects channel: {ch}, frame: {frame}")

        df = get_image_objects(filtered_zstack_timeseries[:,:,:,frame],
                                                 zstack_timeseries[:,:,:,frame],
                                                 percentile_threshold=percentile_threshold,
                                                 min_size_pixels=min_dot_size,
                                                 max_size_pixels=max_dot_size)
        df['frame'] = [frame]*len(df)
        df['channel'] = [channel]*len(df)
        loc_df.append(df)
    loc_df = pd.concat(loc_df,ignore_index=True)

    if verbose==True:
        print(f"Linking trajectories, channel: {ch}")

    linked_df = None
    count = 0
    while linked_df is None:
        try:
            linked_df = create_trajectories(loc_df, 
                        max_dot_size = max_dot_size,
                        min_dot_size = min_dot_size+min_dot_size_increment*count, 
                        search_range=(5,15,15),
                        min_track_length = 15)
        except:
            if verbose==True:
                print(f"Trajectory linking failed. Incrementing the minimum dot size.")

            # increase the threshold on the minimum dot size
            count += 1
            linked_df = None
            
    # if increasing the minimum dot size threshold does not help
    if linked_df is None:

        if verbose==True:
            print(f"Trajectory linking failed. Re-processing image series with higher percentile threshold.")

        # increase the percentile threshold
        percentile_threshold += percentile_threshold_increment
        loc_df, linked_df = get_localizations_iterative(filtered_zstack_timeseries, 
                                        zstack_timeseries,
                                        frames, 
                                        channels,
                                        percentile_threshold=percentile_threshold,
                                        max_dot_size = max_dot_size,
                                        min_dot_size = min_dot_size, 
                                        search_range=[(5,15,15)],
                                        min_track_length = 15,
                                        max_iterations=5,
                                        current_iteration=current_iteration+1)        

        warnings.warn("One iteration failed")
    return loc_df, linked_df

def create_trajectories(loc_df, 
                        max_dot_size = np.inf,
                        min_dot_size = 0, 
                        search_range=(5,15,15),
                        memory=5,
                        min_track_length = 15):
    """
    Use trackpy to create trajectories from dot localizations.
    
    Parameters
    ----------
    loc_df : Pandas DataFrame
        DataFrame of object localizations in 3D and time. 
    
    max_dot_size : int
        Maximum size of the dots in numbers of pixels
        
    min_dot_size : int
        Minimum size of the dots in numbers of pixels
        
    search_range : 3-tuple of ints
        Maximum frame-to-frame displacement (Z,Y,X) of dots
        
    memory : int
        Number of frames a dot may "disappear" and still be called the same trajectory.
        
    min_track_length : int
        Minimum length of trajectories kept
    
    Returns
    -------
    dict : dictionary of particle trajectories
    """
    
    # filter the trajectories based on dot size
    tmp_df = loc_df.copy()
    tmp_df = tmp_df[(tmp_df.dot_size_in_pixels<max_dot_size) & ((tmp_df.dot_size_in_pixels>min_dot_size))]
    
    # create trajectories
    tmp_df = trackpy.link_df(tmp_df, \
            search_range=tuple(search_range), \
            memory=memory)
    
    # filter trajectories based on minimum length
    tmp_df = tmp_df.groupby('particle').filter(lambda x: len(x) > min_track_length)

    return tmp_df


def link_trajectories_across_channels(linked_df, 
                                      min_overlap_length=10,
                                     corrcoeff_min =0.3,
                                     distance_max = 35):
    """
    Links trajectories together across channels
    
    Parameters
    ----------
    linked_df : pandas.DataFrame
        DataFrame of particle trajectories (connected separately for each channel) in 3D and time.         
    
    min_overlap_length : int
        Minimum number of frames of overlap between channels necessary to connect two trajectories
        
    corrcoeff_min : float
        Minimum correlation coefficient between dot pairs (for each dimension separately) in order
        for trajectories to be linked together.
    
    distance_max : int
        Maximum distance in pixels between dots in connected trajectories
        
    Returns
    -------
    joined_df : pandas.DataFrame
        DataFrame of particle trajectories (connected separately for each channel) in 3D and time.  
        This DataFrame may contain some trajectories that do not meet the likely pair criteria. To 
        obtain a clean list, use `linked_trajectory_ids`.
        
    linked_trajectory_ids : pandas.DataFrame
        List of trajectory IDs in `joined_df` meeting the criteria for likely pairs.
    
    
    """
    
    
    fake_linked_df = linked_df.copy()

    # get 4D tracks for Channel 0
    particles0 = list(set(fake_linked_df[fake_linked_df.channel==0].particle))
    tracks_df0 = {}
    for p in particles0:
        X = fake_linked_df[fake_linked_df.particle==p][['frame','x','y','z']].values
        tracks_df0[p] = X

    # get 4D tracks for Channel 1
    particles1 = list(set(fake_linked_df[fake_linked_df.channel==1].particle))
    tracks_df1 = {}
    for p in particles1:
        X = fake_linked_df[fake_linked_df.particle==p][['frame','x','y','z']].values
        tracks_df1[p] = X

    min_id_value = np.min(particles1)
    
    ## Find likely pairs of trajectories
    likely_pairs = []
    for p0i in tracks_df0:
        p0 = tracks_df0[p0i]
        for p1i in tracks_df1:
            p1 = tracks_df1[p1i]
            times, comm0, comm1 = np.intersect1d(p0[:,0],p1[:,0],return_indices=True)
            if len(times) > min_overlap_length:  
                corrcoeffX = np.corrcoef(p0[comm0,1],p1[comm1,1])[0,1]
                corrcoeffY = np.corrcoef(p0[comm0,2],p1[comm1,2])[0,1]
                corrcoeffZ = np.corrcoef(p0[comm0,3],p1[comm1,3])[0,1]
#                 residualsX = np.mean(np.abs(p0[comm0,1]-p1[comm1,1]))
#                 residualsY = np.mean(np.abs(p0[comm0,2]-p1[comm1,2]))
#                 residualsZ = np.mean(np.abs(p0[comm0,3]-p1[comm1,3]))
                residualsX = (np.abs(p0[comm0,1]-p1[comm1,1]))
                residualsY = (np.abs(p0[comm0,2]-p1[comm1,2]))
                residualsZ = (np.abs(p0[comm0,3]-p1[comm1,3]))
                distance = np.sqrt(residualsX**2+residualsY**2+residualsZ**2)

                # link trajectories if they satisfy a minimum correlation coefficient 
                # and are within an expected maximum distance (in nm) of residuals_max
                corrcoeff_good = all([corrcoeff_min < c for c in [corrcoeffX,corrcoeffY,corrcoeffZ]])
                residuals_good = all(distance < distance_max)
                if corrcoeff_good and residuals_good:
                    likely_pairs.append((int(p0i),int(p1i)))

    ################## REFINE THE PAIRING###################
    # cluster trajectories together 
    # find conflicts or overlaps

    clusters = find_clusters(likely_pairs)        
    likelier_pairs = []
    # generate all possible orderings within the cluster
    # for each ordering:
    #     check for overlaps within a channel 
    #      -- if overlap exists, discard this; continue
    #     check for continuous overlap across channels for all segments
    #       -- if overlap is non-contiguous discard; continue
    #     if overlaps are contiguous in time and non-overlapping within a channel
    #       -- score the total length of contig, score average correlation, score average distance; store
    # for each stored/scored ordering
    for ci, c in enumerate(clusters):
        if len(c) < 2:
            continue
        elif len(c)==2:
            likelier_pairs.append(tuple(sorted(c)))
        else:        
            # generate all possible orderings within the cluster
            orderings = partition_select(c)

            # for each ordering
            good_ordering_metrics = {}
            for si, seq in enumerate(orderings):
                # assert the ordering is at least of size 2
                if len(seq) <= 1:
                    #bad_orderings.append(si)
                    continue

                # split the channels and get time points for each group
                group0 = [tracks_df0[s][:,0] for s in seq if s < min_id_value]
                group1 = [tracks_df1[s][:,0] for s in seq if s >= min_id_value]


                # allow for an overlap equal to the number of groups collected -1
                max_overlap0 = 0#len(group0)-1
                max_overlap1 = 0#len(group1)-1

                # CHECK FOR OVERLAPS WITHIN A CHANNEL
                if len(group0) == 0:
                    continue
                elif len(group0)> 1:
                    has_overlap0 = reduce(np.intersect1d, group0)
                else:
                    has_overlap0 = []
                if len(group1) == 0:
                    continue
                elif len(group1)> 1:
                    has_overlap1 = reduce(np.intersect1d, group1)
                else:
                    has_overlap1 = []
                if (len(has_overlap0)>max_overlap0 or len(has_overlap1)>max_overlap1)==True:
                    #bad_orderings.append(si)
                    continue

                # CHECK FOR FULLY CONNECTED TRAJECTORIES
                # count the number of connections for each trajectory
                num_links = {f"group0_{k}":0 for k in range(len(group0))}
                num_links.update({f"group1_{k}":0 for k in range(len(group1))})
                for gi0, g0 in enumerate(group0):
                    for gi1, g1 in enumerate(group1):
                        intersect = np.intersect1d(g0,g1)
                        if len(intersect)> 0:
                            num_links[f"group0_{gi0}"] += 1
                            num_links[f"group1_{gi1}"] += 1
                # if any elements are disconnected, remove this ordering 
                for key in num_links:
                    if num_links[key] == 0:
                        #bad_orderings.append(si)
                        continue

                # WITH THE REMAINING "GOOD ORDERINGS", compute various metrics
                # combine all group0 together
                # combine all group1 together
                times0 = group_arrays(group0)
                times1 = group_arrays(group1)
                times, comm0, comm1 = np.intersect1d(times0,times1,return_indices=True)


                # get X, Y, Z positions for the trajectory
                groupX0 = group_arrays([tracks_df0[s][:,1] for s in seq if s < min_id_value])
                groupX1 = group_arrays([tracks_df1[s][:,1] for s in seq if s >= min_id_value])
                groupY0 = group_arrays([tracks_df0[s][:,2] for s in seq if s < min_id_value])
                groupY1 = group_arrays([tracks_df1[s][:,2] for s in seq if s >= min_id_value])
                groupZ0 = group_arrays([tracks_df0[s][:,3] for s in seq if s < min_id_value])
                groupZ1 = group_arrays([tracks_df1[s][:,3] for s in seq if s >= min_id_value])

                # get only the overlapping subset of trajectory values
                groupX0 = groupX0[comm0]
                groupX1 = groupX1[comm1]
                groupY0 = groupY0[comm0]
                groupY1 = groupY1[comm1]
                groupZ0 = groupZ0[comm0]
                groupZ1 = groupZ1[comm1]

                # get mean correlation coefficient 
                corrcoeff = np.corrcoef(groupX0,groupX0)[0,1]/3 + \
                            np.corrcoef(groupY0,groupY0)[0,1]/3 + \
                            np.corrcoef(groupZ0,groupZ0)[0,1]/3

                # get the mean distance 
                residuals = np.nanmedian(np.sqrt((groupX0-groupX1)**2 + \
                                            (groupY0-groupY1)**2 + \
                                            (groupZ0-groupZ1)**2))

                #print((ci, si))
                # get the overlap length
                good_ordering_metrics[si] = (len(times),residuals,corrcoeff)

            # sort by longest trajectory first
            best_lengths = sorted(good_ordering_metrics, key=lambda k: good_ordering_metrics[k][0],reverse=True)

#             best_trajectory_idx = np.argmin(best_lengths)
            # sort by hightest correlation first
            best_correlations = sorted(good_ordering_metrics, key=lambda k: good_ordering_metrics[k][2],reverse=True)

            # sort by lowest distance first
            best_distance = sorted(good_ordering_metrics, key=lambda k: good_ordering_metrics[k][1],reverse=False)

            # ranked choice vote
            rank_choice =  {k : 0 for k in best_lengths}
            for rank, bl, bc, bd in zip(range(len(best_lengths)),best_lengths,best_correlations,best_distance):
                rank_choice[bl] += rank
                #rank_choice[bc] += rank  # this metric is useless here -- most correlations are close to 1
                #rank_choice[bd] += rank

            # get best ranked trajectory
            min_rank = min(rank_choice.values()) 
            best_trajectory_idx = [key for key in rank_choice if rank_choice[key] == min_rank] 
            best_trajectory_idx = best_trajectory_idx[0]

            # note, these may not be *pairs*, but a tuples of some length 
            likelier_pairs.append(orderings[best_trajectory_idx]) 

    # update likely_pairs with likelier pairs
    likely_pairs = []
    for lp in likelier_pairs:
        lp = sorted(lp)
        if len(lp)==2:
            likely_pairs.append(tuple(lp))
        else:

            #print(lp)
            for lpi in lp[1:]:
                likely_pairs.append((lp[0],lpi) )
    # update particle names
    for p1i, p2i in likely_pairs:
        fake_linked_df['particle'].replace(p2i,p1i,inplace=True)

    linked_particles_ids = list(set([x[0] for x in likely_pairs]))
    return fake_linked_df, linked_particles_ids
    
        
        
def find_clusters( trajectory_pairs ):
    """
    Identifies clusters of trajectories which could be stitched together.
    
    Parameters
    ----------
    trajectory_pairs : list of 2-tuples of ints
        List of 2-tuples of particle IDs, containing potential matched pairs of
        trajectories across channels. This is used for stitching trajectories together.
        
    
    Returns
    -------
    clustered_tuples : n-tuples of ints
        List of n-tuples of particle IDs containing groups of potential trajectory segments.
    
    """
    tuples = trajectory_pairs
    # https://stackoverflow.com/questions/21646703/grouping-elements-from-2-tuples-recursively
    clusterlist=[]
    # clustermap maps an element to the id of the containing
    # cluster within clusterlist
    clustermap = {}

    # here we find the current cluster id for elem, by following the
    # chain within clusterlist, and replace that entire chain
    # with the new cluster id n.   We return the old cluster id.
    def set_cluster_id( elem, n ):
        if elem not in clustermap:
            return None
        k = clustermap[elem]
        # clusters may be replaced by references to other clusters,
        # we follow that chain
        while k < n and isinstance( clusterlist[k], int ):
            k1 = clusterlist[k]
            # this is optional, we make the chain shorter
            # by making this entry point directly to the current cluster
            clusterlist[k] = n
            k = k1
        return k

    for t in tuples:
        # for each tuple we create a new cluster
        thiscluster = set(t)
        n = len( clusterlist ) # the id of thiscluster
        for x in t:
            # we absorb existing clusters into the new one
            # if there is overlap
            k = set_cluster_id(x, n)
            if k is not None and k != n:
                thiscluster.update( clusterlist[k] )
                # we replace the existing cluster
                # with a reference to the new one
                clusterlist[k] = n 
            clustermap[x] = n
        clusterlist.append(thiscluster)

    return [ tuple(x) for x in clusterlist if isinstance( x, set ) ]  
  

def partition_select(vals,min_sz=2):
    """
    This function generates a unique list of all possible partitions
    between the items in `vals`. It also removes all partitions 
    smaller than min_sz in length.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    
    """
    vals = list(vals)
    
    def _partition(collection):
        #https://stackoverflow.com/questions/19368375/set-partitions-in-python
        if len(collection) == 1:
            yield [ collection ]
            return
        first = collection[0]
        for smaller in _partition(collection[1:]):
            # insert `first` in each of the subpartition's subsets
            for n, subset in enumerate(smaller):
                yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
            # put `first` in its own subset 
            #yield [ [ first ] ] + smaller
            yield [[first]] + smaller
         
    ordering = list(_partition(vals))
    flat_list = []
    for sublist in ordering:
        for item in sublist:
            flat_list.append(tuple(item))
    flat_list = list(set(flat_list))
    return [fl for fl in flat_list if len(fl)>= min_sz]

def group_arrays(list_of_arrays):
    return np.array(reduce(lambda x,y: list(x)+list(y),list_of_arrays))


def filter_large_displacements_in_trackpy_trajectories(linked_df,distance_max=35):
    """
    Trajectory quality check. Removes localizations with big frame-to-frame steps.
    Removal is based on absolute dot positions (i.e. within each channel) and not between channels.
    
    
    Parameters
    ----------
    linked_df : pandas.DataFrame
        DataFrame of dot trajectories and localizations.    
    
    Returns
    -------
    new_df_list : pandas.DataFrame
        A cleaned up version of the `linked_df`, without the large displacements.
    """
    
    new_df_list = []
    channels = linked_df.channel.unique()
    for ch in channels:
        for p in linked_df.particle.unique():
            times = linked_df.query(f"(particle=={p}) and (channel=={ch})").frame.values
            x = linked_df.query(f"(particle=={p}) and (channel=={ch})").x.values
            y = linked_df.query(f"(particle=={p}) and (channel=={ch})").y.values
            z = linked_df.query(f"(particle=={p}) and (channel=={ch})").z.values
            mean_intensity = linked_df.query(f"(particle=={p}) and (channel=={ch})").mean_intensity.values
            max_intensity = linked_df.query(f"(particle=={p}) and (channel=={ch})").max_intensity.values
            dot_size_in_pixels = linked_df.query(f"(particle=={p}) and (channel=={ch})").dot_size_in_pixels.values
            particle = linked_df.query(f"(particle=={p}) and (channel=={ch})").particle.values
            channel = linked_df.query(f"(particle=={p}) and (channel=={ch})").channel.values    

            bad_frames = np.r_[0,np.abs(np.diff(x))< distance_max ] == 1 
            if len(bad_frames) - sum(bad_frames)- 1 > 0:
                # find change points
                changepoints = np.r_[np.where(bad_frames==False)[0],np.inf]

                ind0 = []
                ind1 = []

                ind = 0
                count = 0
                for frame in range(len(x)):
                    if frame < changepoints[count]:
                        if  ind ==0:
                            ind0.append(frame)
                        else:
                            ind1.append(frame)

                    if frame == changepoints[count]:
                        count+=1
                        if  ind ==0:
                            ind1.append(frame)
                            ind = 1
                        else:
                            ind0.append(frame)
                            ind = 0
                ind_keep = np.argmax([len(ind0),len(ind1)])
                if ind_keep == 0:
                    x = [x[pos] for pos in ind0]
                    y = [y[pos] for pos in ind0]
                    z = [z[pos] for pos in ind0]
                    mean_intensity = [mean_intensity[pos] for pos in ind0] 
                    max_intensity = [max_intensity[pos] for pos in ind0]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind0]
                    particle = [particle[pos] for pos in ind0]
                    channel = [channel[pos] for pos in ind0]
                    times = [times[pos] for pos in ind0]
                else:
                    x = [x[pos] for pos in ind1]
                    y = [y[pos] for pos in ind1]
                    z = [z[pos] for pos in ind1]
                    mean_intensity = [mean_intensity[pos] for pos in ind1]
                    max_intensity = [max_intensity[pos] for pos in ind1]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind1]
                    particle = [particle[pos] for pos in ind1]
                    channel = [channel[pos] for pos in ind1]
                    times = [times[pos] for pos in ind1]

            ## filter out big displacements in Y
            bad_frames = np.r_[0,np.abs(np.diff(y))< distance_max ] == 1 
            if len(bad_frames) - sum(bad_frames)- 1 > 0:
                # find change points
                changepoints = np.r_[np.where(bad_frames==False)[0],np.inf]

                ind0 = []
                ind1 = []

                ind = 0
                count = 0
                for frame in range(len(x)):
                    if frame < changepoints[count]:
                        if  ind ==0:
                            ind0.append(frame)
                        else:
                            ind1.append(frame)

                    if frame == changepoints[count]:
                        count+=1
                        if  ind ==0:
                            ind1.append(frame)
                            ind = 1
                        else:
                            ind0.append(frame)
                            ind = 0
                ind_keep = np.argmax([len(ind0),len(ind1)])
                if ind_keep == 0:
                    x = [x[pos] for pos in ind0]
                    y = [y[pos] for pos in ind0]
                    z = [z[pos] for pos in ind0]
                    mean_intensity = [mean_intensity[pos] for pos in ind0] 
                    max_intensity = [max_intensity[pos] for pos in ind0]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind0]
                    particle = [particle[pos] for pos in ind0]
                    channel = [channel[pos] for pos in ind0]
                    times = [times[pos] for pos in ind0]
                else:
                    x = [x[pos] for pos in ind1]
                    y = [y[pos] for pos in ind1]
                    z = [z[pos] for pos in ind1]
                    mean_intensity = [mean_intensity[pos] for pos in ind1]
                    max_intensity = [max_intensity[pos] for pos in ind1]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind1]
                    particle = [particle[pos] for pos in ind1]
                    channel = [channel[pos] for pos in ind1]
                    times = [times[pos] for pos in ind1]

            new_df = pd.DataFrame({'particle':particle,
                          'channel':channel,
                          'frame':times,
                          'x':x,
                          'y':y,
                          'z':z,
                          'mean_intensity':mean_intensity,
                          'max_intensity':max_intensity,
                          'dot_size_in_pixels':dot_size_in_pixels})

            new_df_list.append(new_df)
        
    return pd.concat(new_df_list)


def infer_missing_dot_locations(gapped_df,max_gap_length = 10):
    """
    Takes existing joint trajectories dataframes and interpolates to guess where 
    dots may be in frames where localizations are missing.
    
    Parameters
    ----------
    gapped_df : pandas.DataFrame
       DataFrame of trajectories that have been linked across channels.
    
    max_gap_length : int
        Maximum number of frames with which to use within channel interpolation. If 
        the number of missing frames exceeds this value, displacements from the 
        adjacent channel are used.
        
    Returns
    -------
    search_df : pandas.DataFrame
        DataFrame containing best guesses of where a dot is likely to be in a missing frame.
    """
    
    particles = list(set(gapped_df.particle))
    search_df_list = []
    for p in particles:

        x0 = gapped_df.query(f'channel==0 & particle=={p}').x.values
        y0 = gapped_df.query(f'channel==0 & particle=={p}').y.values
        z0 = gapped_df.query(f'channel==0 & particle=={p}').z.values
        t0 = gapped_df.query(f'channel==0 & particle=={p}').frame.values
        i0 = gapped_df.query(f'channel==1 & particle=={p}').mean_intensity.values

        x1 = gapped_df.query(f'channel==1 & particle=={p}').x.values
        y1 = gapped_df.query(f'channel==1 & particle=={p}').y.values
        z1 = gapped_df.query(f'channel==1 & particle=={p}').z.values
        t1 = gapped_df.query(f'channel==1 & particle=={p}').frame.values
        i1 = gapped_df.query(f'channel==1 & particle=={p}').mean_intensity.values

        # if trajectories only exist in one channel skip them
        if len(t1)==0 or len(t0)==0:
            continue

        ##################################################################################
        ## For each particle, fill in the missing gaps
        ##################################################################################
        #max_gap_length = 10
        minT = np.min(np.r_[t0,t1])
        maxT = np.max(np.r_[t0,t1])
        tq = np.arange(minT,maxT+1)

        searchX0 = []
        searchY0 = []
        searchZ0 = []
        searchX1 = []
        searchY1 = []
        searchZ1 = []
        searchT0 = []
        searchT1 = []
        originC0 = []
        originC1 = []

        # identify all non-gaps, store these "known" points
        nongaps1 = list(set(t1).intersection(tq))
        gb = groupby(enumerate(nongaps1), key=lambda x: x[0] - x[1])
        points1 = list([gi[1] for gi in g] for _, g in gb)
        nongaps0 = list(set(t0).intersection(tq))
        gb = groupby(enumerate(nongaps0), key=lambda x: x[0] - x[1])
        points0 = list([gi[1] for gi in g] for _, g in gb)
        for points in points0:
            x0q = np.interp(points,t0,x0)
            y0q = np.interp(points,t0,y0)
            z0q = np.interp(points,t0,z0)
            searchT0.extend(list(points))
            searchX0.extend(list(x0q))
            searchY0.extend(list(y0q))
            searchZ0.extend(list(z0q))
            originC0.extend(['original']*len(points))
        for points in points1:
            x1q = np.interp(points,t1,x1)
            y1q = np.interp(points,t1,y1)
            z1q = np.interp(points,t1,z1)
            searchT1.extend(list(points))
            searchX1.extend(list(x1q))
            searchY1.extend(list(y1q))
            searchZ1.extend(list(z1q))
            originC1.extend(['original']*len(points))

        # identify all gaps and get list of consecutive gaps
        gaps1 = list(set(t1).symmetric_difference(tq))
        gb = groupby(enumerate(gaps1), key=lambda x: x[0] - x[1])
        cgaps1 = list([gi[1] for gi in g] for _, g in gb)
        gaps0 = list(set(t0).symmetric_difference(tq))
        gb = groupby(enumerate(gaps0), key=lambda x: x[0] - x[1])
        cgaps0 = list([gi[1] for gi in g] for _, g in gb)

        # interpolate between small gaps 
        # for large gaaps, get an estimate of frame-to-frame displacement from other channel
        for gap in cgaps0:
            # do interpolation for small gaps
            x0q = np.interp(gap,t0,x0)
            y0q = np.interp(gap,t0,y0)
            z0q = np.interp(gap,t0,z0)
            # if gaps large, or edge cases, use the "other channel" displacements
            if (len(gap) > max_gap_length) or (np.sum(np.diff(x0q)) == 0):

                dx = np.r_[0,np.diff(np.interp(gap,t1,x1))]
                dy = np.r_[0,np.diff(np.interp(gap,t1,y1))]
                dz = np.r_[0,np.diff(np.interp(gap,t1,z1))]

                # is it a left-"edge" case? do things in reverse
                if gap[0] <= t0[0]:
                    # get closest known point
                    x00 = np.interp(gap[-1]+1,t0,x0)
                    y00 = np.interp(gap[-1]+1,t0,y0)
                    z00 = np.interp(gap[-1]+1,t0,z0)
                    # interpolate by the cumulative diferences
                    x0q = x00 -np.cumsum(dx[::-1])
                    y0q = y00 -np.cumsum(dy[::-1])
                    z0q = z00 -np.cumsum(dz[::-1])
                    # reverse the times accordingly
                    gap = gap[::-1]
                else:
                    # get closest known point
                    x00 = np.interp(gap[0]-1,t0,x0)
                    y00 = np.interp(gap[0]-1,t0,y0)
                    z00 = np.interp(gap[0]-1,t0,z0)
                    x0q = x00 + np.cumsum(dx)
                    y0q = y00 + np.cumsum(dy)
                    z0q = z00 + np.cumsum(dz) 
            searchT0.extend(list(gap))
            searchX0.extend(list(x0q))
            searchY0.extend(list(y0q))
            searchZ0.extend(list(z0q))
            originC0.extend(['interp']*len(gap))
        for gap in cgaps1:
            # do interpolation for small gaps
            x1q = np.interp(gap,t1,x1)
            y1q = np.interp(gap,t1,y1)
            z1q = np.interp(gap,t1,z1)
            # if gaps large, or edge cases, use the "other channel" displacements
            if (len(gap) > max_gap_length) or (np.sum(np.diff(x0q)) == 0):

                dx = np.r_[0,np.diff(np.interp(gap,t0,x0))]
                dy = np.r_[0,np.diff(np.interp(gap,t0,y0))]
                dz = np.r_[0,np.diff(np.interp(gap,t0,z0))]

                # is it a left-"edge" case? do things in reverse
                if gap[0] <= t0[0]:
                    # get closest known point
                    x00 = np.interp(gap[-1]+1,t1,x1)
                    y00 = np.interp(gap[-1]+1,t1,y1)
                    z00 = np.interp(gap[-1]+1,t1,z1)
                    # interpolate by the cumulative diferences
                    x1q = x00 -np.cumsum(dx[::-1])
                    y1q = y00 -np.cumsum(dy[::-1])
                    z1q = z00 -np.cumsum(dz[::-1])
                    # reverse the times accordingly
                    gap = gap[::-1]
                else:
                    # get closest known point
                    x00 = np.interp(gap[0]-1,t1,x1)
                    y00 = np.interp(gap[0]-1,t1,y1)
                    z00 = np.interp(gap[0]-1,t1,z1)
                    x1q = x00 + np.cumsum(dx)
                    y1q = y00 + np.cumsum(dy)
                    z1q = z00 + np.cumsum(dz) 
            searchT1.extend(list(gap))
            searchX1.extend(list(x1q))
            searchY1.extend(list(y1q))
            searchZ1.extend(list(z1q))
            originC1.extend(['interp']*len(gap))

        # are there still missing frames?
        assert(len(list(set(searchT1).symmetric_difference(tq)))==0)
        assert(len(list(set(searchT0).symmetric_difference(tq)))==0)

        search_df = pd.DataFrame({'frame':np.r_[searchT0,searchT1],
                                  'channel':np.r_[[0]*len(searchT0),[1]*len(searchT1)],
                                  'particle':np.r_[[p]*len(searchT0),[p]*len(searchT1)],
                                    'x':np.r_[searchX0,searchX1],
                                    'y':np.r_[searchY0,searchY1],
                                    'z':np.r_[searchZ0,searchZ1],
                                    'source':originC0+originC1})
        search_df_list.append(search_df)
    search_df = pd.concat(search_df_list)
#     search_df.to_csv(os.path.join(gapped_pairs_filepath,
#                                      'Inferred_trajectories_'+in_filenames[idx].split('.czi')[0]+'.csv') )
    return search_df


def link_trajectories_across_channels(linked_df, 
                                      min_overlap_length=10,
                                     corrcoeff_min =0.3,
                                     distance_max = 35):
    """
    Links trajectories together across channels
    
    Parameters
    ----------
    linked_df : pandas.DataFrame
        DataFrame of particle trajectories (connected separately for each channel) in 3D and time.         
    
    min_overlap_length : int
        Minimum number of frames of overlap between channels necessary to connect two trajectories
        
    corrcoeff_min : float
        Minimum correlation coefficient between dot pairs (for each dimension separately) in order
        for trajectories to be linked together.
    
    distance_max : int
        Maximum distance in pixels between dots in connected trajectories
        
    Returns
    -------
    joined_df : pandas.DataFrame
        DataFrame of particle trajectories (connected separately for each channel) in 3D and time.  
        This DataFrame may contain some trajectories that do not meet the likely pair criteria. To 
        obtain a clean list, use `linked_trajectory_ids`.
        
    linked_trajectory_ids : pandas.DataFrame
        List of trajectory IDs in `joined_df` meeting the criteria for likely pairs.
    
    
    """
    
    
    fake_linked_df = linked_df.copy()

    # get 4D tracks for Channel 0
    particles0 = list(set(fake_linked_df[fake_linked_df.channel==0].particle))
    tracks_df0 = {}
    for p in particles0:
        X = fake_linked_df[fake_linked_df.particle==p][['frame','x','y','z']].values
        tracks_df0[p] = X

    # get 4D tracks for Channel 1
    particles1 = list(set(fake_linked_df[fake_linked_df.channel==1].particle))
    tracks_df1 = {}
    for p in particles1:
        X = fake_linked_df[fake_linked_df.particle==p][['frame','x','y','z']].values
        tracks_df1[p] = X

    min_id_value = np.min(particles1)
    
    ## Find likely pairs of trajectories
    likely_pairs = []
    for p0i in tracks_df0:
        p0 = tracks_df0[p0i]
        for p1i in tracks_df1:
            p1 = tracks_df1[p1i]
            times, comm0, comm1 = np.intersect1d(p0[:,0],p1[:,0],return_indices=True)
            if len(times) > min_overlap_length:  
                corrcoeffX = np.corrcoef(p0[comm0,1],p1[comm1,1])[0,1]
                corrcoeffY = np.corrcoef(p0[comm0,2],p1[comm1,2])[0,1]
                corrcoeffZ = np.corrcoef(p0[comm0,3],p1[comm1,3])[0,1]
#                 residualsX = np.mean(np.abs(p0[comm0,1]-p1[comm1,1]))
#                 residualsY = np.mean(np.abs(p0[comm0,2]-p1[comm1,2]))
#                 residualsZ = np.mean(np.abs(p0[comm0,3]-p1[comm1,3]))
                residualsX = (np.abs(p0[comm0,1]-p1[comm1,1]))
                residualsY = (np.abs(p0[comm0,2]-p1[comm1,2]))
                residualsZ = (np.abs(p0[comm0,3]-p1[comm1,3]))
                distance = np.sqrt(residualsX**2+residualsY**2+residualsZ**2)

                # link trajectories if they satisfy a minimum correlation coefficient 
                # and are within an expected maximum distance (in nm) of residuals_max
                corrcoeff_good = all([corrcoeff_min < c for c in [corrcoeffX,corrcoeffY,corrcoeffZ]])
                residuals_good = all(distance < distance_max)
                if corrcoeff_good and residuals_good:
                    likely_pairs.append((int(p0i),int(p1i)))

    ################## REFINE THE PAIRING###################
    # cluster trajectories together 
    # find conflicts or overlaps

    clusters = find_clusters(likely_pairs)        
    likelier_pairs = []
    # generate all possible orderings within the cluster
    # for each ordering:
    #     check for overlaps within a channel 
    #      -- if overlap exists, discard this; continue
    #     check for continuous overlap across channels for all segments
    #       -- if overlap is non-contiguous discard; continue
    #     if overlaps are contiguous in time and non-overlapping within a channel
    #       -- score the total length of contig, score average correlation, score average distance; store
    # for each stored/scored ordering
    for ci, c in enumerate(clusters):
        if len(c) < 2:
            continue
        elif len(c)==2:
            likelier_pairs.append(tuple(sorted(c)))
        else:        
            # generate all possible orderings within the cluster
            orderings = partition_select(c)

            # for each ordering
            good_ordering_metrics = {}
            for si, seq in enumerate(orderings):
                # assert the ordering is at least of size 2
                if len(seq) <= 1:
                    #bad_orderings.append(si)
                    continue

                # split the channels and get time points for each group
                group0 = [tracks_df0[s][:,0] for s in seq if s < min_id_value]
                group1 = [tracks_df1[s][:,0] for s in seq if s >= min_id_value]


                # allow for an overlap equal to the number of groups collected -1
                max_overlap0 = 0#len(group0)-1
                max_overlap1 = 0#len(group1)-1

                # CHECK FOR OVERLAPS WITHIN A CHANNEL
                if len(group0) == 0:
                    continue
                elif len(group0)> 1:
                    has_overlap0 = reduce(np.intersect1d, group0)
                else:
                    has_overlap0 = []
                if len(group1) == 0:
                    continue
                elif len(group1)> 1:
                    has_overlap1 = reduce(np.intersect1d, group1)
                else:
                    has_overlap1 = []
                if (len(has_overlap0)>max_overlap0 or len(has_overlap1)>max_overlap1)==True:
                    #bad_orderings.append(si)
                    continue

                # CHECK FOR FULLY CONNECTED TRAJECTORIES
                # count the number of connections for each trajectory
                num_links = {f"group0_{k}":0 for k in range(len(group0))}
                num_links.update({f"group1_{k}":0 for k in range(len(group1))})
                for gi0, g0 in enumerate(group0):
                    for gi1, g1 in enumerate(group1):
                        intersect = np.intersect1d(g0,g1)
                        if len(intersect)> 0:
                            num_links[f"group0_{gi0}"] += 1
                            num_links[f"group1_{gi1}"] += 1
                # if any elements are disconnected, remove this ordering 
                for key in num_links:
                    if num_links[key] == 0:
                        #bad_orderings.append(si)
                        continue

                # WITH THE REMAINING "GOOD ORDERINGS", compute various metrics
                # combine all group0 together
                # combine all group1 together
                times0 = group_arrays(group0)
                times1 = group_arrays(group1)
                times, comm0, comm1 = np.intersect1d(times0,times1,return_indices=True)


                # get X, Y, Z positions for the trajectory
                groupX0 = group_arrays([tracks_df0[s][:,1] for s in seq if s < min_id_value])
                groupX1 = group_arrays([tracks_df1[s][:,1] for s in seq if s >= min_id_value])
                groupY0 = group_arrays([tracks_df0[s][:,2] for s in seq if s < min_id_value])
                groupY1 = group_arrays([tracks_df1[s][:,2] for s in seq if s >= min_id_value])
                groupZ0 = group_arrays([tracks_df0[s][:,3] for s in seq if s < min_id_value])
                groupZ1 = group_arrays([tracks_df1[s][:,3] for s in seq if s >= min_id_value])

                # get only the overlapping subset of trajectory values
                groupX0 = groupX0[comm0]
                groupX1 = groupX1[comm1]
                groupY0 = groupY0[comm0]
                groupY1 = groupY1[comm1]
                groupZ0 = groupZ0[comm0]
                groupZ1 = groupZ1[comm1]

                # get mean correlation coefficient 
                corrcoeff = np.corrcoef(groupX0,groupX0)[0,1]/3 + \
                            np.corrcoef(groupY0,groupY0)[0,1]/3 + \
                            np.corrcoef(groupZ0,groupZ0)[0,1]/3

                # get the mean distance 
                residuals = np.nanmedian(np.sqrt((groupX0-groupX1)**2 + \
                                            (groupY0-groupY1)**2 + \
                                            (groupZ0-groupZ1)**2))

                #print((ci, si))
                # get the overlap length
                good_ordering_metrics[si] = (len(times),residuals,corrcoeff)

            # sort by longest trajectory first
            best_lengths = sorted(good_ordering_metrics, key=lambda k: good_ordering_metrics[k][0],reverse=True)

#             best_trajectory_idx = np.argmin(best_lengths)
            # sort by hightest correlation first
            best_correlations = sorted(good_ordering_metrics, key=lambda k: good_ordering_metrics[k][2],reverse=True)

            # sort by lowest distance first
            best_distance = sorted(good_ordering_metrics, key=lambda k: good_ordering_metrics[k][1],reverse=False)

            # ranked choice vote
            rank_choice =  {k : 0 for k in best_lengths}
            for rank, bl, bc, bd in zip(range(len(best_lengths)),best_lengths,best_correlations,best_distance):
                rank_choice[bl] += rank
                #rank_choice[bc] += rank  # this metric is useless here -- most correlations are close to 1
                #rank_choice[bd] += rank

            # get best ranked trajectory
            min_rank = min(rank_choice.values()) 
            best_trajectory_idx = [key for key in rank_choice if rank_choice[key] == min_rank] 
            best_trajectory_idx = best_trajectory_idx[0]

            # note, these may not be *pairs*, but a tuples of some length 
            likelier_pairs.append(orderings[best_trajectory_idx]) 

    # update likely_pairs with likelier pairs
    likely_pairs = []
    for lp in likelier_pairs:
        lp = sorted(lp)
        if len(lp)==2:
            likely_pairs.append(tuple(lp))
        else:

            #print(lp)
            for lpi in lp[1:]:
                likely_pairs.append((lp[0],lpi) )
    # update particle names
    for p1i, p2i in likely_pairs:
        fake_linked_df['particle'].replace(p2i,p1i,inplace=True)

    linked_particles_ids = list(set([x[0] for x in likely_pairs]))
    return fake_linked_df, linked_particles_ids
    
        
        
def find_clusters( trajectory_pairs ):
    """
    Identifies clusters of trajectories which could be stitched together.
    
    Parameters
    ----------
    trajectory_pairs : list of 2-tuples of ints
        List of 2-tuples of particle IDs, containing potential matched pairs of
        trajectories across channels. This is used for stitching trajectories together.
        
    
    Returns
    -------
    clustered_tuples : n-tuples of ints
        List of n-tuples of particle IDs containing groups of potential trajectory segments.
    
    """
    tuples = trajectory_pairs
    # https://stackoverflow.com/questions/21646703/grouping-elements-from-2-tuples-recursively
    clusterlist=[]
    # clustermap maps an element to the id of the containing
    # cluster within clusterlist
    clustermap = {}

    # here we find the current cluster id for elem, by following the
    # chain within clusterlist, and replace that entire chain
    # with the new cluster id n.   We return the old cluster id.
    def set_cluster_id( elem, n ):
        if elem not in clustermap:
            return None
        k = clustermap[elem]
        # clusters may be replaced by references to other clusters,
        # we follow that chain
        while k < n and isinstance( clusterlist[k], int ):
            k1 = clusterlist[k]
            # this is optional, we make the chain shorter
            # by making this entry point directly to the current cluster
            clusterlist[k] = n
            k = k1
        return k

    for t in tuples:
        # for each tuple we create a new cluster
        thiscluster = set(t)
        n = len( clusterlist ) # the id of thiscluster
        for x in t:
            # we absorb existing clusters into the new one
            # if there is overlap
            k = set_cluster_id(x, n)
            if k is not None and k != n:
                thiscluster.update( clusterlist[k] )
                # we replace the existing cluster
                # with a reference to the new one
                clusterlist[k] = n 
            clustermap[x] = n
        clusterlist.append(thiscluster)

    return [ tuple(x) for x in clusterlist if isinstance( x, set ) ]  
  

def partition_select(vals,min_sz=2):
    """
    This function generates a unique list of all possible partitions
    between the items in `vals`. It also removes all partitions 
    smaller than min_sz in length.
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    
    """
    vals = list(vals)
    
    def _partition(collection):
        #https://stackoverflow.com/questions/19368375/set-partitions-in-python
        if len(collection) == 1:
            yield [ collection ]
            return
        first = collection[0]
        for smaller in _partition(collection[1:]):
            # insert `first` in each of the subpartition's subsets
            for n, subset in enumerate(smaller):
                yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
            # put `first` in its own subset 
            #yield [ [ first ] ] + smaller
            yield [[first]] + smaller
         
    ordering = list(_partition(vals))
    flat_list = []
    for sublist in ordering:
        for item in sublist:
            flat_list.append(tuple(item))
    flat_list = list(set(flat_list))
    return [fl for fl in flat_list if len(fl)>= min_sz]

def group_arrays(list_of_arrays):
    return np.array(reduce(lambda x,y: list(x)+list(y),list_of_arrays))


def filter_large_displacements_in_trackpy_trajectories(linked_df,distance_max=35):
    """
    Trajectory quality check. Removes localizations with big frame-to-frame steps.
    Removal is based on absolute dot positions (i.e. within each channel) and not between channels.
    
    
    Parameters
    ----------
    linked_df : pandas.DataFrame
        DataFrame of dot trajectories and localizations.    
    
    Returns
    -------
    new_df_list : pandas.DataFrame
        A cleaned up version of the `linked_df`, without the large displacements.
    """
    
    new_df_list = []
    channels = linked_df.channel.unique()
    for ch in channels:
        for p in linked_df.particle.unique():
            times = linked_df.query(f"(particle=={p}) and (channel=={ch})").frame.values
            x = linked_df.query(f"(particle=={p}) and (channel=={ch})").x.values
            y = linked_df.query(f"(particle=={p}) and (channel=={ch})").y.values
            z = linked_df.query(f"(particle=={p}) and (channel=={ch})").z.values
            mean_intensity = linked_df.query(f"(particle=={p}) and (channel=={ch})").mean_intensity.values
            max_intensity = linked_df.query(f"(particle=={p}) and (channel=={ch})").max_intensity.values
            dot_size_in_pixels = linked_df.query(f"(particle=={p}) and (channel=={ch})").dot_size_in_pixels.values
            particle = linked_df.query(f"(particle=={p}) and (channel=={ch})").particle.values
            channel = linked_df.query(f"(particle=={p}) and (channel=={ch})").channel.values    

            bad_frames = np.r_[0,np.abs(np.diff(x))< distance_max ] == 1 
            if len(bad_frames) - sum(bad_frames)- 1 > 0:
                # find change points
                changepoints = np.r_[np.where(bad_frames==False)[0],np.inf]

                ind0 = []
                ind1 = []

                ind = 0
                count = 0
                for frame in range(len(x)):
                    if frame < changepoints[count]:
                        if  ind ==0:
                            ind0.append(frame)
                        else:
                            ind1.append(frame)

                    if frame == changepoints[count]:
                        count+=1
                        if  ind ==0:
                            ind1.append(frame)
                            ind = 1
                        else:
                            ind0.append(frame)
                            ind = 0
                ind_keep = np.argmax([len(ind0),len(ind1)])
                if ind_keep == 0:
                    x = [x[pos] for pos in ind0]
                    y = [y[pos] for pos in ind0]
                    z = [z[pos] for pos in ind0]
                    mean_intensity = [mean_intensity[pos] for pos in ind0] 
                    max_intensity = [max_intensity[pos] for pos in ind0]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind0]
                    particle = [particle[pos] for pos in ind0]
                    channel = [channel[pos] for pos in ind0]
                    times = [times[pos] for pos in ind0]
                else:
                    x = [x[pos] for pos in ind1]
                    y = [y[pos] for pos in ind1]
                    z = [z[pos] for pos in ind1]
                    mean_intensity = [mean_intensity[pos] for pos in ind1]
                    max_intensity = [max_intensity[pos] for pos in ind1]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind1]
                    particle = [particle[pos] for pos in ind1]
                    channel = [channel[pos] for pos in ind1]
                    times = [times[pos] for pos in ind1]

            ## filter out big displacements in Y
            bad_frames = np.r_[0,np.abs(np.diff(y))< distance_max ] == 1 
            if len(bad_frames) - sum(bad_frames)- 1 > 0:
                # find change points
                changepoints = np.r_[np.where(bad_frames==False)[0],np.inf]

                ind0 = []
                ind1 = []

                ind = 0
                count = 0
                for frame in range(len(x)):
                    if frame < changepoints[count]:
                        if  ind ==0:
                            ind0.append(frame)
                        else:
                            ind1.append(frame)

                    if frame == changepoints[count]:
                        count+=1
                        if  ind ==0:
                            ind1.append(frame)
                            ind = 1
                        else:
                            ind0.append(frame)
                            ind = 0
                ind_keep = np.argmax([len(ind0),len(ind1)])
                if ind_keep == 0:
                    x = [x[pos] for pos in ind0]
                    y = [y[pos] for pos in ind0]
                    z = [z[pos] for pos in ind0]
                    mean_intensity = [mean_intensity[pos] for pos in ind0] 
                    max_intensity = [max_intensity[pos] for pos in ind0]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind0]
                    particle = [particle[pos] for pos in ind0]
                    channel = [channel[pos] for pos in ind0]
                    times = [times[pos] for pos in ind0]
                else:
                    x = [x[pos] for pos in ind1]
                    y = [y[pos] for pos in ind1]
                    z = [z[pos] for pos in ind1]
                    mean_intensity = [mean_intensity[pos] for pos in ind1]
                    max_intensity = [max_intensity[pos] for pos in ind1]
                    dot_size_in_pixels = [dot_size_in_pixels[pos] for pos in ind1]
                    particle = [particle[pos] for pos in ind1]
                    channel = [channel[pos] for pos in ind1]
                    times = [times[pos] for pos in ind1]

            new_df = pd.DataFrame({'particle':particle,
                          'channel':channel,
                          'frame':times,
                          'x':x,
                          'y':y,
                          'z':z,
                          'mean_intensity':mean_intensity,
                          'max_intensity':max_intensity,
                          'dot_size_in_pixels':dot_size_in_pixels})

            new_df_list.append(new_df)
        
    return pd.concat(new_df_list)


def infer_missing_dot_locations(gapped_df,max_gap_length = 10):
    """
    Takes existing joint trajectories dataframes and interpolates to guess where 
    dots may be in frames where localizations are missing.
    
    Parameters
    ----------
    gapped_df : pandas.DataFrame
       DataFrame of trajectories that have been linked across channels.
    
    max_gap_length : int
        Maximum number of frames with which to use within channel interpolation. If 
        the number of missing frames exceeds this value, displacements from the 
        adjacent channel are used.
        
    Returns
    -------
    search_df : pandas.DataFrame
        DataFrame containing best guesses of where a dot is likely to be in a missing frame.
    """
    
    particles = list(set(gapped_df.particle))
    search_df_list = []
    for p in particles:

        x0 = gapped_df.query(f'channel==0 & particle=={p}').x.values
        y0 = gapped_df.query(f'channel==0 & particle=={p}').y.values
        z0 = gapped_df.query(f'channel==0 & particle=={p}').z.values
        t0 = gapped_df.query(f'channel==0 & particle=={p}').frame.values
        i0 = gapped_df.query(f'channel==1 & particle=={p}').mean_intensity.values

        x1 = gapped_df.query(f'channel==1 & particle=={p}').x.values
        y1 = gapped_df.query(f'channel==1 & particle=={p}').y.values
        z1 = gapped_df.query(f'channel==1 & particle=={p}').z.values
        t1 = gapped_df.query(f'channel==1 & particle=={p}').frame.values
        i1 = gapped_df.query(f'channel==1 & particle=={p}').mean_intensity.values

        # if trajectories only exist in one channel skip them
        if len(t1)==0 or len(t0)==0:
            continue

        ##################################################################################
        ## For each particle, fill in the missing gaps
        ##################################################################################
        #max_gap_length = 10
        minT = np.min(np.r_[t0,t1])
        maxT = np.max(np.r_[t0,t1])
        tq = np.arange(minT,maxT+1)

        searchX0 = []
        searchY0 = []
        searchZ0 = []
        searchX1 = []
        searchY1 = []
        searchZ1 = []
        searchT0 = []
        searchT1 = []
        originC0 = []
        originC1 = []

        # identify all non-gaps, store these "known" points
        nongaps1 = list(set(t1).intersection(tq))
        gb = groupby(enumerate(nongaps1), key=lambda x: x[0] - x[1])
        points1 = list([gi[1] for gi in g] for _, g in gb)
        nongaps0 = list(set(t0).intersection(tq))
        gb = groupby(enumerate(nongaps0), key=lambda x: x[0] - x[1])
        points0 = list([gi[1] for gi in g] for _, g in gb)
        for points in points0:
            x0q = np.interp(points,t0,x0)
            y0q = np.interp(points,t0,y0)
            z0q = np.interp(points,t0,z0)
            searchT0.extend(list(points))
            searchX0.extend(list(x0q))
            searchY0.extend(list(y0q))
            searchZ0.extend(list(z0q))
            originC0.extend(['original']*len(points))
        for points in points1:
            x1q = np.interp(points,t1,x1)
            y1q = np.interp(points,t1,y1)
            z1q = np.interp(points,t1,z1)
            searchT1.extend(list(points))
            searchX1.extend(list(x1q))
            searchY1.extend(list(y1q))
            searchZ1.extend(list(z1q))
            originC1.extend(['original']*len(points))

        # identify all gaps and get list of consecutive gaps
        gaps1 = list(set(t1).symmetric_difference(tq))
        gb = groupby(enumerate(gaps1), key=lambda x: x[0] - x[1])
        cgaps1 = list([gi[1] for gi in g] for _, g in gb)
        gaps0 = list(set(t0).symmetric_difference(tq))
        gb = groupby(enumerate(gaps0), key=lambda x: x[0] - x[1])
        cgaps0 = list([gi[1] for gi in g] for _, g in gb)

        # interpolate between small gaps 
        # for large gaaps, get an estimate of frame-to-frame displacement from other channel
        for gap in cgaps0:
            # do interpolation for small gaps
            x0q = np.interp(gap,t0,x0)
            y0q = np.interp(gap,t0,y0)
            z0q = np.interp(gap,t0,z0)
            # if gaps large, or edge cases, use the "other channel" displacements
            if (len(gap) > max_gap_length) or (np.sum(np.diff(x0q)) == 0):

                dx = np.r_[0,np.diff(np.interp(gap,t1,x1))]
                dy = np.r_[0,np.diff(np.interp(gap,t1,y1))]
                dz = np.r_[0,np.diff(np.interp(gap,t1,z1))]

                # is it a left-"edge" case? do things in reverse
                if gap[0] <= t0[0]:
                    # get closest known point
                    x00 = np.interp(gap[-1]+1,t0,x0)
                    y00 = np.interp(gap[-1]+1,t0,y0)
                    z00 = np.interp(gap[-1]+1,t0,z0)
                    # interpolate by the cumulative diferences
                    x0q = x00 -np.cumsum(dx[::-1])
                    y0q = y00 -np.cumsum(dy[::-1])
                    z0q = z00 -np.cumsum(dz[::-1])
                    # reverse the times accordingly
                    gap = gap[::-1]
                else:
                    # get closest known point
                    x00 = np.interp(gap[0]-1,t0,x0)
                    y00 = np.interp(gap[0]-1,t0,y0)
                    z00 = np.interp(gap[0]-1,t0,z0)
                    x0q = x00 + np.cumsum(dx)
                    y0q = y00 + np.cumsum(dy)
                    z0q = z00 + np.cumsum(dz) 
            searchT0.extend(list(gap))
            searchX0.extend(list(x0q))
            searchY0.extend(list(y0q))
            searchZ0.extend(list(z0q))
            originC0.extend(['interp']*len(gap))
        for gap in cgaps1:
            # do interpolation for small gaps
            x1q = np.interp(gap,t1,x1)
            y1q = np.interp(gap,t1,y1)
            z1q = np.interp(gap,t1,z1)
            # if gaps large, or edge cases, use the "other channel" displacements
            if (len(gap) > max_gap_length) or (np.sum(np.diff(x0q)) == 0):

                dx = np.r_[0,np.diff(np.interp(gap,t0,x0))]
                dy = np.r_[0,np.diff(np.interp(gap,t0,y0))]
                dz = np.r_[0,np.diff(np.interp(gap,t0,z0))]

                # is it a left-"edge" case? do things in reverse
                if gap[0] <= t0[0]:
                    # get closest known point
                    x00 = np.interp(gap[-1]+1,t1,x1)
                    y00 = np.interp(gap[-1]+1,t1,y1)
                    z00 = np.interp(gap[-1]+1,t1,z1)
                    # interpolate by the cumulative diferences
                    x1q = x00 -np.cumsum(dx[::-1])
                    y1q = y00 -np.cumsum(dy[::-1])
                    z1q = z00 -np.cumsum(dz[::-1])
                    # reverse the times accordingly
                    gap = gap[::-1]
                else:
                    # get closest known point
                    x00 = np.interp(gap[0]-1,t1,x1)
                    y00 = np.interp(gap[0]-1,t1,y1)
                    z00 = np.interp(gap[0]-1,t1,z1)
                    x1q = x00 + np.cumsum(dx)
                    y1q = y00 + np.cumsum(dy)
                    z1q = z00 + np.cumsum(dz) 
            searchT1.extend(list(gap))
            searchX1.extend(list(x1q))
            searchY1.extend(list(y1q))
            searchZ1.extend(list(z1q))
            originC1.extend(['interp']*len(gap))

        # are there still missing frames?
        assert(len(list(set(searchT1).symmetric_difference(tq)))==0)
        assert(len(list(set(searchT0).symmetric_difference(tq)))==0)

        search_df = pd.DataFrame({'frame':np.r_[searchT0,searchT1],
                                  'channel':np.r_[[0]*len(searchT0),[1]*len(searchT1)],
                                  'particle':np.r_[[p]*len(searchT0),[p]*len(searchT1)],
                                    'x':np.r_[searchX0,searchX1],
                                    'y':np.r_[searchY0,searchY1],
                                    'z':np.r_[searchZ0,searchZ1],
                                    'source':originC0+originC1})
        search_df_list.append(search_df)
    search_df = pd.concat(search_df_list)
#     search_df.to_csv(os.path.join(gapped_pairs_filepath,
#                                      'Inferred_trajectories_'+in_filenames[idx].split('.czi')[0]+'.csv') )
    return search_df