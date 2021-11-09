## Import required modules
import matplotlib.pyplot as plt # for plotting
import matplotlib # for plotting
import numpy as np # for manipulating arrays
import os # for making/deleting directories
import bioformats # for reading image series
import javabridge # for interfacing with java (required for bioformats)
from tifffile import xml2dict # for parsing the metadata from bioformats
import pickle # for saving python objects and other data
from scipy.optimize import curve_fit # for making fits to the PSF
from scipy.ndimage import gaussian_laplace, gaussian_filter # for dot localization (image filtering)
from skimage import measure # for segmenting images
from skimage.morphology import remove_small_objects, closing, disk # for morphological filtering of images
from skimage.segmentation import clear_border # for filtering images
from skimage.filters import threshold_otsu
import pandas as pd # for creating and manipulating tabulated data
from collections import Iterable
from itertools import product
import copy
import scipy

# settings for making nice pdfs
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.sans-serif'] = "DejaVu Sans"
plt.rcParams['font.family'] = "sans-serif"

javabridge.start_vm(class_path=bioformats.JARS) # start java virtual machine

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

def filter_zstack_DoG(zstack,dog_sigma1 = 1.5,dog_sigma2 = 15,absolute_value=True):
    """
    Applies Difference of Gaussian (DoG) filtering on a single z-stack. 

    Parameters
    ----------
    zstack : numpy.ndarray [sizeY by sizeX by sizeZ]
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

def get_image_threshold(image,method,**kwargs):
    """
    Returns a threshold value for binarizing an image for morphological filtering and dot localization.
    
    Parameters
    ----------
    image : numpy.ndarray [sizeY by sizeX]
    
    method : str {'otsu','percentile'}  
    
    kwargs : For method 'otsu'  
             `nbins` : int (optinal) 
             number of bins used for otsu method
             
             For method 'percentile' 
             `percentile_threshold` : float
             value ranging from 0 to 100
    
    Returns
    -------   
    threshold : float
                Value of threshold determined by the specified method. By default, it is the
                99th percentile of pixel intensities of the image.
                
    """   
    method = method.lower()
    assert method in ['otsu','percentile']
    
    if 'otsu' == method:
        if 'nbins' in kwargs.keys():
            threshold = threshold_otsu(image,kwargs['nbins'])
        else:
            threshold = threshold_otsu(image)
    else: #'percentile' == method:
        if 'percentile_threshold' in kwargs.keys():
            threshold = np.percentile(image,kwargs['percentile_threshold'])
        else:
            threshold = np.percentile(image,99)
        
    return threshold

def localize_dots_XY_projection(filtered_zstack, min_object_area=50,\
                                intensity_threshold=None, projectionAxis=2):
    """
    Roughly localizes dots in maximum projection image using morphological filtering.

    Parameters
    ----------
    filtered_zstack : numpy.ndarray [sizeY by sizeX by sizeZ]
                      Z-stack containing 3 spatial dimentions.
             
    min_object_area : float, optional
                      Minimum area (in pixels) of the object being localized.
                      
                      
    intensity_threshold : float, optional
                          Threshold value by which to binarize the image. By default, this value will
                          be the 99th percentile of pixel intensity values of the maximum projection 
                          image. For other ways to choose the `intensity_threshold` value, we refer to: 
                          skimage.filters (e.g. threshold_otsu).
                
    projectionAxis : {2,1,0}, optional
                     Value of the dimension along which to compute the maximum intensity projection. 
                     The default is 2 (i.e. removes the Z-dimension).

    Returns
    -------
    centroids : list of ints
                List of integer pixel values close to the centroid of each located "dot" in the
                maximum intensity projection image
    
    (blobs, blobs_labels,blob_regions) : numpy.ndarray, numpy.ndarray, list of RegionProperties
            `blobs` is the thresholded, morphologically filtered maximum intensity projection image. 
            `blobs_labels` is the segmentation of the image after connecting proximal pixels.
            `blob_metrics` is an object containing a list of measured attribues for each unique
             region of `blobs_labels`; `blob_metrics` is the output of skimage.measure.regionprops().
    """
    
    max_proj = np.max(filtered_zstack,axis=projectionAxis) # get maximum intensity projection
    
    if intensity_threshold is None:
        intensity_threshold = np.percentile(max_proj,99)  
        
    blobs = max_proj  >  intensity_threshold # binarize image based on global threshold
    
    # filter objects based on size
    blobs = remove_small_objects(blobs, min_size=min_object_area) 
    
    # remove objects touching the edges of the image
    blobs = clear_border(blobs) 
    
    # "closing" operation to connect proximal pixels
    # blobs =  closing(blobs > intensity_threshold, disk(2)) 

    # get segmentation of the image from connected pixels 
    blobs_labels = measure.label(blobs, background=0) 

    # measure things for each unique feature identified in blobs_labels
    blob_metrics = measure.regionprops(blobs_labels, max_proj )

    # get centroids of objects. i.e. (x,y) coordinates
    # note that the values are actually returned as (y,x) coordinates
    centroids = [tuple(np.array(x.weighted_centroid,dtype=int)) for x in blob_metrics]
    
    return centroids, (blobs, blobs_labels,blob_metrics)


def fit_Gaussian_3D_PSF(zstack, dot_positions_xy, window_size=10,\
                        do_classification=False,do_gaussian_fitting=False,verbose=False):
    """
    Fits specified dots in zstack to 3D Gaussian function.

    Parameters
    ----------
    zstack : numpy.ndarray [sizeY by sizeX by sizeZ]
            Original Z-stack from the image series.
             
    dot_positions_xy : list of 2-tuples of ints
            List of approximate (X,Y) positions of dots in the Z-stack.
    
    window_size : int, optional 
            Length of area used to crop features out of the z-stack. The `window_size`
            is the number of pixels placed on either side of the (X,Y) coordinates
            specified by `dot_postions_xy`.
                  
    do_classification : {T,F}
            Classifies the number of modes (i.e. number of unique features) in each cropped image.  
                  
    do_gaussian_fitting : {T,F}
            If True, a true 3D PSF is fit to the data, otherwise, maximum intensity & x,y,z positions are
            returned, and guesses for the variances. 
                  
    Returns
    -------
    dot_fits_dict : dict            
                    Contains 3D PSF parameter fit values, and other metrics used 
                    for quality control of the fit and feature localization.
                    
                    Attributes of `dot_fits_dict`. 
                    'max_projection_xy_data' : maximum intensity projection of the data (XY plane)
                    'max_projection_xz_data' : maximum intensity projection of the data (XZ plane)
                    'max_projection_yz_data' : maximum intensity projection of the data (YZ plane)
                    'max_projection_xy_fit' : maximum intensity projection of the fit (XY plane)
                    'max_projection_xz_fit' : maximum intensity projection of the fit (XZ plane)
                    'max_projection_yz_fit' : maximum intensity projection of the fit (YZ plane)
                    'I0_fit' : maximum intensity of the dot (from fit)
                    'wxy_fit' : standard deviation of the dot along the x and y dimensions (from fit)
                    'wz_fit' : standard deviation of the dot along the z dimension (from fit)
                    'x0_fit' : x dimension best fit value for dot center  
                    'y0_fit' : y dimension best fit value for dot center  
                    'z0_fit' : z dimension best fit value for dot center 
                    'pcov' : covariance matrix for the parameters 
                            (I0_fit,wxy_fit,wz_fit,x0_fit,y0_fit,z0_fit)
                    'num_modes' : number of modes identified in `max_projection_{}_data` image 
    """

    dot_fits_dict = {}
    win = window_size
    
    for di, (xc, yc) in enumerate(dot_positions_xy):
        
        # skip points too close to the frame edge
        sizeX = zstack.shape[0]
        sizeY = zstack.shape[1]
        sizeZ = zstack.shape[2]
        if (xc < win) or (xc >= sizeX-win) or (yc < win) or (yc >= sizeY-win):
            continue
        
        # crop out the "dot" from the zstack
        dot_volume =  zstack[xc-win:xc+win,yc-win:yc+win,:]
        
        # flatten the voxels around the dot for fitting purposes
        flat_vol =  np.ndarray.flatten(dot_volume)  
    
        # define the 3D PSF kernel (for plotting)
        def _gauss3D(I0,wxy,wz,x0,y0,z0,background):
            xx = np.arange(xc-win,xc+win)
            yy = np.arange(yc-win,yc+win)
            zz = np.arange(0,sizeZ)
            xmesh,ymesh,zmesh = np.meshgrid(xx, yy,zz, sparse=True)
            divxy = 2*wxy**2
            divz = 2*wz**2
            prefactor = (2*np.pi)**1.5*wxy**2*wz
            return I0*np.exp(-((xmesh-x0)**2+(ymesh-y0)**2)/divxy-(zmesh-z0)**2/divz)/prefactor

        # define the 3D PSF kernel (for fitting)
        def _gauss3D_fit(self,I0,wxy,wz,x0,y0,z0,background):
            xx = np.arange(xc-win,xc+win)
            yy = np.arange(yc-win,yc+win)
            zz = np.arange(0,sizeZ)
            xmesh,ymesh,zmesh = np.meshgrid(xx, yy,zz, sparse=True)
            divxy = 2*wxy**2
            divz = 2*wz**2
            prefactor = (2*np.pi)**1.5*wxy**2*wz
            gauss_ker = I0*np.exp(-((xmesh-x0)**2+(ymesh-y0)**2)/divxy-(zmesh-z0)**2/divz)/prefactor+background
            return np.ndarray.flatten(gauss_ker)

        # generate initial guess of fit values for the curve fitting algorithm 
        I0_guess = np.max(dot_volume)
        wxy_guess = 2
        wz_guess = 0.5
        
        # refine original "centroid" coordinates with a better guess
        yc_rel, xc_rel, zc_rel = np.unravel_index(np.argmax(dot_volume, axis=None), dot_volume.shape)
        yc_guess = yc + yc_rel - window_size
        xc_guess = xc + xc_rel - window_size
        zc_guess = zc_rel 
        
        if do_gaussian_fitting == True:
            # add background parameter to the fit
            background_guess = np.median(dot_volume)
            initial_guess = [I0_guess,wxy_guess,wz_guess,xc_guess,yc_guess,zc_guess,background_guess]

            # place bounds on the fitting parameters
            unc = 2 # pixel uncertainty on the centroid position
            maxI = np.max(dot_volume)
            minI = np.min(dot_volume)
            lower_bounds = [minI,0,0,xc_guess-unc,yc_guess-unc,zc_guess-unc,minI]
            upper_bounds = [maxI,window_size,window_size,\
                            xc_guess+unc,yc_guess+unc,zc_guess+unc,maxI]

            # get the fit parameters
            try:
                (I0_fit,wxy_fit,wz_fit,x0_fit,y0_fit,z0_fit,background_fit), pcov = \
                            curve_fit(_gauss3D_fit,flat_vol, flat_vol,p0=initial_guess,\
                                     bounds=(lower_bounds,upper_bounds))
            except:
                if verbose == True: 
                    print('failed at dot {}'.format(di))
                continue
        else:
            I0_fit = I0_guess
            wxy_fit = wxy_guess
            wz_fit = wz_guess
            x0_fit = xc_guess
            y0_fit = yc_guess
            z0_fit = zc_guess
            background_fit = 0
            pcov = []
            
        # generate the fit volume
        fit_psf = _gauss3D(I0_fit,wxy_fit,wz_fit,x0_fit,y0_fit,z0_fit,background_fit)

        max_projection_xy_data = np.max(dot_volume,axis=2)
        max_projection_xz_data = np.max(dot_volume,axis=0)
        max_projection_yz_data = np.max(dot_volume,axis=1)
        max_projection_xy_fit = np.max(fit_psf,axis=2)
        max_projection_xz_fit = np.max(fit_psf,axis=0)
        max_projection_yz_fit = np.max(fit_psf,axis=1)
            
        # write maximum projection data and fits to dictionary
        dot_fits_dict[di] = {'max_projection_xy_data':max_projection_xy_data,\
                          'max_projection_xz_data':max_projection_xz_data, \
                          'max_projection_yz_data':max_projection_yz_data, \
                          'max_projection_xy_fit':max_projection_xy_fit, \
                          'max_projection_xz_fit':max_projection_xz_fit, \
                          'max_projection_yz_fit':max_projection_yz_fit, \
                          'I0_fit':I0_fit,\
                          'wxy_fit':wxy_fit,\
                          'wz_fit':wz_fit,\
                          'x0_fit':x0_fit,\
                          'y0_fit':y0_fit,\
                          'z0_fit':z0_fit,\
                          'pcov':pcov,\
                          'num_modes': {}}  
        
        # classify the number of modes in each maximum projection data image
        if do_classification == True:
            num_modes = {}
            for img_key in ['max_projection_xy_data','max_projection_xz_data','max_projection_yz_data']:
                img = dot_fits_dict[di][img_key]
                dot_fits_dict[di]['num_modes'].update({img_key : count_dots_from_threshold(img)}) 
            
    return dot_fits_dict

def do_one_frame(filename,frame, channel=0, img_info=None, dog_sigma1=1.5, dog_sigma2=3, \
                 min_object_area=50, intensity_threshold_method='percentile', 
                 window_size=10, classify_dots=True,do_gaussian_fitting=False, load_file_path=None, \
                 save_intermediates_file_path=None, return_intermediates=False,verbose=False,**kwargs ):
    """
    Localizes dots and performs 3D PSF fitting on a single frame (z-stack)  

    Parameters
    ----------
    filename : str
               Name of the file from which to retrieve the z-stack.

    frame : int
               The temporal slice of the image series from which to retrieve the z-stack.
            
    channel : int, optional
               The channel from which to retrieve the z-stack.
               
    img_info : tuple of ints, optional 
                Pre-retrieved metadata for increased speed in batch processing.
                5-tuple containing lengths of the `X`, `Y`, `Z` (spatial), `T` (temporal) 
                dimensions of the image series, and the number of channels, `num_channels`.
                See output of get_CZI_metadata(). 
          
    
    dog_sigma1 : float, optional
                Standard deviation of the first Gaussian distribution of the DoG filter. 
                `dog_sigma1` should be close in size to the "dots" being tracked.
                See filter_zstack_DoG().
             
    dog_sigma2 : float, optional
                Standard deviation of the second Gaussian distribution of the DoG filter. 
                `dog_sigma2` should be larger than `dog_sigma1`. See filter_zstack_DoG(). 
                
    min_object_area : float, optional
                Minimum area (in pixels) of the object being localized. 
                See localize_dots_XY_projection().                      
                      
    intensity_threshold_method : str, optional
                Method of selecting the threshold value by which to binarize the filtered z-stack 
                image. By default, the method is 'percentile', and will use the 99th percentile 
                of pixel intensity values. For other methods, see get_image_threshold().

    window_size : int, optional 
                Length of area used to crop features out of the z-stack. The `window_size`
                is the number of pixels placed on either side of localized dot centroids.
                See fit_Gaussian_3D_PSF()
                  
    classify_dots : {T,F}
                Counts the number of dots found in each cropped feature (of window size
                defined by `window_size`). 
                
    load_file_path : str, optional
                Path to the file from which to retrieve the z-stack.
                    
    save_intermediates_file_path : str, optional 
                Path to a folder in which to save intermediate results from the analysis.
                Intermediates saved will include `dot_fits_dict`, `blobs`, `filtered_zstack`.
                If the specified folder does not exist, it is created.
    
    return_intermediates : {T,F}, optional
                Option to return not only `fits_df` but also the intermediates including 
                `dot_fits_dict`, `blobs`, `blobs_labels`, `blob_metrics` and `filtered_zstack`.
    
    verbose : {T,F}, optional
                Prints to standard output the steps being performed.
                
    **kwargs : optional
                Pass key word arguments. For example to get_image_threshold() to specify 
                parameters for the thresholding method (e.g. if `intensity_threshold_method`
                is 'percentile', one can optionally pass `percentile_threshold=90` to threshold
                at the 90th percentile instead of the default of 99th percentile).
            
    Returns
    -------            
    fits_df : pandas DataFrame
                DataFrame containing information on the X,Y,Z PSF localization, frame number,
                channel and intensity of each localized dot in the z-stack.
                
    Additionally Returns (if `return_intermediates`== True): 
    --------------------
    
    dot_fits_dict : dict            
                Contains 3D PSF parameter fit values, and other metrics used 
                for quality control of the fit and feature localization.

                Attributes of `dot_fits_dict`. 
                'max_projection_xy_data' : maximum intensity projection of the data (XY plane)
                'max_projection_xz_data' : maximum intensity projection of the data (XZ plane)
                'max_projection_yz_data' : maximum intensity projection of the data (YZ plane)
                'max_projection_xy_fit' : maximum intensity projection of the fit (XY plane)
                'max_projection_xz_fit' : maximum intensity projection of the fit (XZ plane)
                'max_projection_yz_fit' : maximum intensity projection of the fit (YZ plane)
                'I0_fit' : maximum intensity of the dot (from fit)
                'wxy_fit' : standard deviation of the dot along the x and y dimensions (from fit)
                'wz_fit' : standard deviation of the dot along the z dimension (from fit)
                'x0_fit' : x dimension best fit value for dot center  
                'y0_fit' : y dimension best fit value for dot center  
                'z0_fit' : z dimension best fit value for dot center 
                'pcov' : covariance matrix for the parameters 
                        (I0_fit,wxy_fit,wz_fit,x0_fit,y0_fit,z0_fit) 
                'num_modes' : dict; key is `max_projection_{}_data`, value is # modes found in image 
    
    zstack : numpy.ndarray [sizeY by sizeX by sizeZ]
             Z-stack of the image series for a single channel (containing 3 spatial dimentions)
    
    filtered_zstack : numpy.ndarray
                      Absolute value of Difference of Gaussian filtered z-stack. 
    
    centroids : list of ints
                List of integer pixel values close to the centroid of each located "dot" in the
                maximum intensity projection image
    
    blobs_labels : numpy.ndarray
                `blobs_labels` is the segmentation of the image after thresholding, morphologically 
                filtering and connecting proximal pixels of `filtered_zstack_max_projection` image. 
                
    blob_metrics : object
                Metrics for each `blob_labels` region can be obtained from skimage.measure.regionprops(). 
                
                     
    Saved to disk (if `save_intermediates_file_path` is provided)
    -------------
    
    fits_df : pandas DataFrame 
                (see `fits_df` above)
    
    dot_fits_dict : dict 
                (see `dot_fits_dict` above)            
    
    filtered_zstack_max_projection : numpy.ndarray [sizeY by sizeX]
                Maximum projection of the filtered Z-stack onto the X,Y plane.
    
    centroids : list of ints 
                (see `centroids` above)
    
    blobs_labels : numpy.ndarray 
                (see `blob_labels` above)
    
    """      
    
    
    # loads a z-stack from a single channel from the specified frame 
    if verbose: print("\nLoading: {}\nFrame {} Channel {}".format(filename,frame,channel))
    zstack = get_CZI_zstack(filename,frame,channel,filepath=load_file_path,img_info=img_info)
    
    # apply DoG filter 
    if verbose: print("1) Appling Difference of Gaussian filter to z-stack.")
    filtered_zstack = filter_zstack_DoG(zstack,dog_sigma1=dog_sigma1,dog_sigma2=dog_sigma2)
    
    # obtain image threshold
    if verbose: print("2) Obtaining image threshold using method: {}".format(intensity_threshold_method))
    thresh = get_image_threshold(np.max(filtered_zstack,2), intensity_threshold_method, **kwargs)
    
    # localize dots [(X,Y) coordinates] from the maximum projection of `filtered_zstack`
    if verbose: print("3) Morphological filtering and localizing dots in 2D.")
    loc_output = localize_dots_XY_projection(filtered_zstack, min_object_area=min_object_area,\
                      intensity_threshold=thresh, projectionAxis=2)
    centroids, (blobs, blobs_labels,blob_metrics) = loc_output 
    
    # do 3D PSF fitting
    if verbose: print("4) 3D PSF fitting and localizing dots in 3D.")
    dot_fits_dict = fit_Gaussian_3D_PSF(zstack, centroids, window_size=window_size,\
                                        do_classification=classify_dots, \
                                        do_gaussian_fitting=do_gaussian_fitting, verbose=False)
    
    
    
    # unpack `dot_fits_dict` fit values into a pandas DataFrame               
    if verbose: print("5) Generating pandas DataFrame from the PSF fits")                  
    I0_fits = []
    wxy_fits = []
    wz_fits = []
    x0_fits = []
    y0_fits = []
    z0_fits = []
    num_modes = []
    for key in dot_fits_dict.keys():
        I0_fits.append(dot_fits_dict[key]['I0_fit'])
        wxy_fits.append(dot_fits_dict[key]['wxy_fit'])
        wz_fits.append(dot_fits_dict[key]['wz_fit'])
        x0_fits.append(dot_fits_dict[key]['x0_fit'])
        y0_fits.append(dot_fits_dict[key]['y0_fit'])
        z0_fits.append(dot_fits_dict[key]['z0_fit'])
        modes_list = [dot_fits_dict[key]['num_modes'][k] for k in dot_fits_dict[key]['num_modes'].keys()]
        num_modes.append(np.mean(modes_list)) # allow for 1 false-positive
        
    colnames = ['channel','frame','x','y','z','intensity','avg_num_modes']
    channels = [int(channel)]*len(I0_fits)
    frames = [int(frame)]*len(I0_fits)
                      
    fits_df = pd.DataFrame([channels,frames,x0_fits,y0_fits,z0_fits,I0_fits,num_modes],index=colnames).transpose()
    
    # save intermediates
    if not save_intermediates_file_path is None:
        if verbose: print("6) Saving intermediates.")             
        # check if folder exists, if not make it
        if not os.path.exists(save_intermediates_file_path):
            os.makedirs(save_intermediates_file_path)
                      
        # generate file name and save Data Frame to csv filel
        data_frame_filename = os.path.join(save_intermediates_file_path,\
                                           'frame{}_channel{}_dotFitsDict.csv'.format(frame,channel))
        fits_df.to_csv(data_frame_filename)
                      
        # generate file name and save other information to pickled object
        analysis_intermediates = {'filtered_zstack_max_projection': np.max(filtered_zstack,2),\
                                 #'blobs': blobs,\
                                 'blobs_labels': blobs_labels,\
                                 #'blob_metrics': blob_metrics,\
                                 'centroids': centroids,\
                                 'dot_fits_dict': dot_fits_dict}
        intermediates_filename = os.path.join(save_intermediates_file_path,\
                               'frame{}_channel{}_AnalysisIntermediates.pkl'.format(frame,channel))
        pickle.dump(analysis_intermediates, open(intermediates_filename,'wb'))
        
    if return_intermediates == True: 
        return fits_df, dot_fits_dict, zstack, filtered_zstack, centroids, blobs_labels, blob_metrics
    else: 
        return fits_df 


def do_all_frames_all_channels(save_output_path=None,**params_dict):
    """
    Parameters
    ----------
    params_dict : dict
                Dictionary containing 3 mandatory components:
                1) filename : str
                              Name of the image series to be analyzed
                2) filepath : str
                              Full file path to the image series
                3) channel_params_dict : dict (or list of dicts)
                              A dictionary of keyword arguments that will be passed to do_one_frame().
                              If `channel` is anonymous (i.e. if it is not specified in the arguments), 
                              the parameters passed  to do_one_frame() will be the same for all channels. 
                              If `channel_params_dict` contains a list of dictionaries (e.g. specifying 
                              parameters for each channel), the channel-specific parameters will be passed
                              to do_one_frame(); if channel-specific parameters are unspecified, the 
                              function attempts to use parameters from the last anonymous channel. Otherwise,
                              function defaults are used.
                      
    save_output_path : str, optional
                Path specifying where to save `df_all` to disk
                
    Example
    -------
    E.g. params_dict = {'filename': 'my_file.czi',
                        'filepath': './file_location/',
                        'channel_params_dict': [{'channel': 0,
                                              'dog_sigma1': 1.5, 
                                              'dog_sigma2': 15, 
                                              'min_object_area': 50, 
                                              'intensity_threshold_method': 'percentile', 
                                              'percentile_threshold': 99,
                                              'window_size': 10,
                                              'save_intermediates_file_path': './tmp'},

                                               {'channel': 1,
                                               'dog_sigma1': 1.5, 
                                               'dog_sigma2': 3, 
                                               'min_object_area': 35, 
                                               'intensity_threshold_method': 'percentile', 
                                               'percentile_threshold': 99,
                                               'window_size': 10,
                                               'save_intermediates_file_path': './tmp'},
                                                ]
                         }
    Returns
    -------
    df_all : pandas DataFrame
            DataFrame containing 3D PSF fit information for all localized dots in the image series.
            This DataFrame is structured such that it can be passed to .
    
    Saving intermediate steps
    --------------------------
    Note. Intermediate outputs are saved to disk if `save_intermedites_file_path` is specified. 
          Saved outputs are structured as specified in do_one_frame().
    """

    """
    To do: fix depracated: channel_params_dict = params_dict['channel_params_dict'][channel]
    """
    
    df_list = []

    filename = params_dict['filename']
    filepath = params_dict['filepath']

    # get metadata
    img_info, _ = get_CZI_metadata(filename,filepath)
    (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info

    # iterate over channels and frames
    for (channel,frame) in product(range(num_channels), range(sizeT)):

        # retrieve fitting parameters for specified channel
        channel_params_dict = params_dict['channel_params_dict']
        this_channel_params_dict = None
        
        # search list for channel-specific parameters dicts
        if type(channel_params_dict)==list:
            for d in channel_params_dict:
                if 'channel' in d:
                    if d['channel'] == channel:
                        # use channel-specific parameters
                        this_channel_params_dict = d
                        break
                else:
                    # use anonymous channel parameters 
                    this_channel_params_dict = d  
        
        # search for channel-specific parameters         
        elif type(channel_params_dict) == dict:
            if 'channel' in channel_params_dict:
                if channel_params_dict['channel'] == channel:
                    this_channel_params_dict = channel_params_dict
            else:
                # use anonymous channel parameters
                this_channel_params_dict = channel_params_dict
                
        # if there is no dict for an anonymous channel, use default parameters
        if this_channel_params_dict is None:
            this_channel_params_dict = {}
    
        # do analysis for one frame
        df = do_one_frame(filename,frame, channel, \
                          img_info, load_file_path=filepath, \
                          verbose = True, **this_channel_params_dict)   

        # append output to list
        df_list.append(df)
    df_all = pd.concat(df_list)
    
    # generate file name and save `df_all` to csv file
    if not save_output_path is None:
        if not os.path.exists(save_output_path):
                os.makedir(save_output_path)
        fits_df.to_csv(os.path.join(save_output_path,'Combined_dotFitsDict.csv'))
        
    return df_all

def batch_parameter_sweep(function_to_call=None,**batch_arguments):
    """
    Makes `params_dict` dictionaries for all combinations of the parameters passed 
    to `batch_arguments`. Subsequently, calls `function_to_call` usings `params_dict`.  
    
    Usage
    ------
    E.g. The following will run "do_first_and_last()" for 'file1' and 'file2' for all 
         combinations of 'dog_sigma2' and 'dog_sigma1'.
    
        batch_parameter_sweep(do_first_and_last, filename='['file1','file2'], \
                              dog_sigma2=[3,15], dog_sigma1=[1.5,2])
    
    Parameters
    ----------
    function_to_call : function
                    Any function can use `params_dict` as an input argument
                    such as: do_all_frames_all_channels()
                    
    **batch_arguments : dict
                    Dictionary of keyword, value pairs. Keywords should be arguments for
                    `channel_params_dict` and the corresponding values may be iterables. 
    
    Returns
    -------
    params_dict_list : list of dicts
                    List of `params_dict` dictionaries - one dictionary for each permutation
                    of values passed to `batch_arguments`
    
    function_to_call_output : obj
                    The output object will depend on whatever is the output of `function_to_call`
    
    """ 
    
    # if filepath is specified, return error
    if 'filepath' in batch_arguments.keys():
        raise Exception('Do not specify ''filepath'' as an argument. Include the file path in filename.')
        
    special_keys = ['filename','channel']    
    
    # if any "value" in batch_arguments is not iterable, make it iterable
    # do not treat strings as iterable objects -> put non-iterable items into a list
    for key, value in batch_arguments.items():
        if not isinstance(value, Iterable):
            batch_arguments[key] = [value]
        elif isinstance(value,str):
            batch_arguments[key] = [value]
    
    # get all keyword : (iterable) value pairs that are not 'filename' and 'filepath'
    batch_dict = {x: batch_arguments[x] for x in batch_arguments if x not in special_keys}
    
    # get all keyword : (iterable) value pairs that are not 'filename' and 'channel'
    names_dict = {x: batch_arguments[x] for x in batch_arguments if x in special_keys}
    
    # generate all permutations of keyword : value pairs from the iterables in batch_dict
    channel_params_dict_list = [dict(zip(batch_dict.keys(), a)) for a in product(*batch_dict.values())]
    
    
    tmp_list = []
    if 'channel' in names_dict:
        for di, dict_item in enumerate(channel_params_dict_list):
            tmp_list.append([])
            # copy each channel_params_dict in channel_params_dict_list           
            # now add 'channel' information
            for ci, ch in enumerate(names_dict['channel']):
                tmp_list[di].append(copy.deepcopy(dict_item))
                tmp_list[di][ci].update({'channel':ch})
        channel_params_dict_list = tmp_list
    
    # create a list of params_dict
    params_dict_list = []
    for fullname in names_dict['filename']:
         for channel_params in channel_params_dict_list:
                
            filepath, filename  = os.path.split(fullname)
            params_dict_list.append({'filename':filename, \
                                 'filepath':filepath, \
                                 'channel_params_dict': channel_params})            
    
    function_to_call_output = None
    if not function_to_call is None:
        for params_dict in params_dict_list:
            function_to_call_output = function_to_call(params_dict)
            
    return params_dict_list


def do_first_and_last(plot_order=('frame','channel'),these_channels_only=None,\
                      figure_save_path=None, verbose=False, **params_dict):
    """
    Performs "do_one_frame()" on the first and last frame of the data series for the channels 
    specified by `channel_params_dict` and outputs intermediates for easy visualization. If no 
    channels are specified, all channels are used
    
    
    Parameters
    ----------
    plot_order : ('frame','channel') or ('channel', 'frame')
               Default value ('frame','channel') plots the 'first' and 'last' frame, grouped by 'channel'. 
               The value ('channel','frame') plots the channels in order, grouped by frame. 
    
    these_channels_only : list of ints, optional
                The list of integers in `these_channels_only` specify the specific channels for which 
                ouput plots are desired. Plots are generated for the first and last frame.
                If `these_channels_only` is unspecified, plots are generated for all channels.
                
    figure_save_path : str, optional
                Path to a directory where the ouput plots will be saved to disk. If the directory
                does not exist, it is created. If `figure_save_path` is not set, plots are not saved. 
                
    verbose : {T,F}, optional
                Toggles 'on' or 'off' the verbose option of do_one_frame(). See do_one_frame()
                
    **params_dict : dict or keyword/value pairs, containing 3 mandatory components
                1) filename : str
                              Name of the image series to be analyzed
                2) filepath : str
                              Full file path to the image series
                3) channel_params_dict : dict (or list of dicts)
                              A dictionary of keyword arguments that will be passed to do_one_frame().
                              If `channel` is anonymous (i.e. if it is not specified in the arguments), 
                              the parameters passed  to do_one_frame() will be the same for all channels. 
                              If `channel_params_dict` contains a list of dictionaries (e.g. specifying 
                              parameters for each channel), the channel-specific parameters will be passed
                              to do_one_frame(); if channel-specific parameters are unspecified, the 
                              function attempts to use parameters from the last anonymous channel. Otherwise,
                              function defaults are used.
                      
                If a dict is passed (instead of keyword/value pairs), `params_dict` values must be unpacked 
                (i.e. pass `**params_dict` to the function instad of simply `params_dict`).
                
    Example formatting
    ------------------
    E.g. params_dict = {'filename': 'my_file.czi',
                        'filepath': './file_location/',
                        'channel_params_dict': [{'channel': 0,
                                              'dog_sigma1': 1.5, 
                                              'dog_sigma2': 15, 
                                              'min_object_area': 50, 
                                              'intensity_threshold_method': 'percentile', 
                                              'percentile_threshold': 99,
                                              'window_size': 10,
                                              'save_intermediates_file_path': './tmp'},

                                               {'channel': 1,
                                               'dog_sigma1': 1.5, 
                                               'dog_sigma2': 3, 
                                               'min_object_area': 35, 
                                               'intensity_threshold_method': 'percentile', 
                                               'percentile_threshold': 99,
                                               'window_size': 10,
                                               'save_intermediates_file_path': './tmp'},
                                                ]
                         }   
    See also: batch_parameter_sweep() as a method to easily generate this dict.
    
                
    Output
    -------
    
    Figure 1: {filename}_LocaDotPlots.pdf
                Shows the maximum intensity projection of onto the xy plane of 1) the raw z-stack, 
                2) the filtered z-stack, 3) The dot segmentation (with 'dot' IDs).
                
    Figure 2 and up : {filename}_DotsFitPlots_Channel{channel}_Frame_{frame}.pdf
                Each figure will show the maximum projection of a small window around each dot.
                For each dot, xy, xz and yz projections are shown for 1) the raw data and 
                2) the PSF fits. 
    
    """ 

    filename = params_dict['filename']
    filepath = params_dict['filepath']

    # get metadata
    img_info, _ = get_CZI_metadata(filename,filepath)
    (sizeX,sizeY,sizeZ,sizeT,num_channels) = img_info
    
    
    # prepare figure for plotting 
    numPlots = len(list(product(range(num_channels),[0,sizeT-1])))
    xwidth = 3
    ywidth = numPlots
    width_inches = 3

    fig, gs = _gridspec_inches(wcols = np.array([width_inches]*xwidth),\
                              hrows =np.array([width_inches]*ywidth),\
                             hspace=0.35,wspace=0.25)

    # iterate over channels and frames
    count = 0
    
    if these_channels_only is 'None':
        constrain_to_channels = range(num_channels)
    else:
        constrain_to_channels = these_channels_only
        
    if plot_order == ('channel','frame'):
        channel_frame_pairs = product(constrain_to_channels,[0,sizeT-1])
    else:
        channel_frame_pairs = [(y,x) for (x,y) in list(product([0,sizeT-1],constrain_to_channels))]
    
    for (channel,frame) in channel_frame_pairs:

        # retrieve fitting parameters for specified channel
        channel_params_dict = params_dict['channel_params_dict']
        this_channel_params_dict = None
        
        # search list for channel-specific parameters dicts
        if type(channel_params_dict)==list:
            for d in channel_params_dict:
                if 'channel' in d:
                    if d['channel'] == channel:
                        # use channel-specific parameters
                        this_channel_params_dict = d
                        break
                else:
                    # use anonymous channel parameters 
                    this_channel_params_dict = d  
        
        # search for channel-specific parameters         
        elif type(channel_params_dict) == dict:
            if 'channel' in channel_params_dict:
                if channel_params_dict['channel'] == channel:
                    this_channel_params_dict = channel_params_dict
            else:
                # use anonymous channel parameters
                this_channel_params_dict = channel_params_dict
                
        # if there is no dict for an anonymous channel, use default parameters
        if this_channel_params_dict is None:
            this_channel_params_dict = {}

        # do analysis for one frame
        fits_df, dot_fits_dict, zstack, filtered_zstack, centroids, blobs_labels, _ = \
                        do_one_frame(filename,frame, img_info=img_info, load_file_path=filepath, \
                          return_intermediates=True,verbose = verbose, **this_channel_params_dict)   
        
        # plot maximum projection of the z-stack
        max_proj = np.max(zstack,2)
        filt_max_proj = np.max(filtered_zstack,2)

        # save segmentation of images
        plt.figure(fig.number)
        plt.subplot(gs[count]); count += 1;
        plt.imshow(max_proj, vmin=np.percentile(max_proj,2),\
                   vmax=np.percentile(max_proj,99.5),cmap='coolwarm')
        plt.title('Channel {} Frame {} (Raw)'.format(channel,frame))
        plt.subplot(gs[count]); count += 1;
        plt.imshow(filt_max_proj, vmin=np.percentile(filt_max_proj,2),\
                   vmax=np.percentile(filt_max_proj,99.5),cmap='coolwarm')
        plt.title('(Filtered)'.format(channel,frame))
        plt.subplot(gs[count]); count += 1;
        plt.imshow(blobs_labels>0,vmax=1,cmap='gray')
        plt.title('(Localization)'.format(channel,frame))
        # superimpose the dots on the segmented image
        for dot_key, dot_dict in dot_fits_dict.items():
            plt.text(dot_dict['y0_fit'], dot_dict['x0_fit'],"{}".format(dot_key),color='y')
            #plt.plot(dot_dict['y0_fit'], dot_dict['x0_fit'],'o',markersize=3)
        
        ## show dots and fits
        numDots = len(dot_fits_dict)
        dot_fits_keys = [key for key in dot_fits_dict[0].keys() if 'max' in key]
        numKeys = len(dot_fits_keys)
        xwidth = 3
        ywidth = int(numKeys*numDots/3)

        fig_dots, gs_dots = _gridspec_inches(wcols = np.array([width_inches]*xwidth),\
                          hrows =np.array([width_inches]*ywidth),\
                         hspace=0.25,wspace=0.25)

        plt.figure(fig_dots.number)
        dot_subplot_count = 0
        for dot_id in dot_fits_dict.keys():
            dot_dict = dot_fits_dict[dot_id]
            for key in dot_fits_keys:
                plt.subplot(gs_dots[dot_subplot_count])
                if 'z' in key:
                    plt.imshow(dot_dict[key].T)
                else:
                    plt.imshow(dot_dict[key])
                if key in dot_dict['num_modes']:
                    plt.title("{}\n# dots found: {}".format(key,dot_dict['num_modes'][key]))
                else: 
                    plt.title(key)
                if np.mod(dot_subplot_count,3)==0:
                    plt.ylabel("Dot ID: {}".format(dot_id))
                dot_subplot_count += 1


        if not figure_save_path is None:
            # make directory if it does not yet exist
            if not os.path.exists(figure_save_path):
                os.mkdir(figure_save_path)
            # save figure
            plt.figure(fig_dots.number)
            dot_figure_name = "{}_DotsFitPlots_Channel{}_Frame_{}.pdf" \
                                .format(filename[:-4],channel,frame)
            plt.savefig(os.path.join(figure_save_path,dot_figure_name), bbox_inches = "tight" ) 
            
    if not figure_save_path is None:
        # make directory if it does not yet exist
        if not os.path.exists(figure_save_path):
            os.mkdir(figure_save_path)
        # save figure
        plt.figure(fig.number) # call the correct figure
        dot_figure_name = "{}_LocaDotPlots.pdf" \
                            .format(filename[:-4],channel,frame)
        plt.savefig(os.path.join(figure_save_path,dot_figure_name), bbox_inches = "tight" )             
            
    return fits_df

# internal helper function to help with plotting
def _gridspec_inches(
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


def count_dots_from_threshold(img,threshold_percentage=98,\
                                 min_object_area = 2,return_segmentation=False):
    """
    Segments images based on a simple threshold and returns the number of observed spots
    
    Parameters
    ----------
    
    img : numpy.ndarray
            Image on which to perform segmentation
            
    threshold_percentage : float
            Percentage threshold on the `img` intensity values used to binarize the image
            
    min_object_area : int
            Minimum number of pixels used to call/identify object
            
    return_segmentation : {T,F}
            If True, the function returns 1) the number of dots, 2) metrics about each dot.
            If False, the function only returns the number of dots
            
    Returns
    -------
    
    1) len(blobs_metrics) : int
            Number of dots identified from the image segmentation
            
    2) `blobs_metrics` : object
            Object from skimage.measure.regionprops() applied on the thresholded `img`.
    """
    
    blobs=  img>np.percentile(img,threshold_percentage) 

    # filter objects based on size    
    blobs = remove_small_objects(blobs, min_size=min_object_area) 

    # remove objects touching the edges of the image
    blobs = clear_border(blobs) 

    # get segmentation of the image from connected pixels 
    blobs_labels = measure.label(blobs, background=0) 

    # measure things for each unique feature identified in blobs_labels
    blob_metrics = measure.regionprops(blobs_labels, img)   
    
    if return_segmentation == True:
        return len(blob_metrics), blobs
    else:
        return len(blob_metrics)
    
    
    
def log_zerocross(img,log_sigma=2,num_std_threshold=0):
    """
    Returns the edges detected by the Laplacian of Gaussian method applied on the image.
    
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


def _get_weighted_centroid(zstack,xc,yc, win,dot_intensity_percentile,\
                           img_info, min_dot_size = 50, max_dot_size=500):
    """
    Parameters
    ----------
    (see )
    
    Output
    ------
    centerX, centerY, centerZ : float
    
    mean_intensity : float
    
    mean_surrounding : float
    
    sum_intensity :float
    
    dot_size_pixels : int
    """    
    sizeX,sizeY,sizeZ,sizeT,num_channels = img_info
    # crop out the "dot" from the zstack
    dot_volume =  zstack[:,yc-win:yc+win,xc-win:xc+win]            
    binary_volume = np.array((dot_volume > 
                              np.percentile(dot_volume,dot_intensity_percentile))
                             ,dtype=int)


    # create a mask for the dot using connected pixels components on binarized volume
    blobs, num_features = scipy.ndimage.measurements.label(binary_volume)

    # filter out for size
    freq = np.bincount(np.ndarray.flatten(blobs))
    freq[freq>max_dot_size] = 0 
    freq[freq<min_dot_size] = 0
    dot_label = np.argmax(freq)    
    mask = blobs
    mask[mask != dot_label] = 0

    # calculate the centroid position
    weights = np.ndarray.flatten(mask*dot_volume)
    if np.sum(weights) == 0: 
        return None # (if there are no dots, skip)
    else:
        # prepare meshgrid for finding centroid positions of the dots
        zz, yy, xx = np.meshgrid(range(zstack.shape[0]), 
                                 np.arange(yc-win,yc+win), 
                                 np.arange(xc-win,xc+win),indexing='ij')  
        # calculate the centroid position
        centerX = np.average(np.ndarray.flatten(xx),weights=weights)
        centerY = np.average(np.ndarray.flatten(yy),weights=weights)
        centerZ = np.average(np.ndarray.flatten(zz),weights=weights)                

    # calculate the total intensity/ mean intensity of the dot, and surrounding background
    mean_intensity = np.nanmean(dot_volume[mask!=0])
    mean_surrounding = np.nanmean(dot_volume[mask==0])
    sum_intensity = np.nansum(dot_volume[mask!=0])
    dot_size_pixels = np.count_nonzero(np.ndarray.flatten(mask))
    
    return centerX, centerY, centerZ, mean_intensity, mean_surrounding, sum_intensity, dot_size_pixels

def get_weighted_centroid_from_dot_volume(dot_volume,x_low,y_low,z_low, win,dot_intensity_percentile,\
                           img_info, min_dot_size = 50, max_dot_size=500):
    """
    Parameters
    ----------
    
    xc, yc, zc : float
        Real dot position in the Z-stack (not relative position)
    
    Output
    ------
    """    
    sizeX,sizeY,sizeZ,sizeT,num_channels = img_info # metadata
    
    winZ, winY, winX = win # window size around localized dot for weighted average
    
    # crop out the "dot" from the zstack & render binary         
    binary_volume = np.array((dot_volume > 
                              np.percentile(dot_volume,dot_intensity_percentile))
                             ,dtype=int)

    # create a mask for the dot using connected pixels components on binarized volume
    blobs, num_features = scipy.ndimage.measurements.label(binary_volume)

    # filter out dots based on size
    freq = np.bincount(np.ndarray.flatten(blobs))
    freq[freq>max_dot_size] = 0 
    freq[freq<min_dot_size] = 0
    dot_label = np.argmax(freq)    
    mask = blobs
    mask[mask != dot_label] = 0

    # calculate the centroid position
    weights = np.ndarray.flatten(mask*dot_volume)
    if np.sum(weights) == 0: 
        return None # (if there are no dots, skip)
    else:
        # prepare meshgrid for finding centroid positions of the dots
        zz, yy, xx = np.meshgrid(np.arange(z_low,z_low+winZ), 
                                 np.arange(y_low,y_low+winY), 
                                 np.arange(x_low,x_low+winX),indexing='ij')  
        # calculate the centroid position
        centerX = np.average(np.ndarray.flatten(xx),weights=weights)
        centerY = np.average(np.ndarray.flatten(yy),weights=weights)
        centerZ = np.average(np.ndarray.flatten(zz),weights=weights)                

    # calculate the total intensity/ mean intensity of the dot, and surrounding background
    mean_intensity = np.nanmean(dot_volume[mask!=0])
    mean_surrounding = np.nanmean(dot_volume[mask==0])
    sum_intensity = np.nansum(dot_volume[mask!=0])
    dot_size_pixels = np.count_nonzero(np.ndarray.flatten(mask))
    
    return (centerX, centerY, centerZ, mean_intensity, mean_surrounding, sum_intensity, dot_size_pixels)


import re
import itertools
def natural_sort(l):
    # found here: https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)
