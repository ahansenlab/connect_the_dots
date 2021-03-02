## Filtering functions for 3D or 4D time series data. 
# Written by Hugo B. Brandao (hbrandao@g.harvard.edu)
# (c) 2021, Hugo B. Brandao

## Import required modules
import numpy as np # for manipulating arrays
from scipy.optimize import curve_fit # for making fits to the PSF
from scipy.ndimage import gaussian_laplace, gaussian_filter # for dot localization (image filtering)
from skimage import measure # for segmenting images
from skimage.morphology import remove_small_objects, closing, disk, dilation # for morphological filtering of images
from skimage.segmentation import clear_border # for filtering images
from skimage.filters import threshold_otsu
import pandas as pd # for creating and manipulating tabulated data
from collections import Iterable
import itertools
from itertools import product, groupby
import copy
import scipy
import trackpy # library containing tracking algorithms
import pywt # wavelet transform library
import re # regex 
import warnings

from pathlib import Path


def natural_sort(l):
    """
    Takes in a list of strings and returns the list sorted in "natural" order. 
    (e.g. [test1, test10, test11, test2, test20] -> [test1, test2, test10, test11, test20])
    Source: https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
    
    Parameters
    ----------
    l : list of str
        Unsorted list of strings
        
    Returns
    -------
    sorted_l : list of str
        Sorted list of strings
   
    """    
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


def filter_intensity_jumps(i0,winsz=20,num_quantiles=2,quantile=0.5):
    """
    Identify locations in the trajectory where the intensity jumps significantly
    
    Parameters
    ----------
    i0 : numpy.ndarray
        Dot intensity time series.
        
    winsz : int
        Window size for creating the trendline
        
    num_std : float
        Number of `quantile` deviations used for spotting outliers
        
    quantile : float
        Value from 0 to 1 used to get a cutoff value for calling outliers.
        Quantile is computed on the absolute difference of `i0` and the trend line.
        It is a rolling quantile of size `winsz`.
        
    Results
    -------
    good_points : list of booleans
        List indicating which points of `i0` are expected and which may be outliers.
    
    """
    trend = np.convolve(i0,np.ones(winsz)/winsz)[0:len(i0)]
    resid = pd.Series(np.abs(i0-trend))
    trend_std = resid.rolling(winsz).quantile(quantile)
    upper_error = i0 > (trend+num_quantiles*trend_std)
    lower_error = i0 < (trend-num_quantiles*trend_std)
    good_points = (upper_error | lower_error)==False
    return good_points

def create_mask(dot_volume,xc,yc,zc,radius_xy=4,radius_z=3):
    """
    Creates a 3D ball-shaped mask at the specified pixel coordinate.
    
    Parameters
    ----------
    dot_volume : numpy.ndarray
        3D voxel around a dot.
    
    xc, yc, zc : float
        Estimated location of the dot within the voxel.
        
    radius_xy : int
        Radius (in units of pixels) used to create the mask for the X,Y dimensions
        
    radius_z : int
        Radius (in units of pixels) used to create the mask for the Z dimension
        
    
    Returns
    -------
    mask : numpy.ndarray
        Mask used for calculating the intensity weighted centroid.
    
    """    
    mask = np.zeros_like(dot_volume)    
    for i in range(mask.shape[0]):
        if np.sqrt((i-zc)**2) > radius_z:
            continue
 
        for j in range(mask.shape[1]):
            if np.sqrt((j-yc)**2) > radius_xy:
                continue
        
            for k in range(mask.shape[2]):
                if np.sqrt((k-xc)**2+(j-yc)**2) <= radius_xy:
                    mask[i,j,k] = 1                
    return mask

def create_mask_symmetric(dot_volume,xc,yc,zc,radius_xy=4,radius_z=3,min_rxy=2,min_rz=1):
    """
    Creates a 3D ball-shaped mask at the specified pixel coordinates.
    The mask radius near the boundaries are adjusted to maintain symmetry. 
    
    Parameters
    ----------
    dot_volume : numpy.ndarray
        3D voxel around a dot.
    
    xc, yc, zc : float
        Estimated location of the dot within the voxel.
        
    radius_xy : int
        Radius (in units of pixels) used to create the mask for the X,Y dimensions
        
    radius_z : int
        Radius (in units of pixels) used to create the mask for the Z dimension
                
    min_rxy : int
        Minimum radius size for the adaptive masking in X or Y
        
    min_rz : int
        Minimum radius size for the adaptive masking in Z. 
    
    Returns
    -------
    mask : numpy.ndarray
        Mask used for calculating the intensity weighted centroid.
    
    """    
    mask = np.zeros_like(dot_volume)    
    
    # check boundary case for z
    dist_edgeZ = np.min([np.abs(mask.shape[0]-zc-1),zc])
    if dist_edgeZ <= radius_z:
        if dist_edgeZ >= min_rz:
            radius_z = dist_edgeZ
        else:
            return mask

    # check boundary case for y
    dist_edgeY = np.min([np.abs(mask.shape[1]-yc-1),yc])
    if dist_edgeY <= radius_xy:
        if dist_edgeY >= min_rxy:
            radius_xy = dist_edgeY
        else:
            return mask        
        
    # check boundary case for x
    dist_edgeX = np.min([np.abs(mask.shape[2]-xc-1),xc])
    if dist_edgeX <= radius_xy:
        if dist_edgeX >= min_rxy:
            radius_xy = dist_edgeX
        else:
            return mask
        
    for i in range(mask.shape[0]):
        if np.sqrt((i-zc)**2) > radius_z:
            continue
 
        for j in range(mask.shape[1]):
            if np.sqrt((j-yc)**2) > radius_xy:
                continue
        
            for k in range(mask.shape[2]):
                if np.sqrt((k-xc)**2+(j-yc)**2) <= radius_xy:
                    mask[i,j,k] = 1                
    return mask