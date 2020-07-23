from setuptools import setup
#from distutils.core import setup

setup(name="connect_the_dots",
			version='0.1',
			description='Particle tracking software common',
			author='Hugo B. Brandao',
			author_email='hbrandao@g.harvard.edu',
			license='MIT',
			packages=['connect_the_dots'])


"""
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
"""