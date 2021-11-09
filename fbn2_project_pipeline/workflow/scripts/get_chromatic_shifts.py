import os
import pandas as pd
import datetime
import numpy as np
import trackpy
import matplotlib.pyplot as plt
#from connect_the_dots.io import get_CZI_metadata
from scipy.signal import savgol_filter
trackpy.quiet()

##################################################################################
# PARAMETERS 
##################################################################################
num_dates_to_average = 7 # the number of bead files to average over

# some parameters - need to put into configuration file eventually
verbose = False
bead_intensity_min = 100 # minimum volume in pixel counts
bead_intensity_max = 2000 # maximum volume in pixel counts
search_range = [3,4,4] # for linking beads across channels
min_len = 5 # minimum trajectory length for beads
memory = 0 # number of allowed missing frames for bead trajectories
    
sav_gol_window = 91 # filtering window length in X and Y
sav_gol_windowZ = 21 # filtering window length in Z
sav_gol_order = 1 # order of SG filter
step_size = 0.5 # step size in pixels for saving the chromatic shift


# superblocks -- need to get from configuration file
super_blocks = [('2020-08-27','2020-10-19'), # Colibri 7 installed
                ('2020-10-20','2020-11-10'), # Colibri 7 to 5 swap
                ('2020-11-11','2021-02-07'), # Definite focus issue fixed
                ('2021-02-08','2021-05-06'), # Piezo adjusted
                ('2021-05-07','2021-12-31'), # PSF check date
               ]

##################################################################################

# GET METADATA 
# get pixel sizes and image dimensions from sample metadta
# data_info , sample_metadata = get_CZI_metadata(snakemake.input.sample_czi)
fx = 0.08507776712474871 #sample_metadata['OME']['Image']['Pixels']['PhysicalSizeX']
fy = 0.08507776712474871 #sample_metadata['OME']['Image']['Pixels']['PhysicalSizeY']
fz = 0.16 #sample_metadata['OME']['Image']['Pixels']['PhysicalSizeZ']
ux = 'um' #sample_metadata['OME']['Image']['Pixels']['PhysicalSizeXUnit']
uy = 'um' #sample_metadata['OME']['Image']['Pixels']['PhysicalSizeYUnit']
uz = 'um' #sample_metadata['OME']['Image']['Pixels']['PhysicalSizeZUnit']
(sizeX,sizeY,sizeZ,sizeT,num_channels) = (584, 584, 30, 365, 2) #data_info

traj = snakemake.input.traj

# get date from the sample movie
traj = traj.split('/')[-1]
try:
    YY = int(traj.split('_')[0])
    MM = int(traj.split('_')[1])
    DD = int(traj.split('_')[2])  
except:
    YY = int(traj[0:4])
    MM = int(traj[5:7])
    DD = int(traj[8:10])    
sample_date =  datetime.datetime(YY,MM,DD)

# get dates from the beads movies
beads_files = snakemake.input.bead_localizations
        
# get dates from the file name (assumes YY_MM_DD)
date_list = []
for f in beads_files:
    f = f.split('/')[-1]
    try:
        YY = int(f.split('_')[0])
        MM = int(f.split('_')[1])
        DD = int(f.split('_')[2])
        date = datetime.date(YY, MM, DD)
    except:
        try:
            YY = int(f[0:4])
            MM = int(f[5:7])
            DD = int(f[8:10])
            date = datetime.date(YY, MM, DD)                        
        except:        
            date = ''
    date_list.append(date)
    
# create datetime object for the dates
df = pd.DataFrame({'date':date_list})
df['date'] = pd.to_datetime(df['date'])  
df['beads_files'] = beads_files


# remove rows that have the NaT timestamp
df = df.loc[df.date.notnull()]

# get the nearest date
def get_nearest_date(bead_datetime_list, sample_datetime):
    return min(bead_datetime_list, key=lambda x: abs(x - sample_datetime))

def dates_to_average(idx,map_number_to_date,group_size=3,verbose=False):
    
    try:
        assert group_size%2 == 1
    except:
        print("group_size must be an odd number")
        return 0
        
    if group_size > len(map_number_to_date):
        group_size = len(map_number_to_date)
        
    if idx >= len(map_number_to_date):
        print(f"Date out of range for idx={idx}")
        return 0
    
    
    inds = np.array(sorted(map_number_to_date.keys()))
    lower = inds > idx+group_size//2
    upper = inds < idx-group_size//2

    # for non-edge cases
    if (idx >= group_size//2) and (idx < len(inds)-group_size//2):
        indices = [i for i, (l, u) in enumerate(zip(lower, upper)) if (l==False) and (u==False)]
        if verbose==True:
            print("Middle")

    elif (idx >= len(inds)-group_size//2):
        indices = [i for i in np.arange(len(inds)-group_size,len(inds))]
        if verbose==True:
            print('Lower')    
    else:
        indices = [i for i in np.arange(0,group_size)]
        if verbose==True:
            print("Upper")  
    
    dates = [map_number_to_date[inds[i]] for i in indices]
    return dates


# get the superblock - filter out the bead dates outside the given superblock
block_df = []
for start_date, end_date in super_blocks:
    YY = int(start_date.split('-')[0])
    MM = int(start_date.split('-')[1])
    DD = int(start_date.split('-')[2])        
    start_datetime = datetime.datetime(YY,MM,DD)
    
    YY = int(end_date.split('-')[0])
    MM = int(end_date.split('-')[1])
    DD = int(end_date.split('-')[2])       
    end_datetime = datetime.datetime(YY,MM,DD)
        
    # is this the correct date superblock?
    if start_datetime <= sample_date <= end_datetime:
        
        # create filtered dataframe with only dates within sample superblock 
        block_df = df[(df['date'] >= start_date) & (df['date'] <= end_date)]          
        break
      
# get the list of unique dates, and sort
sorted_unique_dates = sorted(set(block_df['date']))

# create mappings from date to a running index for grouping
map_date_to_number = {}
map_number_to_date = {}
for di, d in enumerate(sorted_unique_dates):
    map_date_to_number[d] = di
    map_number_to_date[di] = d
                                    
# get the bead date that is closest to the sample date
nearest_bead_date = get_nearest_date(block_df.date, sample_date)
nearest_bead_date_idx = map_date_to_number[nearest_bead_date]



bead_dates_group = dates_to_average(nearest_bead_date_idx,
                                    map_number_to_date,
                                    group_size=num_dates_to_average,
                                    verbose=False)


files_to_average = block_df[block_df['date'].isin(bead_dates_group)].beads_files

if verbose == True:
    print(f"Averaging over {len(files_to_average)} files")
if len(files_to_average)==0:
    if verbose == True:
        print(f"Skipping this date range -- no files found.")
    assert(len(files_to_average)>0)


## get the chromatic shifts for from the selected beads files
DX = []
DY = []
DZ = []
X1 = []
Y1 = []
Z1 = []
X0 = []
Y0 = []
Z0 = []
for file in files_to_average.values:

    # load pre-processed bead locations
    df_particles = pd.read_csv(file)

    # filter based on bead intensity
    df_particles = df_particles.query(f"dot_size_in_pixels>{bead_intensity_min}" +
                                      f"& dot_size_in_pixels<{bead_intensity_max}")

    colour_scheme = ['r','b','g','m','c']
    windowZYXdims = search_range

    # load csv files with segmentation positions ouputs
    df_particles['timestep'] = df_particles['frame']


    # get special positions
    c = df_particles['channel'].values
    t = df_particles['timestep'].values
    x = df_particles['x'].values
    y = df_particles['y'].values
    z = df_particles['z'].values

    # make new, reduced data frame
    df_short = pd.DataFrame({'channel': c, 'frame': t,'x': x, 'y': y, 'z': z})

    # do for each channel
    channels = set(df_short['channel'].values)

    # link trajectories across channels:
    # artificially shift the frames and create trajectories
    linked_df = df_short.copy()
    linked_df['frame'] = linked_df['frame']*len(channels)+linked_df['channel']
    try:
        linked_df = trackpy.link_df(linked_df, \
                    search_range=search_range, \
                    memory=memory,)
    except:
        continue
    # reverse shift the frame numbers
    linked_df['frame'] = (linked_df['frame']-linked_df['channel'])/len(channels)

    # get unique IDs for sets of linked dots
    particle_set = set(linked_df.particle)

    for p in particle_set:

        ch = 0
        x_vals0 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].x.values
        y_vals0 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].y.values
        z_vals0 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].z.values
        t_vals0 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].frame.values

        ch = 1
        x_vals1 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].x.values
        y_vals1 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].y.values
        z_vals1 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].z.values
        t_vals1 = linked_df[(linked_df.particle==p) & (linked_df.channel==ch)].frame.values

        if (len(x_vals0) == len(x_vals1)) and (len(x_vals0) > min_len) and np.all(t_vals0==t_vals1):
            X0.extend(x_vals0)
            Y0.extend(y_vals0)
            Z0.extend(z_vals0)
            X1.extend(x_vals1)
            Y1.extend(y_vals1)
            Z1.extend(z_vals1)
            DX.extend( x_vals1-x_vals0)
            DY.extend( y_vals1-y_vals0)
            DZ.extend( z_vals1-z_vals0)

X0 = np.array(X0)
Y0 = np.array(Y0)
Z0 = np.array(Z0)
X1 = np.array(X1)
Y1 = np.array(Y1)
Z1 = np.array(Z1)
dx = np.array(DX)
dy = np.array(DY)
dz = np.array(DZ)


# add the trend lines
maxX = sizeX 
maxY = sizeY 
maxZ = sizeZ 
minX = 0
minY = 0
minZ = 0


# do X  axis
xq = np.arange(minX,maxX,step_size)
xx = np.array(X0)
yy = np.array(dx)
ss = np.argsort(xx)
yq = np.interp(xq,xx[ss],yy[ss])
yq = savgol_filter(yq,sav_gol_window,sav_gol_order)
chromatic_shift_X = np.c_[xq,yq]

# do Y axis
xq = np.arange(minY,maxY,step_size)
xx = np.array(Y0)
yy = np.array(dy)
ss = np.argsort(xx)
yq = np.interp(xq,xx[ss],yy[ss])
yq = savgol_filter(yq,sav_gol_window,sav_gol_order)
chromatic_shift_Y = np.c_[xq,yq]

# do Z axis
xq = np.arange(minZ,maxZ,step_size)
xx = np.array(Z0)
yy = np.array(dz)
ss = np.argsort(xx)
yq = np.interp(xq,xx[ss],yy[ss])
yq = savgol_filter(yq,sav_gol_windowZ,sav_gol_order)
chromatic_shift_Z = np.c_[xq,yq]

if verbose == True:
    print("Saving chromatic shift")
# save the chromatic shift information
np.savez(snakemake.output.chromatic_shift,
         chromatic_shift_X=chromatic_shift_X,
         chromatic_shift_Y=chromatic_shift_Y,
         chromatic_shift_Z=chromatic_shift_Z,
         fx=fx,fy=fy,fz=fz,ux=ux,uy=uy,uz=uz)
