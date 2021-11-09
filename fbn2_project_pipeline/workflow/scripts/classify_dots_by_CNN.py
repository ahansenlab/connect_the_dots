## Import required modules
import numpy as np 
import pandas as pd
import pickle
from tqdm.keras import TqdmCallback
import keras
from connect_the_dots.filtering import create_mask_symmetric

verbose = True
if verbose == True:
    print(f"Doing file {snakemake.output.traj_ml}")

    
def data_to_keras_input(dot_volume_ch0, centroids_ch0, dot_volume_ch1, centroids_ch1,data_shape=(12,12,8)):
    new_vol = np.zeros((data_shape[0],data_shape[1],data_shape[2],4))
    z, x, y = centroids_ch0[0]
    
    for vi, (vol, loc) in enumerate([(dot_volume_ch0,centroids_ch0),(dot_volume_ch1,centroids_ch1)]):        
        # prepare the raw voxel
        z, x, y = loc[0]
        raw_vol = vol
        lenZ, lenX, lenY  = raw_vol.shape
        vol = np.moveaxis(vol, 0, 2)
        vol_min = np.min(vol)
        vol_max = np.max(vol)
        raw_vox = (vol-vol_min)/(vol_max-vol_min)*2 - 1
        new_vol[:lenX, :lenY, :lenZ, vi] = raw_vox
        # prepare the localization voxel
        mask = create_mask_symmetric(raw_vol,x,y,z)
        mask = np.moveaxis(mask, 0, 2)
        new_vol[:lenX, :lenY, :lenZ, 2 + vi] = raw_vox*mask

    return np.array([new_vol]), np.c_[centroids_ch0,centroids_ch1]/np.r_[data_shape,data_shape]

# load models for the CNN classification
replication_model = keras.models.load_model(snakemake.input.ml_replication_model)
good_bad_model = keras.models.load_model(snakemake.input.ml_good_bad_model)



# Load data
S = pickle.load(open(snakemake.input.dot_volumes,'rb')) # load dot volume
S_df = pd.read_csv(snakemake.input.traj) # load trajectory 
S_df['ML_class'] = ['Unclassified']*len(S_df) # new column for classification of Good/Bad dots
S_df['ML_isReplicated'] = ['Unclassified']*len(S_df) # new column for classification of Replication or not


# for each particle, use CNN to classify the dot
for p in sorted(list(S.keys())): 
    keep_frame = []
    discard_frame = []
    replicated_frame = []
    unreplicated_frame = []        
    no_dot = []

    max_frame = max([max(list(sorted(S[p][c].keys()))) for c in list(sorted(S[p].keys()))])
    min_frame = min([min(list(sorted(S[p][c].keys()))) for c in list(sorted(S[p].keys()))])
    both_spots_localized = np.ones(max_frame+1, dtype=bool)
    both_spots_localized[:min_frame] = False

    channels = list(sorted(S[p].keys()))

    for c in channels :
        for frame in list(sorted(S[p][c].keys())):
            if np.any(np.isnan(S[p][c][frame][1])):
                both_spots_localized[frame] = False

    # get frames that exist in both channels
    frames_dict = {c:list(sorted(S[p][c].keys())) for c in channels} # get list of frames for each channel
    consensus_frames = sorted(list(set(frames_dict[0]).intersection(set(frames_dict[1])))) 

    count = 0
    for f, frame in enumerate(consensus_frames):
        
        try:
            dot_volume_ch0, centroids_ch0, centroids_ch0_true = S[p][0][frame]
            dot_volume_ch1, centroids_ch1, centroids_ch1_true = S[p][1][frame]
        except:
            print(f"Missing localization on frame {frame}")
            continue                                    

        is_good = True
        is_replicated = True
        if (type(centroids_ch0)==tuple or type(centroids_ch1)==tuple)==False:
            # prepare the dot volume for ML 
            try:
                volume_data, loc_data =  data_to_keras_input(dot_volume_ch0, centroids_ch0, dot_volume_ch1, centroids_ch1)
                is_good = good_bad_model.predict([volume_data, loc_data])[0][0]>0.5 # is good
                is_replicated = replication_model.predict([volume_data, loc_data])[0][0]>0.5 # is replicated
            except:
                continue

        # keep track of replication
        if is_replicated == True:
            replicated_frame.append(frame) 
        else:
            unreplicated_frame.append(frame) 

        # keep track of replication
        if (is_good == True) and (both_spots_localized[frame]==True):
            keep_frame.append(frame) 
        elif (is_good==False) and (both_spots_localized[frame]==True):
            discard_frame.append(frame)
        elif (both_spots_localized[frame]==False):
            no_dot.append(frame)

    S_df['ML_class'].mask((S_df.particle==p) & S_df.frame.isin(no_dot),'NoDot',inplace=True)
    S_df['ML_class'].mask((S_df.particle==p) & S_df.frame.isin(discard_frame),'Bad',inplace=True)
    S_df['ML_class'].mask((S_df.particle==p) & S_df.frame.isin(keep_frame),'Good',inplace=True)

    S_df['ML_isReplicated'].mask((S_df.particle==p) & S_df.frame.isin(replicated_frame),'Replicated',inplace=True)
    S_df['ML_isReplicated'].mask((S_df.particle==p) & S_df.frame.isin(unreplicated_frame),'Unreplicated',inplace=True)   

S_df.to_csv(snakemake.output.traj_ml,index=False)
