
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/hbrandao/anaconda3/lib/python3.7/site-packages', '/home/hbrandao/libs/data_analysis_Fbn2/snakemake/workflow/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X\x93\x00\x00\x00results/Rad21/2021_01_19_Fbn2_Rad21F1M_noIAA_488nm_0p5_561nm_0p05_SC8_2X_30z_250nm_780V_365T_20s_movie1-01__processed_low..czi.tracks_corrected.csvq\x06X\x8e\x00\x00\x00results/Rad21/2021_01_19_Fbn2_Rad21F1M_noIAA_488nm_0p5_561nm_0p05_SC8_2X_30z_250nm_780V_365T_20s_movie1-01__processed_low..czi.dot_volumes.pklq\x07X>\x00\x00\x00ml_models/Replication_model_2021_04_14_hugo_classifications.h5q\x08XH\x00\x00\x00ml_models/4D_model_2021_04_08_mixedInputClean_v2_doubleClassification.h5q\te}q\n(X\x06\x00\x00\x00_namesq\x0b}q\x0c(X\x04\x00\x00\x00trajq\rK\x00N\x86q\x0eX\x0b\x00\x00\x00dot_volumesq\x0fK\x01N\x86q\x10X\x14\x00\x00\x00ml_replication_modelq\x11K\x02N\x86q\x12X\x11\x00\x00\x00ml_good_bad_modelq\x13K\x03N\x86q\x14uX\x12\x00\x00\x00_allowed_overridesq\x15]q\x16(X\x05\x00\x00\x00indexq\x17X\x04\x00\x00\x00sortq\x18eh\x17cfunctools\npartial\nq\x19cbuiltins\ngetattr\nq\x1acsnakemake.io\nNamedlist\nq\x1bX\x0f\x00\x00\x00_used_attributeq\x1c\x86q\x1dRq\x1e\x85q\x1fRq (h\x1e)}q!X\x05\x00\x00\x00_nameq"h\x17sNtq#bh\x18h\x19h\x1e\x85q$Rq%(h\x1e)}q&h"h\x18sNtq\'bh\rh\x06h\x0fh\x07h\x11h\x08h\x13h\tubX\x06\x00\x00\x00outputq(csnakemake.io\nOutputFiles\nq))\x81q*X\x8c\x00\x00\x00results/Rad21/2021_01_19_Fbn2_Rad21F1M_noIAA_488nm_0p5_561nm_0p05_SC8_2X_30z_250nm_780V_365T_20s_movie1-01__processed_low..czi.tracks_ML.csvq+a}q,(h\x0b}q-X\x07\x00\x00\x00traj_mlq.K\x00N\x86q/sh\x15]q0(h\x17h\x18eh\x17h\x19h\x1e\x85q1Rq2(h\x1e)}q3h"h\x17sNtq4bh\x18h\x19h\x1e\x85q5Rq6(h\x1e)}q7h"h\x18sNtq8bh.h+ubX\x06\x00\x00\x00paramsq9csnakemake.io\nParams\nq:)\x81q;}q<(h\x0b}q=h\x15]q>(h\x17h\x18eh\x17h\x19h\x1e\x85q?Rq@(h\x1e)}qAh"h\x17sNtqBbh\x18h\x19h\x1e\x85qCRqD(h\x1e)}qEh"h\x18sNtqFbubX\t\x00\x00\x00wildcardsqGcsnakemake.io\nWildcards\nqH)\x81qI(X\x05\x00\x00\x00Rad21qJXp\x00\x00\x002021_01_19_Fbn2_Rad21F1M_noIAA_488nm_0p5_561nm_0p05_SC8_2X_30z_250nm_780V_365T_20s_movie1-01__processed_low..cziqKe}qL(h\x0b}qM(X\x0c\x00\x00\x00sample_groupqNK\x00N\x86qOX\x0b\x00\x00\x00sample_nameqPK\x01N\x86qQuh\x15]qR(h\x17h\x18eh\x17h\x19h\x1e\x85qSRqT(h\x1e)}qUh"h\x17sNtqVbh\x18h\x19h\x1e\x85qWRqX(h\x1e)}qYh"h\x18sNtqZbX\x0c\x00\x00\x00sample_groupq[hJX\x0b\x00\x00\x00sample_nameq\\hKubX\x07\x00\x00\x00threadsq]K\x01X\t\x00\x00\x00resourcesq^csnakemake.io\nResources\nq_)\x81q`(K\x01K\x01e}qa(h\x0b}qb(X\x06\x00\x00\x00_coresqcK\x00N\x86qdX\x06\x00\x00\x00_nodesqeK\x01N\x86qfuh\x15]qg(h\x17h\x18eh\x17h\x19h\x1e\x85qhRqi(h\x1e)}qjh"h\x17sNtqkbh\x18h\x19h\x1e\x85qlRqm(h\x1e)}qnh"h\x18sNtqobhcK\x01heK\x01ubX\x03\x00\x00\x00logqpcsnakemake.io\nLog\nqq)\x81qr}qs(h\x0b}qth\x15]qu(h\x17h\x18eh\x17h\x19h\x1e\x85qvRqw(h\x1e)}qxh"h\x17sNtqybh\x18h\x19h\x1e\x85qzRq{(h\x1e)}q|h"h\x18sNtq}bubX\x06\x00\x00\x00configq~}q\x7f(X\n\x00\x00\x00img_folderq\x80X?\x00\x00\x00/mnt/md0/Hansen Lab Dropbox/DataStorage/Imaging/Fbn2/processed/q\x81X\x11\x00\x00\x00samples_file_nameq\x82X\x12\x00\x00\x00config/samples.csvq\x83X\x0f\x00\x00\x00beads_file_nameq\x84X\x18\x00\x00\x00config/samples_beads.csvq\x85X\x14\x00\x00\x00ml_replication_modelq\x86X>\x00\x00\x00ml_models/Replication_model_2021_04_14_hugo_classifications.h5q\x87X\x11\x00\x00\x00ml_good_bad_modelq\x88XH\x00\x00\x00ml_models/4D_model_2021_04_08_mixedInputClean_v2_doubleClassification.h5q\x89X\x0c\x00\x00\x00samples_listq\x8a}q\x8b(X\x05\x00\x00\x00Rad21q\x8c}q\x8d(X\r\x00\x00\x00name_includesq\x8e]q\x8f(X\x04\x00\x00\x00.cziq\x90X\x05\x00\x00\x00rad21q\x91eX\r\x00\x00\x00name_excludesq\x92]q\x93(X\x05\x00\x00\x00beadsq\x94X\x04\x00\x00\x00fastq\x95X\x03\x00\x00\x00bfpq\x96X\n\x00\x00\x00vermicelliq\x97X\x06\x00\x00\x00advaitq\x98X\x07\x00\x00\x00neuronsq\x99X\n\x00\x00\x00nocodazoleq\x9aX\x03\x00\x00\x0060sq\x9bX\x04\x00\x00\x005minq\x9ceuX\x04\x00\x00\x00CTCFq\x9d}q\x9e(X\r\x00\x00\x00name_includesq\x9f]q\xa0(X\x04\x00\x00\x00.cziq\xa1X\x04\x00\x00\x00ctcfq\xa2eX\r\x00\x00\x00name_excludesq\xa3]q\xa4(X\x05\x00\x00\x00beadsq\xa5X\x04\x00\x00\x00fastq\xa6X\x03\x00\x00\x00bfpq\xa7X\n\x00\x00\x00vermicelliq\xa8X\x06\x00\x00\x00advaitq\xa9X\x07\x00\x00\x00neuronsq\xaaX\n\x00\x00\x00nocodazoleq\xabX\x03\x00\x00\x0060sq\xacX\x04\x00\x00\x005minq\xadeuX\x04\x00\x00\x00WAPLq\xae}q\xaf(X\r\x00\x00\x00name_includesq\xb0]q\xb1(X\x04\x00\x00\x00.cziq\xb2X\x04\x00\x00\x00waplq\xb3eX\r\x00\x00\x00name_excludesq\xb4]q\xb5(X\x05\x00\x00\x00beadsq\xb6X\x04\x00\x00\x00fastq\xb7X\x03\x00\x00\x00bfpq\xb8X\n\x00\x00\x00vermicelliq\xb9X\x06\x00\x00\x00advaitq\xbaX\x07\x00\x00\x00neuronsq\xbbX\n\x00\x00\x00nocodazoleq\xbcX\x03\x00\x00\x0060sq\xbdX\x04\x00\x00\x005minq\xbeeuX\x03\x00\x00\x00C27q\xbf}q\xc0(X\r\x00\x00\x00name_includesq\xc1]q\xc2(X\x04\x00\x00\x00.cziq\xc3X\x03\x00\x00\x00c27q\xc4eX\r\x00\x00\x00name_excludesq\xc5]q\xc6(X\x05\x00\x00\x00beadsq\xc7X\x04\x00\x00\x00fastq\xc8X\x03\x00\x00\x00bfpq\xc9X\n\x00\x00\x00vermicelliq\xcaX\x06\x00\x00\x00advaitq\xcbX\x07\x00\x00\x00neuronsq\xccX\n\x00\x00\x00nocodazoleq\xcdX\x03\x00\x00\x0060sq\xceX\x04\x00\x00\x005minq\xcfeuX\x03\x00\x00\x00C65q\xd0}q\xd1(X\r\x00\x00\x00name_includesq\xd2]q\xd3(X\x04\x00\x00\x00.cziq\xd4X\x03\x00\x00\x00c65q\xd5eX\r\x00\x00\x00name_excludesq\xd6]q\xd7(X\x05\x00\x00\x00beadsq\xd8X\x04\x00\x00\x00fastq\xd9X\x03\x00\x00\x00bfpq\xdaX\n\x00\x00\x00vermicelliq\xdbX\x06\x00\x00\x00advaitq\xdcX\x07\x00\x00\x00neuronsq\xddX\n\x00\x00\x00nocodazoleq\xdeX\x03\x00\x00\x0060sq\xdfX\x04\x00\x00\x005minq\xe0euX\x03\x00\x00\x00C36q\xe1}q\xe2(X\r\x00\x00\x00name_includesq\xe3]q\xe4(X\x04\x00\x00\x00.cziq\xe5X\x03\x00\x00\x00c36q\xe6eX\r\x00\x00\x00name_excludesq\xe7]q\xe8(X\x05\x00\x00\x00beadsq\xe9X\x04\x00\x00\x00fastq\xeaX\x03\x00\x00\x00bfpq\xebX\n\x00\x00\x00vermicelliq\xecX\x06\x00\x00\x00advaitq\xedX\x07\x00\x00\x00neuronsq\xeeX\n\x00\x00\x00nocodazoleq\xefX\x03\x00\x00\x0060sq\xf0X\x04\x00\x00\x005minq\xf1euuX\x1b\x00\x00\x00chromatic_shift_superblocksq\xf2]q\xf3(X\x1b\x00\x00\x00(\'2020-08-27\',\'2020-10-19\')q\xf4X\x1b\x00\x00\x00(\'2020-10-20\',\'2020-11-10\')q\xf5X\x1b\x00\x00\x00(\'2020-11-11\',\'2021-02-07\')q\xf6X\x1b\x00\x00\x00(\'2021-02-08\',\'2021-12-31\')q\xf7euX\x04\x00\x00\x00ruleq\xf8X \x00\x00\x00classify_replicated_trajectoriesq\xf9X\x0f\x00\x00\x00bench_iterationq\xfaNX\t\x00\x00\x00scriptdirq\xfbXA\x00\x00\x00/home/hbrandao/libs/data_analysis_Fbn2/snakemake/workflow/scriptsq\xfcub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/home/hbrandao/libs/data_analysis_Fbn2/snakemake/workflow/scripts/classify_dots_by_CNN.py';
######## snakemake preamble end #########
## Import required modules
import numpy as np 
import pandas as pd
import pickle
from tqdm.keras import TqdmCallback
import keras
from connect_the_dots.filtering import create_mask_symmetric

verbose = True

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
        if verbose == True:
            print(f"Doing frame {frame}")
        
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
