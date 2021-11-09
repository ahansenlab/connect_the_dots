import pandas as pd
import tracklib as tl
import tracklib
import numpy as np
import pandas as pd

# identify sample group
group = snakemake.wildcards.sample_group 
data_groups = [group] # dummy variable

fx = snakemake.config['pixel_size_X_microns']
fy = snakemake.config['pixel_size_Y_microns']
fz = snakemake.config['pixel_size_Z_microns']
ft = snakemake.config['time_step_seconds']


# identify the subset of post-ML trajectory files pertaining to this group
files_list = [f for f in snakemake.input.traj_ml if f"results/{group}/" in f]

tag_set_file = snakemake.output.tagged_set_file

do_good_only = True
do_truncation = False 
do_good_avg = True
do_tier3 = snakemake.params.do_tier_3 


for gi in range(len(data_groups)):
    new_df_list = []
    running_count = 0
   
    if do_tier3 == False:
        summary_file_list = [f for f in snakemake.input.summary_file if f"{group}" in f]
        summary_file = summary_file_list[gi]        
        summary_df = pd.read_csv(summary_file)
    else:
        summary_df = pd.read_excel(snakemake.input.summary_file,sheet_name=f"{group}")
        
    files = files_list
    
    movie_idx = 0
    for traj_file in sorted(files):
        movie_idx += 1
        tracks_df = pd.read_csv(traj_file)
        has_name = summary_df.Filename.str.contains(traj_file.split('/')[-1].split('.czi')[0])

        if sum(has_name)==0:
            print(f"Missing file:  \n{traj_file.split('/')[-1].split('.czi')[0]}")
            continue            
            
        particles = tracks_df.particle.unique()
        for p in particles:

            has_particle = (summary_df.Particle==p)

            # check that some file exists
            if sum(has_name & has_particle)== 0:
                print(f"Missing file particle #{p} from file: \n{traj_file.split('/')[-1].split('.czi')[0]}")
                continue
            
            if do_tier3 == False:
                # skip files if ML_replicated is "True"
                if summary_df[has_name & has_particle].ML_replicated.values[0]==True:
                    continue                       
            else:
                # skip files if the Keep column is not "T"
                if not summary_df[has_name & has_particle].Keep.values[0]=='T':
                    continue                   
                
            T0 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==0)].frame.values
            X0 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==0)].x.values*fx
            Y0 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==0)].y.values*fy
            Z0 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==0)].z.values*fz
            G0 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==0)].ML_class.values
            T1 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==1)].frame.values
            X1 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==1)].x.values*fx
            Y1 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==1)].y.values*fy
            Z1 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==1)].z.values*fz
            G1 = tracks_df[(tracks_df.particle==p) & (tracks_df.channel==1)].ML_class.values

            if len(T0) > 0 and len(T1)> 0:
                
                if do_good_only==True:

                    # now, do only good frame
                    commG0 = np.array([int(x) for x in np.where(G0=='Good')[0]])
                    commG1 = np.array([int(x) for x in np.where(G1=='Good')[0]])
                    if len(commG0)==0 or len(commG1)==0:
                        continue

                    T0 = T0[commG0]
                    X0 = X0[commG0]
                    Y0 = Y0[commG0]
                    Z0 = Z0[commG0]
                    T1 = T1[commG1]
                    X1 = X1[commG1]
                    Y1 = Y1[commG1]
                    Z1 = Z1[commG1]

                overlap,comm0,comm1 = np.intersect1d(T0,T1,return_indices=True)
                
                D = np.sqrt( (X0[comm0]-X1[comm1])**2 + (Y0[comm0]-Y1[comm1])**2 + (Z0[comm0]-Z1[comm1])**2 )

                if len(T0)==0 or len(T1)==0:
                    continue

                running_count += 1
                #date = traj_file.split('/')[-1].split('_Fbn2')[0]
                date = "_".join(traj_file.split('/')[-1].split('_')[0:3])
                movie_number = date+"_"+traj_file.split('20s_')[1].split('_proc')[0]
                new_df = pd.DataFrame({'id':[running_count]*len(comm0),'t':T0[comm0],
                                       'x':X0[comm0],'y':Y0[comm0],'z':Z0[comm0], 
                                       'x2':X1[comm1],'y2':Y1[comm1],'z2':Z1[comm1],
                                       'date':[date]*len(comm0),'movie_index':[movie_number]*len(comm0)})
                new_df_list.append(new_df)     

    new_df = pd.concat(new_df_list)
    new_df.to_csv(tag_set_file,sep='\t',index=False)


