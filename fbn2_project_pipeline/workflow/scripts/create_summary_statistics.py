import pandas as pd
        
# identify sample group
group = snakemake.wildcards.sample_group 

# identify the subset of files pertaining to this group
files = [f for f in snakemake.input.traj_ml if f"results/{group}/" in f]

# make summary spreadsheet
idx_list = []
f_list = []
p_list = []
minframe_list = []
maxframe_list = []
g_count = []
b_count = []
n_count = []
u_count = []
bad_list = []

gr_count = []
gu_count = []
br_count = []
bu_count = []

for fi, file in enumerate(sorted(files)):
    
    traj_dict = pd.read_csv(file)
    short_name = (file.split('/')[-1]).split('.track')[0]
    
    particles = traj_dict.particle.unique()
    for p in sorted(particles):

        min_frame = traj_dict[(traj_dict.particle==p)].frame.min()
        max_frame = traj_dict[(traj_dict.particle==p)].frame.max()

        bad_frames = traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Bad') ].frame.unique()
        good_frames =  traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Good') ].frame.unique()
        nodot_frames =  traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='NoDot') ].frame.unique()
        nocat_frames =  traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Unclassified') ].frame.unique()

        good_rep_frames = traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Good') & 
                                    (traj_dict['ML_isReplicated']=='Replicated') ].frame.unique()
        good_unrep_frames = traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Good') & 
                                    (traj_dict['ML_isReplicated']=='Unreplicated') ].frame.unique()        
        bad_rep_frames = traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Bad') & 
                                    (traj_dict['ML_isReplicated']=='Replicated') ].frame.unique()        
        bad_unrep_frames = traj_dict[(traj_dict.particle==p) & (traj_dict['ML_class']=='Bad') & 
                                    (traj_dict['ML_isReplicated']=='Unreplicated') ].frame.unique()        


        idx_list.append(fi)
        f_list.append(short_name)
        p_list.append(p)        
        minframe_list.append(min_frame)
        maxframe_list.append(max_frame)
        g_count.append(len(good_frames) )
        b_count.append(len(bad_frames))
        n_count.append(len(nodot_frames))
        u_count.append(len(nocat_frames))
        bad_list.append(sorted(bad_frames))
        gr_count.append(len(good_rep_frames))
        gu_count.append(len(good_unrep_frames))
        br_count.append(len(bad_rep_frames))
        bu_count.append(len(bad_unrep_frames))


summary_df = pd.DataFrame({'idx':idx_list,'Filename':f_list,
             'Particle':p_list,
              'Start':minframe_list,
              'End':maxframe_list,              
             'NumGood':g_count,
             'NumBad':b_count,
             'NumNone':n_count,
             'NumUncategorized':u_count,
             'NumGoodUnreplicated':gu_count,
             'NumBadReplicated':br_count,
             'NumGoodReplicated':gr_count,
             'NumBadUnreplicated':bu_count,
             'Bad_list':bad_list})

summary_df['ML_replicated']  = (summary_df.NumBadReplicated/
                               (summary_df.NumBadReplicated+summary_df.NumGoodUnreplicated)>0.33)
summary_df.to_csv(snakemake.output.summary_sheets,index=False)

