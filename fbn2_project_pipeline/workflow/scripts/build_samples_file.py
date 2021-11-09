import pandas as pd
from pathlib import Path
import os

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


img_folder = snakemake.config["img_folder"]
samples_list = snakemake.config["samples_list"]
sample_groups = snakemake.config["sample_groups"]

filenames_list = []
filepath_list = []
group_list = []

# map sample keys to sample groups
sample_name_to_group = {}
for k in sample_groups.keys():    
    for kk in sample_groups[k]:
        sample_name_to_group[kk] = k   

for k in samples_list.keys():
    # get all the files matching the given keys
    name_includes = samples_list[k]["name_includes"]
    name_excludes = samples_list[k]["name_excludes"]
    names, paths = search_for(img_folder,name_includes,name_excludes)
    if len(names)>0:
        filenames_list.extend(names)
        filepath_list.extend(paths)
        group = sample_name_to_group[k]
        group_list.extend([group]*len(names))

df = pd.DataFrame({'group': group_list,
                    'file_names':filenames_list,
                   'file_paths':filepath_list,})
df.to_csv(snakemake.output[0],index=False)



