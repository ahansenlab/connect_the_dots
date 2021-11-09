How to run the pipeline (2021/04/26): 


## Installation
1. Run:
```pip install snakemake
```

## Optional steps to set up the pipeline run

1. Find some appropriate imaging parameters to run `connect_the_dots`. 
    - Run the notebook called **scripts/00_get_pipeline_parameters.ipynb** to find the imaging parameters for the pipeline
    - Most of the parameter changes that you will need to make are inside of the function called *get_localizations_iterative()*. 
    - Once you are happy with the parameters, you can update the **snakefile_fixed_mask** file with your modifications

2.  Next, set up the configuration files (which can be found in the config/ folder)
    - Set up the **config_beads.yaml** file (for beads only)
    - Set up the **config.yaml file** (for samples)
    
    
## Run the pipeline for samples files initial processing. 

0. Change to the correct folder: `cd /home/hbrandao/libs/data_analysis_Fbn2/snakemake/workflow`

1. To make or update the samples.csv file, run:
`snakemake -s prep_samples_snakefile --cores 1 --forceall`
   - Check that the **config/samples.csv** file has updated and contains all the files you want

2. Run: `snakemake -s snakefile_fixed_mask --resources mem_gb=240 --cores 1 --keep-going` 


3. Make the MIP movie with trajectories overlayed. Run the notebook **scripts/create_movies_with_trajectories.ipynb**.


## Run the pipeline for beads files and all downstream analyses (including ML and chromatic shift correction)

4. To make the samples_beads.csv file, run:
`snakemake -s prep_beads_snakefile --cores 1 --forceall`
   - Check that the **config/samples_beads.csv** file has updated and contains all the files you want

5. Run `snakemake -s snakefile_beads --resources mem_gb=240 --cores 40 --keep-going`


6. Optionally generate the snakemake report file. Run `snakemake -s snakefile_beads --report`