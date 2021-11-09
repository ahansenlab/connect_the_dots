from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# container: "docker://continuumio/miniconda3"

##### Load config and sample sheets #####
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples_df = pd.read_csv(config["samples_file_name"], dtype = str)
validate(samples, schema="../schemas/samples.schema.yaml")

"""

##### wildcard constraints #####
wildcard_constraints:
    grou_name = "|".join(samples_df["group"]),
    file_names = "|".join(samples_df["file_names"])
    file_paths = "|".join(samples_df["file_names"])



####### helpers ###########
def get_samples(wildcards):
    """Get raw .czi files from the samples sheet."""
    os.path.join()
    #u = samples_df.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
    #return [ f"{u.fq1}", f"{u.fq2}" ]
    
"""