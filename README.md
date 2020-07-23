# Connect the dots
3D "dot" tracking for pairs of tagged chromatin loci. 

## Installation instructions

There are a few dependencies that need to be installed before one can run `connect_the_dots`.

1. Install `BioFormats` 

    - https://github.com/ome/bioformats   
    - After downloading, install via the command line: run `ant jars`
    - https://docs.openmicroscopy.org/bio-formats/5.8.2/developers/building-bioformats.html#source-obtain-and-build
    
    
2. Install Java-Development-Kit

    - https://openjdk.java.net/install/
    

3. Install `javabridge`
    
    - IMPORTANT: Make sure you have correctly set the JAVA_HOME variable so that javabridge knows where to look for Java (see issue here: https://github.com/LeeKamentsky/python-javabridge/issues/152)
    - `Boto3` is a required dependency, so run: `pip install boto3`
    - Finally, run `pip install javabridge`


4. Install `python-bioformats`

    - run `pip install python-bioformats`
    
    
5. Install `tifffile`

    - This is used to parse metadata
    - run `pip install tifffile`   


6. Now, you should be all set!

    - run `python setup.py develop`, and you are ready to go!

## References (to cite)

    - Bioformats: https://www.ncbi.nlm.nih.gov/pubmed/20513764


