"""
 ldsample.py - Module to load sample metadata
 
"""

import config
import os
import pandas as pd
import glob
import warnings

# Initialize dictionaries
SAMPLES = {}       # Keys: sample names, Values: sequencing file paths
SAMPLES_CTL = {}   # Keys: sample names, Values: corresponding input control sample names

def load_samples(filename=None):
    if filename is None:
        filename = "samples.txt" # Default
        
    # Read in sample files from filename using pandas
    sample_metadata = os.path.join(config.METADATA, filename)
    
    if not os.path.exists(sample_metadata):
        raise FileNotFoundError(
                f"Cannot open sample metadata file at '{sample_metadata}'." 
                )
            
    # Read the sample list into a DataFrame
    samples_df = pd.read_csv(sample_metadata, sep='\t')
    
    # Clear dictionaries to avoid retaining data from previous runs
    SAMPLES.clear()
    SAMPLES_CTL.clear()
    
    
    # Iterate over the DataFrame rows
    for index, row in samples_df.iterrows():
        sample_name = row['sample.name']
        control_sample = row['sample.control']
        
        # Initialize filepath as None
        filepath = None
        
        # Search for files with possible extensions
        for ext in config.EXTENSIONS:
            # Construct the search pattern
            search_pattern = os.path.join(config.DATADIR, f"{sample_name}{ext}")
            # Use glob to find matching files
            matching_files = glob.glob(search_pattern)
            if matching_files:
                filepath = matching_files[0]  # Take the first matching file
                break  # Exit the loop once a file is found
        
        if filepath:
            SAMPLES[sample_name] = filepath
        else:
            # No file found with either extension
            warnings.warn(
                f"No sequencing file found for sample '{sample_name}' in directory '{config.DATADIR}'"
            )
        
        if control_sample != "-":
            SAMPLES_CTL[sample_name] = control_sample
