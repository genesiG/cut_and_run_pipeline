#!/usr/bin/env python3

"""
step_4b_bigwigtobedgraph.py - (Optional) script to convert bigWig files into bedGraph

"""

### Import modules
import os
import time
import pandas as pd

import config
###

# Working directory and samples file with scaling factors
work_dir = os.path.join(config.BEDGRAPHDIR)
log_dir = os.path.join(work_dir, "log")
bw_dir = os.path.join(config.IMPORTABLE_DATA, "bigwig")
bw_sufix = "normalized.bw"
scaling_file = "scale_UNC_vs_DMSO_SF.txt"


# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

        
def run_bigwigtobedgraph(sample_name):
    """
    Creates and submits a batch script to generate scaled bigWig files using bamCoverage.
    """
    
    # Remove suffix from filename
    basename = sample_name.replace('.bw','.bedgraph')
    
    # Job name based on sample name
    job_name = f"{basename}.bigwigtobedgraph"

    # Input files
    input_bw = os.path.join(bw_dir, f"{basename}.{bw_sufix}")
    
    # Output files
    output_bedgraph = os.path.join(work_dir, f"{basename}.bedgraph")

    # Log files
    log_out_file = os.path.join(log_dir, f"{job_name}.out")
    log_err_file = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    # Build the batch command
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 8,32  # number of cores to request
#BSUB -M 64400
#BSUB -R "rusage [mem=64400] span[hosts=1]"    # Reserve the different CPU cores under the same node

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq


bigWigToBedGraph {input_bw} {output_bedgraph} 

"""

    # Write the batch command to a file
    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"Submitting job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    """
    Loads scaling factors from SF.txt and processes each sample.
    """
    # Read the scaling factors file
    sf_file = os.path.join(config.SCALINGDIR, scaling_file)
    sf_df = pd.read_csv(sf_file, sep='\t')

    # Iterate over each sample in the DataFrame
    for _, row in sf_df.iterrows():
        sample_name = row['ID']
        # Remove suffix from filename
        sample_name = sample_name.replace('.final.sorted.bam','')
        
        # Submit job
        run_bigwigtobedgraph(sample_name)

if __name__ == "__main__":
    main()