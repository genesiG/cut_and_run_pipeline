#!/usr/bin/env python3

"""
fastqc.py - Script to run FASTQC program for quality control analysis of the fastq files
          - run the script before and/or after trimmomatic.py
"""

import os
import time
import config
from utils import ldsample

# Load samples
ldsample.load_samples()  # This loads the samples into ldsample.SAMPLES

def fastqc(sample_name):
    """Run fastqc on a single sample."""
    data_id = sample_name
    work_dir = config.QCDIR1
    job_name = f"{data_id}.fastqc1"

    if config.USE_TRIMMOMATIC:
        in_file = os.path.join(config.TRIMDIR, f"{data_id}.fastq.gz")
        output_dir = os.path.join(work_dir, "trimmed")
    else:
        in_file = os.path.join(config.DATADIR, f"{data_id}.fastq.gz")
        output_dir = os.path.join(work_dir, "untrimmed")

    # Output file paths
    log_dir = os.path.join(work_dir, "log")
    log_out_file = os.path.join(log_dir, f"{data_id}.out")
    log_err_file = os.path.join(log_dir, f"{data_id}.error")
    batch_file = os.path.join(log_dir, f"{data_id}.batch")

    # Create directories if they don't exist
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    # Build the batch command script
    batch_cmd = f"""
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -M 40960    # Memory in MB
#BSUB -n 4       # Number of cores
#BSUB -R "rusage [mem=10240] span[hosts=1]"    # Reserve the different CPU cores under the same node 

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

fastqc -o {output_dir} --extract {in_file}

rm {output_dir}/*fastqc.zip
"""

    # Write the batch command to a file
    with open(batch_file, 'w', encoding="utf-8") as batch:
        batch.write(batch_cmd)

    # Submit the batch job
    print(f"bsub < {batch_file}")
    os.system(f"bsub < {batch_file}")

    # Sleep for 1 second between submissions
    time.sleep(1)

def main():
    """Run fastqc on all samples."""
    # Iterate over each sample and run fastqc
    for sample_name in ldsample.SAMPLES.keys():
        fastqc(sample_name)

if __name__ == "__main__":
    main()
