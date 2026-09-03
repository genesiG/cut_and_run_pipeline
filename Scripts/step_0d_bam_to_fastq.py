#!/usr/bin/env python3

"""
step_0d_bam_to_fastq.py - Optional step to convert BAM files back to FASTQ format.
                        - Can be used when the original FASTQ files are not available.


EXAMPLE USAGE:
    python3 step_0d_bam_to_fastq.py
"""

import os
import time
import config
from utils import ldsample

def convert_bam_to_fastq(sample_name):
    """
    Submits a batch job to align reads and do post-processing.
    """

    # Create the main alignment directory and log directory
    work_dir = os.path.join(config.ORIGINAL_DATA, "fastq_files")
    original_bam_dir = config.ORIGINALBAMDIR
    log_dir = os.path.join(work_dir, "log")

    # Define job name based on sample name
    data_id = sample_name
    job_name = f"{data_id}.bamtofastq"

    # Input files
    input_file = os.path.join(original_bam_dir, f"{data_id}.bam")

    # Output files
    output_file = os.path.join(work_dir, f"{data_id}.fastq")

    if config.IS_PAIRED_END:
        out_1 = os.path.join(work_dir, f"{data_id}.R1.fastq")
        out_2 = os.path.join(work_dir, f"{data_id}.R2.fastq")

        command_line = f"samtools fastq -1 {out_1} -2 {out_2}"
        gzip_cmd = f"gzip {out_1} && gzip {out_2}"
    else:
        command_line = f"samtools fastq {input_file} > {output_file}"
        gzip_cmd = f"gzip {output_file}"

    # Log paths
    log_out_file = os.path.join(log_dir, f"{data_id}.out")
    log_err_file = os.path.join(log_dir, f"{data_id}.error")
    batch_file = os.path.join(log_dir, f"{data_id}.batch")

    # Create directories if they do not exist
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Construct the batch command
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -n 8,32  # number of cores to request
#BSUB -M 64400
#BSUB -R "rusage [mem=64400] span[hosts=1]"    # Reserve the different CPU cores under the same node

conda activate chipseq

{command_line}
{gzip_cmd}
"""

    # Write the batch command to a file
    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    # Submit the batch job
    print(f"bsub < {batch_file}")
    os.system(f"bsub < {batch_file}")

    # Sleep for 1 second to avoid flooding the scheduler
    time.sleep(1)

def main():
    """
    Main function that:
    1. Loads SAMPLES_CTL mapping from config
    2. Iterates over each sample
    3. Calls convert_bam_to_fastq()
    """
    # Read the samples.txt file
    ldsample.load_samples()

    # Iterate over the samples
    for sample_name in ldsample.SAMPLES.keys():
        convert_bam_to_fastq(sample_name)


if __name__ == "__main__":
    main()
