#!/usr/bin/env python3

"""
align_bowtie2.py - Uses Bowtie2 to align reads to reference genome

NOTE:
    - In CUT&RUN, read lengths are shorter than in ChIP-seq.
    - Therefore, trimming and alignment settings must be adjusted to account for this. 
    - For this reason, trimming is NOT AVISED in CUT&RUN and during alignment. 
    - It is also advised that the dovetailed read mates are considered for alignment.
    - https://github.com/hbc/knowledgebase/blob/master/chipseq/cutandrun.md

Converts SAM to BAM, filters and sorts the resulting file.

Steps:
1. Decide where to read fastq files from (TRIMDIR if trimming was done, else DATADIR).
2. For paired-end data, retrieve forward and reverse fastq from SAMPLES_CTL.
3. Build a batch script for Bowtie2 + samtools commands.
4. Submit the script with bsub.
5. Sleep briefly to avoid cluster flooding.

EXAMPLE USAGE:
    python3 step_2_align_bowtie2.py
"""

import os
import time
import config
from utils import ldsample

NUM_CORES = 16
MAX_MEM = 96000
MEM_PER_CORE = MAX_MEM/NUM_CORES

def run_alignment(sample_name, paired_sample=None):
    """
    Submits a batch job to align reads and do post-processing.
    """

    # Create the main alignment directory and log directory
    work_dir = config.ALIGNDIR
    bam_dir = config.BAMDIR
    original_bam_dir = config.ORIGINALBAMDIR
    processed_bam_dir = config.PROCESSEDBAMDIR
    log_dir = os.path.join(work_dir, "log")

    # Define job name based on sample name
    data_id = sample_name
    job_name = f"{data_id}.align"

    # Set fastq directory
    if config.USE_TRIMMOMATIC:
        datadir = config.TRIMDIR
    else:
        datadir = config.DATADIR

    # Input files
    if config.IS_PAIRED_END:
        in_1 = os.path.join(datadir, f"{data_id}.fastq.gz")
        in_2 = os.path.join(datadir, f"{paired_sample}.fastq.gz")

        sample_parameters = f"""-1 {in_1} \\
          -2 {in_2} \\
          --no-discordant \\
          --no-mixed"""

    else:
        in_file = os.path.join(datadir, f"{data_id}.fastq.gz")
        sample_parameters = f"-U {in_file}"

    # Output files
    output_prefix = os.path.join(work_dir, sample_name)
    output_sam = f"{output_prefix}.sam"
    output_bam = f"{output_prefix}.bam"


    # Output config
    output_prefix = os.path.join(work_dir, data_id)
    alignment_prefix = os.path.join(original_bam_dir, data_id)

    # Save BAM files to importable data dir
    output_sam = f"{alignment_prefix}.sam"
    output_bam = f"{alignment_prefix}.bam"

    # Log paths
    log_out_file = os.path.join(log_dir, f"{data_id}.log")
    log_err_file = os.path.join(log_dir, f"{data_id}.error")
    batch_file = os.path.join(log_dir, f"{data_id}.batch")

    # CUT and RUN specific parameter
    if config.EXPERIMENT == "cutandrun":
        align_param = "--dovetail"
    elif config.EXPERIMENT == "atacseq":
        # accept fragments up to 2000bp as concordant
        align_param = "--dovetail -X 2000"
    else:
        align_param = ""

    # Alingment command
    if config.IS_RAW_DATA_FASTQ:
        command_line = f"""\
          ## Align using Bowtie2 (output SAM file)
          bowtie2 -p {NUM_CORES} \\
          --local --very-sensitive-local \\
          -q \\
          -t \\
          -x {config.GNM_IDX} \\
          {sample_parameters} \\
          {align_param} \\
          -S {output_sam}
                    
          ## Convert SAM to BAM, filter, sort, and index
          samtools view \\
          -h \\
          -b \\
          -S \\
          {output_sam} > {output_bam}
          
          samtools sort -m 2G -@ {NUM_CORES} {output_bam} -o {output_bam}.tmp
          mv {output_bam}.tmp {output_bam}
          samtools index {output_bam}"""
    else:
        command_line = ""

    # Create directories if they do not exist
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(bam_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(original_bam_dir, exist_ok=True)
    os.makedirs(processed_bam_dir, exist_ok=True)

    # Check if outputs already exist to avoid redundant jobs
    if os.path.exists(output_bam) and os.path.exists(f"{output_bam}.bai") and os.path.getsize(output_bam) > 0:
        print(f"Skipping {job_name} - output files already exist.")
        return

    # Construct the batch command
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -n {NUM_CORES}  # number of cores to request
#BSUB -M {MAX_MEM}
#BSUB -R "rusage [mem={MEM_PER_CORE}] span[hosts=1]"    # Reserve the different CPU cores under the same node


source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

{command_line}

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

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
    3. Calls run_alignment()
    """
    # If you keep your sample->paired-sample mapping in config.SAMPLES_CTL:
    if config.IS_PAIRED_END:
        # Read the paired_samples.txt file
        ldsample.load_samples("paired_samples.txt")

        # Iterate over the samples
        for sample_name, paired_sample in ldsample.SAMPLES_CTL.items():
            if paired_sample:
                run_alignment(sample_name, paired_sample)
            else:
                print(f"WARNING: Paired sample for {sample_name} is missing. \
                      Check the 'paired_samples.txt' file at {config.METADATA}")

    else:
        # Read the samples.txt file
        ldsample.load_samples()

        # Iterate over the samples
        for sample_name in ldsample.SAMPLES.keys():
            run_alignment(sample_name)


if __name__ == "__main__":
    main()
