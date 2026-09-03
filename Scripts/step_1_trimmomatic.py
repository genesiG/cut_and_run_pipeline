#!/usr/bin/env python3

""" 
trimmomatic.py  - (Optional) Script to perform adapter trimming using Trimmomatic
                
NOTE:
    - In CUT&RUN, read lengths are shorter than in ChIP-seq.
    - Therefore, trimming and alignment settings must be adjusted to account for this. 
    - For this reason, trimming is NOT AVISED advised in CUT&RUN and during alignment. 
    - It is also advised that the dovetailed read mates are considered for alignment.
    - https://github.com/hbc/knowledgebase/blob/master/chipseq/cutandrun.md

INSTRUCTIONS:
    - First, run 'step_0c_prepare_samples_file.py' to generate the "paired_samples.txt" file
    - Run this script in the terminal as "python3 step_1_trimmomatic.py"
    - (Optional) check read quality using 'step_0e_fastqc.py' 

"""

### Import modules
import os
import time
import config
from utils import ldsample

### Define function to search for sequencing files
def find_input_file(data_dir, sample_name):
    extensions = config.EXTENSIONS
    for ext in config.EXTENSIONS:
        file_path = os.path.join(data_dir, f"{sample_name}{ext}")
        if os.path.exists(file_path):
            return file_path
    raise FileNotFoundError(f"Input file for sample {sample_name} not found in {data_dir} with extensions {extensions}")
###

### Define function to perform trimming of sequencing reads
def trim(sample_name, paired_sample = None):
    data_dir = config.DATADIR
    work_dir = config.TRIMDIR
    log_dir = os.path.join(config.TRIMDIR, "log")
    job_name = f"{sample_name}.trim"

    if config.EXPERIMENT in ("atacseq", "cutandtag"):
        adapter_file = "NexteraPE-PE.fa"
    else:
        adapter_file = "TruSeq3-PE-2.fa" if paired_sample else "TruSeq3-SE.fa"

    if paired_sample:
        # Input files
        in_1 = find_input_file(data_dir, sample_name)
        in_2 = find_input_file(data_dir, paired_sample)

        # Output files
        o_1 = os.path.join(work_dir, os.path.basename(in_1))
        o_2 = os.path.join(work_dir, os.path.basename(in_2))
        o_1_unpaired = os.path.join(work_dir,
                                    f"{os.path.splitext(os.path.basename(in_1))[0]}_unpaired{os.path.splitext(in_1)[1]}")
        o_2_unpaired = os.path.join(work_dir,
                                    f"{os.path.splitext(os.path.basename(in_2))[0]}_unpaired{os.path.splitext(in_2)[1]}")

        # Command line parameters
        seq_parameter = "PE"
        adapter_parameter = f"{adapter_file}:2:30:10:8:True"
        sample_parameters = f"{in_1} {in_2} {o_1} {o_1_unpaired} {o_2} {o_2_unpaired}"

        # Verify read length
        read_len_par = f"""
# Verify minimum read length
echo "\\nRead length verification:"

# For R1
zcat {o_1} | awk 'NR%4==2 {{print length}}' | sort | uniq -c | awk '$2 < 25 {{print "WARNING: Reads shorter than 25bp found (forward pair):", $0}}'
# For R2
zcat {o_2} | awk 'NR%4==2 {{print length}}' | sort | uniq -c | awk '$2 < 25 {{print "WARNING: Reads shorter than 25bp found (reverse pair):", $0}}'
"""

    else:
        # Construct input and output file paths for single-end sequencing
        # Input and output files for single-end sequencing
        in_file = find_input_file(data_dir, sample_name)
        o_file = os.path.join(work_dir, os.path.basename(in_file))

        # Command line parameters
        seq_parameter = "SE"
        adapter_parameter = "{adapter_file}:1:30:10:1"
        sample_parameters = f"{in_file} {o_file}"

        read_len_par = ""


    trim_params = ("LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25"
                   if config.EXPERIMENT == "atacseq"
                   else "MINLEN:25")

    # Log files
    log_out_file = os.path.join(log_dir, f"{sample_name}.log")
    log_err_file = os.path.join(log_dir, f"{sample_name}.error")
    batch_file = os.path.join(log_dir, f"{sample_name}.batch")

    # Check if outputs already exist to avoid redundant jobs
    outputs_exist = False
    if paired_sample and os.path.exists(o_1) and os.path.exists(o_2):
        outputs_exist = True
    elif not paired_sample and os.path.exists(o_file):
        outputs_exist = True

    if outputs_exist:
        print(f"Skipping {job_name} - output files already exist.")
        return

    # Create directories if they don't exist
    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    # Build the batch command script
    batch_cmd = f"""
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 8  # number of cores to request
#BSUB -M 32000
#BSUB -R "rusage [mem=4000] span[hosts=1]"    # Reserve the different CPU cores under the same node

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running Trimmomatic with parameters:"
echo "Adapter file: {adapter_file}"
echo "Trim params: {trim_params}"

trimmomatic \\
{seq_parameter} \\
-threads 16 \\
{sample_parameters} \\
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic-0.40-0/adapters/{adapter_parameter} \\
{trim_params}


{read_len_par}

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    # Write the batch command to a file
    with open(batch_file, 'w', encoding="utf-8") as batch:
        batch.write(batch_cmd)

    # Submit the batch job
    print(f"Submitting job {job_name}")
    os.system(f"bsub < {batch_file}")

    # Sleep for 1 second between submissions
    time.sleep(1)
###

### Define main function loop
def main():
    if config.IS_PAIRED_END:
        # Read the paired_samples.txt file
        ldsample.load_samples("paired_samples.txt")

        # Iterate over the samples
        for sample_name, paired_sample in ldsample.SAMPLES_CTL.items():
            if paired_sample:
                trim(sample_name, paired_sample)
            else:
                print(f"WARNING: Paired sample for {sample_name} is missing. Check the 'paired_samplex.txt' file at {config.METADATA}")

    else:
        # Read the samples.txt file
        ldsample.load_samples()

        # Iterate over the samples
        for sample_name in ldsample.SAMPLES.keys():
            trim(sample_name)
###

if __name__ == "__main__":
    main()
