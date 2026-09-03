#!/usr/bin/env python3

"""
step_4a_chipseqspikeinfree.py - Script to submit normalization R script job
Steps:
1. Set up log directory using config
2. Build batch script for normalization job using NORMALIZATION_RUNS
3. Submit to cluster using bsub
"""

### Import modules
import os
import time
import config

### Set up parameters
SAMPLES_FILE = "sample_metadata.txt"
SUFFIX = (f".qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
          if getattr(config, 'REMOVE_DUPLICATES', False) else f".qc.sort.markdup.mapq{config.MAPQ}.final.bam")

### Set up directories
work_dir = os.path.join(config.SCALINGDIR)
output_dir = os.path.join(work_dir, "chipseqspikeinfree")
SCRIPT_PATH = os.path.join(config.CODEDIR, "utils", "ChIPseqSpikeInFree.R")
INSTALL_SCRIPT = os.path.join(config.CODEDIR, "utils", "SpikeInFreeInstall.R")

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

INSTALL_LOG_DIR = os.path.join(output_dir, "log")
os.makedirs(INSTALL_LOG_DIR, exist_ok=True)

def submit_install_job():
    """
    Submit a single, BLOCKING bsub job (-K flag) that installs
    ChIPseqSpikeInFree once into the shared conda environment.

    Using -K causes os.system() to wait until the job completes before
    returning, so all subsequent normalization jobs are guaranteed to find
    the package already installed.
    """
    job_name   = "spikein_free_install"
    log_out    = os.path.join(INSTALL_LOG_DIR, "install.log")
    log_err    = os.path.join(INSTALL_LOG_DIR, "install.error")
    batch_file = os.path.join(INSTALL_LOG_DIR, "install.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J "{job_name}"
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 1
#BSUB -M 4000
#BSUB -R "rusage [mem=4000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

Rscript {INSTALL_SCRIPT}
"""

    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    print("Submitting package installation job (blocking until complete) ...")
    # -K: bsub waits for job completion before returning
    ret = os.system(f"bsub -K < {batch_file}")
    if ret != 0:
        raise RuntimeError(
            "ChIPseqSpikeInFree installation job failed. "
            f"Check logs in {INSTALL_LOG_DIR}"
        )
    print("Installation complete. Submitting normalization jobs ...\n")

def submit_normalization_job(run_name, params):
    """Creates and submits batch script for normalization analysis"""
    target = params["targets"]
    bin_size = params["bin_size"]
    cutoff = params["cutoff"]
    max_turns = params["max_turns"]

    log_dir = (os.path.join(output_dir, "log")
               if bin_size == 10000 else os.path.join(output_dir, f"log_{bin_size}bp_bins"))
    os.makedirs(log_dir, exist_ok=True)

    job_name = f"chipseqspikeinfree_{run_name}"

    log_out = os.path.join(log_dir, f"{run_name}.log")
    log_err = os.path.join(log_dir, f"{run_name}.error")
    batch_file = os.path.join(log_dir, f"{run_name}.batch")

    # Build batch command
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J "{job_name}"
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 4
#BSUB -M 64000
#BSUB -R "rusage [mem=16000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Run R script
Rscript {SCRIPT_PATH} \\
    {bin_size} \\
    "{SAMPLES_FILE}" \\
    "{SUFFIX}" \\
    {cutoff} \\
    {max_turns} \\
    "{target}" \\
    "{output_dir}" \\
    {str(config.IS_PAIRED_END).upper()} \\
    "{run_name}" \\

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    # Write batch file
    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    print(f"Submitting normalization job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    """Main function to submit jobs."""
    # Step 1: install package once, blocking, before any parallel runs
    submit_install_job()
    # Step 2: fan out all normalization runs in parallel
    for run_name, params in config.NORMALIZATION_RUNS.items():
        submit_normalization_job(run_name, params)

if __name__ == "__main__":
    main()
