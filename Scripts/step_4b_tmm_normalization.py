#!/usr/bin/env python3

"""
step_4b_tmm_normalization.py - Script to submit TMM normalization job
Steps:
1. Set up parameters for TMM normalization
2. Build batch script to run R normalization pipeline
3. Submit to cluster using bsub
"""

import os
import time
import config

### Set up parameters

# Build a comma-separated targets string from NORMALIZATION_RUNS for scaleFactors.R.
# Each entry's 'targets' regex selects BAM files for that run; one entry per antibody
# means all groups (DMSO, EEDi, EPZ) are normalised together, capturing global changes.
TARGETS = ",".join(
    params["targets"]
    for params in config.NORMALIZATION_RUNS.values()
    if params.get("targets")
)
# Specify whether to correct for efficiency bias (TRUE/FALSE)
EFFICIENCY_BIAS = "FALSE"
# Bin width for composition bias normalization (bp)
BIN_WIDTH = 10000
# Suffix for BAM files
SUFFIX = f".qc.sort.rmdup.mapq{config.MAPQ}.final.bam" if config.REMOVE_DUPLICATES else f".qc.sort.markdup.mapq{config.MAPQ}.final.bam"

### Set up directories
work_dir = os.path.join(config.SCALINGDIR, "tmm")
log_dir = os.path.join(work_dir, "log")

os.makedirs(config.SCALINGDIR, exist_ok=True)
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

# Paths to R scripts
GET_PE_SIZES_R = os.path.join(os.path.dirname(__file__), "utils", "getPEsizes.R")
SCALE_FACTORS_R = os.path.join(os.path.dirname(__file__), "utils", "scaleFactors.R")

def submit_tmm_job():
    """Creates and submits batch script for TMM normalization"""
    job_name = "tmm_normalization"
    bin_wd = int(BIN_WIDTH/1000) # Convert to kb for naming
    log_out = os.path.join(log_dir, f"tmm_{bin_wd}_kb_bins.log")
    log_err = os.path.join(log_dir, f"tmm_{bin_wd}_kb_bins.error")
    batch_file = os.path.join(log_dir, f"tmm_{bin_wd}_kb_bins.batch")

    # Build batch command
    batch_cmd = f"""#!/bin/bash
#BSUB -J "{job_name}"
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 8
#BSUB -M 64000
#BSUB -R "rusage [mem=8000] span[hosts=1]"

# Load required modules
source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Step 1: Get fragment sizes
echo "Running getPEsizes.R"
Rscript {GET_PE_SIZES_R}

# Step 2: Calculate TMM scaling factors
echo "Running scaleFactors.R with parameters:"
echo "TARGETS={TARGETS}"
echo "EFFICIENCY_BIAS={EFFICIENCY_BIAS}"
echo "BIN_WIDTH={BIN_WIDTH}"
Rscript {SCALE_FACTORS_R} "{TARGETS}" {EFFICIENCY_BIAS} {BIN_WIDTH} && echo "TMM normalization complete"

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    # Write batch file
    with open(batch_file, 'w', encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"Submitting TMM normalization job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    """Main function to submit job"""
    submit_tmm_job()

if __name__ == "__main__":
    main()
