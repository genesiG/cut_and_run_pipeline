#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7d_k27me2_classification.py — Submit a batch job to classify H3K27me2 regions
into Retained and Lost using csaw and spike-in normalization factors.

Usage (run from project root on a login node):
  python Scripts/step_7d_k27me2_classification.py
"""

import argparse
import os
import sys
import time

import config

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
NUM_CORES    = 4
MAX_MEM_MB   = 32768
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES

def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit bsub job: H3K27me2 Retained vs Lost classification."
    )
    parser.add_argument(
        "--regions",
        default="Analysis_Data/peaks/csaw/H3K27me2_DMSO.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed",
        help="Path to pre-specified BED regions.",
    )
    parser.add_argument(
        "--sf_file",
        default="Analysis_Data/normalization/spikein/K27me2_spikein_SF.txt",
        help="Path to spike-in SF file.",
    )
    return parser.parse_args()

def main():
    args = parse_args()

    out_dir = os.path.join(config.ANALYSIS_DATA, "k27me2_classification")
    os.makedirs(out_dir, exist_ok=True)

    r_script = os.path.join(config.CODEDIR, "utils", "k27me2_lfc_spikein.R")

    job_name   = "k27me2_classification"
    log_dir    = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    batch_cmd = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM_MB}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq
export RETICULATE_PYTHON=$(which python3)

ulimit -v unlimited 2>/dev/null || true

Rscript {r_script} --regions {args.regions} --sf_file {args.sf_file} --out_dir {out_dir}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"\nSubmitting bsub job: {job_name}")
    print(f"  Memory : {MAX_MEM_MB} MB ({NUM_CORES} cores × {MEM_PER_CORE} MB)")
    print(f"  Output : {out_dir}")
    print(f"  Log    : {log_out}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\nDone. Monitor with: bjobs -J {job_name}")

if __name__ == "__main__":
    main()
