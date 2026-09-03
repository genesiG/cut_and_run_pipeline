#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7c_inhibitor_resistance.py — Submit a batch job to run the inhibitor
resistance analysis and generate plots/heatmaps using step_7c_inhibitor_resistance.R.

The script:
  1. Builds and submits a bsub batch script that calls:
       Rscript Scripts/step_7c_inhibitor_resistance.R --caller <caller>
  2. Allocates appropriate memory (64 GB) to prevent OOM errors during
     EnrichedHeatmap and csaw edgeR matrix processing.

Usage (run from project root on a login node):
  python Scripts/step_7c_inhibitor_resistance.py --caller csaw
"""

import argparse
import os
import sys
import time

import config

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
NUM_CORES    = 1
MAX_MEM_MB   = 65536 // 2
MEM_PER_CORE = MAX_MEM_MB


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit bsub job: Inhibitor resistance DA and heatmaps."
    )
    parser.add_argument(
        "--caller",
        required=True,
        choices=["macs", "csaw", "seacr"],
        help="Peak caller used (csaw, macs, or seacr).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    caller = args.caller

    out_dir = os.path.join(config.ANALYSIS_DATA, "inhibitor_resistance", caller)
    os.makedirs(out_dir, exist_ok=True)

    # --- Build bsub script --------------------------------------------------
    r_script = os.path.join(config.CODEDIR, "utils", "inhibitor_resistance_analysis.R")

    job_name   = f"inhibitor_resistance.{caller}"
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

Rscript {r_script} --caller {caller}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    # --- Submit job ---------------------------------------------------------
    print(f"\\nSubmitting bsub job: {job_name}")
    print(f"  Memory : {MAX_MEM_MB} MB ({NUM_CORES} cores × {MEM_PER_CORE} MB)")
    print(f"  Output : {out_dir}")
    print(f"  Log    : {log_out}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\\nDone. Monitor with: bjobs -J {job_name}")


if __name__ == "__main__":
    main()
