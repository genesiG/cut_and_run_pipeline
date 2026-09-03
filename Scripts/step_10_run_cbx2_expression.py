#!/usr/bin/env python3

"""
step_10_run_cbx2_expression.py
Master Launcher for CBX2 Expression Analysis

Submits Phase 1 and Phase 2. Phase 3 (deepTools) is submitted in step_10a.
"""

import os
import time
import argparse
import config

NUM_CORES = 2
MAX_MEM_MB = 32768
MEM_PER_CORE = 16384

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", default="Analysis_Data/cbx2_expression")
    return parser.parse_args()

def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    log_dir = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    
    job_name = "cbx2_expression_prep_plot"
    log_out = os.path.join(log_dir, f"{job_name}.log")
    log_err = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")
    
    r_prep = os.path.join(config.CODEDIR, "utils", "cbx2_expression_prep.R")
    r_plots = os.path.join(config.CODEDIR, "utils", "cbx2_expression_plots.R")
    
    batch_cmd = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM_MB}
#BSUB -R "rusage [mem={MEM_PER_CORE}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq
export RETICULATE_PYTHON=$(which python3)

echo "=== Running Phase 1: Data Prep ==="
Rscript {r_prep} --out_dir {out_dir} --workdir {config.WORKDIR} --codedir {config.CODEDIR}

if [ $? -ne 0 ]; then
    echo "Phase 1 failed!"
    exit 1
fi

echo "=== Running Phase 2: Plots ==="
Rscript {r_plots} --out_dir {out_dir}

if [ $? -ne 0 ]; then
    echo "Phase 2 failed!"
    exit 1
fi

echo "=== Phase 1 & 2 Completed Successfully ==="
"""
    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)
        
    print(f"Submitting {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    
if __name__ == "__main__":
    main()
