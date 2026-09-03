#!/usr/bin/env python3

"""
step_4_qc.py  –  Submit a bsub job that generates normalization QC plots.

Runs Scripts/utils/qc_normalization_plots.R on the cluster with elevated
cores and memory, producing MA plots, PCA, and Spearman correlation heatmaps
for every normalization run defined in config.NORMALIZATION_RUNS.

Steps:
1. Set up log directory inside the normalization plots output folder.
2. Build a batch script pointing at the utility R script.
3. Submit via bsub and sleep briefly.

Usage:
    python3 Scripts/step_4_qc.py

Outputs: Analysis_Data/normalization/qc/
Logs:    Analysis_Data/normalization/qc/log/
"""

### Import modules
import os
import time
import config
###

### Constants
NUM_CORES    = 4
MAX_MEM      = 64000
MEM_PER_CORE = MAX_MEM // NUM_CORES

### Paths
work_dir   = os.path.join(config.SCALINGDIR, "qc")
log_dir    = os.path.join(work_dir, "log")
SCRIPT_PATH = os.path.join(config.CODEDIR, "utils", "qc_normalization_plots.R")

os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir,  exist_ok=True)
###


def submit_qc_job() -> None:
    """Build and submit a single bsub job that runs qc_normalization_plots.R."""
    job_name   = f"{config.PROJECT_NAME}.step4_qc"
    log_out    = os.path.join(log_dir, "step_4_qc.log")
    log_err    = os.path.join(log_dir, "step_4_qc.error")
    batch_file = os.path.join(log_dir, "step_4_qc.batch")

    if not os.path.exists(SCRIPT_PATH):
        raise FileNotFoundError(
            f"Utility script not found: {SCRIPT_PATH}\n"
            "Run step_4_qc.py only after Scripts/utils/qc_normalization_plots.R exists."
        )

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Run from project root so config.py resolves all relative paths correctly.
cd {config.WORKDIR}

echo "=== step_4_qc: normalization comparison plots ==="
echo "Project : {config.PROJECT_NAME}"
echo "Script  : {SCRIPT_PATH}"
echo "Output  : {work_dir}"
echo "Started : $(date)"
echo ""

Rscript {SCRIPT_PATH}

echo ""
echo "Finished: $(date)"
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    print(f"Submitting job: {job_name}")
    print(f"  Script : {SCRIPT_PATH}")
    print(f"  Output : {work_dir}")
    print(f"  Log    : {log_out}")
    print(f"  Error  : {log_err}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main() -> None:
    """Submit the normalization QC bsub job."""
    submit_qc_job()


if __name__ == "__main__":
    main()
