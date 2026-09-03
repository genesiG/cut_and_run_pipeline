#!/usr/bin/env python3

"""
qc_2_multiqc.py - Aggregate FastQC outputs into a single MultiQC HTML report

Steps:
1. Define input (config.QCDIR1/fastqc) and output (config.QCDIR1/multiqc) directories.
2. Build a single bsub batch script that runs MultiQC on all FastQC outputs.
3. Submit the batch script to the cluster.
"""

### Import modules
import os
import time

import config
###

# Working directory and log directory
work_dir = os.path.join(config.QCDIR1, "multiqc")
log_dir  = os.path.join(work_dir, "log")

# FastQC input directory (output from step_0e_fastqc.py)
if config.USE_TRIMMOMATIC:
    fastqc_dir = os.path.join(config.QCDIR1, "trimmed")
else:
    fastqc_dir = os.path.join(config.QCDIR1, "untrimmed")

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir,  exist_ok=True)

NUM_CORES   = 4
MAX_MEM     = 16000
MEM_PER_CORE = MAX_MEM / NUM_CORES


def run_multiqc():
    """
    Creates and submits a batch script to:
    1. Run MultiQC on the FastQC output directory.
    2. Write the HTML report to the multiqc output directory.
    """

    job_name     = f"{config.PROJECT_NAME}.multiqc"
    log_out_file = os.path.join(log_dir, "multiqc.log")
    log_err_file = os.path.join(log_dir, "multiqc.error")
    batch_file   = os.path.join(log_dir, "multiqc.batch")

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running MultiQC on FastQC outputs in: {fastqc_dir}"

multiqc \\
    {fastqc_dir} \\
    --outdir {work_dir} \\
    --filename {config.PROJECT_NAME}_multiqc_report.html \\
    --force \\
    --verbose

echo "MultiQC report written to: {work_dir}/{config.PROJECT_NAME}_multiqc_report.html"

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"Submitting job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """
    Checks that the FastQC input directory exists, then submits
    a single MultiQC aggregation job.
    """
    if not os.path.isdir(fastqc_dir):
        print(f"WARNING: FastQC directory not found at {fastqc_dir}. "
              f"Run step_0e_fastqc.py first.")
        return

    run_multiqc()


if __name__ == "__main__":
    main()
