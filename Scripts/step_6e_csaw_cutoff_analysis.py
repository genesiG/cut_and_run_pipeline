#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_5c_opt_callpeaks_csaw.py - Cutoff-optimisation wrapper for csaw

Mirrors step_5c_callpeaks_csaw.py but calls optimize_csaw_cutoffs.R instead
of callpeaks_csaw.R.  For each (sample, target) pair it submits a bsub batch
job that performs an 11 x 9 grid search to find the (c_small, c_large) pair
that maximises significant DBRs at the configured FDR threshold.

Outputs per sample (in Analysis_Data/csaw/cutoff_opt/{sample_clean}/):
  *_cutoff_results.csv        — full grid results table
  *_pareto_heatmap.pdf        — 2-D Pareto heatmap
  *_opt_filter_histograms.pdf — filter histogram at optimal cutoff
  *_opt_bcv.pdf               — BCV plot at optimal cutoff
  *_opt.bed                   — peaks at optimal cutoff only
"""

import os
import time
import glob
import re
import config

### Configuration
BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if config.REMOVE_DUPLICATES
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

CTRL_LABEL = "IgG" if config.EXPERIMENT in ("cutandrun", "cutandtag") else "Input"

CSAWDIR  = getattr(config, "CSAW_OUTDIR", os.path.join(config.PEAKDIR, "csaw"))
CUTOFF_DIR = os.path.join(CSAWDIR, "cutoff_analysis")
LOG_DIR  = os.path.join(CUTOFF_DIR, "log")
os.makedirs(LOG_DIR, exist_ok=True)
###


def clean_sample(sample: str) -> str:
    """Strip trailing underscores / spaces from a sample name for file naming."""
    return re.sub(r"[_\s]+$", "", sample)


def find_bams(pattern: str) -> list[str]:
    return sorted(glob.glob(os.path.join(config.PROCESSEDBAMDIR, pattern)))


def submit_opt_job(sample: str, target: str,
                   chip_bams: list[str], ctrl_bams: list[str]) -> None:
    sample_clean = clean_sample(sample)
    sample_id    = f"{sample_clean}_{target}"
    job_name     = f"{sample_id}_csaw_opt"

    log_out    = os.path.join(LOG_DIR, f"{job_name}.log")
    log_err    = os.path.join(LOG_DIR, f"{job_name}.error")
    batch_file = os.path.join(LOG_DIR, f"{job_name}.batch")

    r_script       = os.path.join(config.CODEDIR, "utils", "optimize_csaw_cutoffs.R")
    chip_bams_str  = " ".join(chip_bams)
    ctrl_bams_str  = " ".join(ctrl_bams)

    # 8 cores as agreed — BiocParallel uses all available workers
    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 1
#BSUB -M 131072
#BSUB -R "rusage [mem=131072] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq
export RETICULATE_PYTHON=$(which python3)

ulimit -v unlimited 2>/dev/null || true

Rscript {r_script} \\
  --chip_bams {chip_bams_str} \\
  --ctrl_bams {ctrl_bams_str} \\
  --sample_id {sample_id} \\
  --workdir   {config.WORKDIR} \\
  --codedir   {config.CODEDIR} \\
  --norm_method TMM

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as f:
        f.write(batch_cmd)

    print(f"\nSubmitting {job_name}")
    print(f"  ChIP BAMs ({len(chip_bams)}):")
    for b in chip_bams:
        print(f"    {os.path.basename(b)}")
    print(f"  Ctrl BAMs ({len(ctrl_bams)}):")
    for b in ctrl_bams:
        print(f"    {os.path.basename(b)}")

    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    samples = config.SAMPLES
    targets = config.ANTIBODY_TARGETS

    for sample in samples:
        for target in targets:
            print(f"\n=== Processing: {sample} / {target} ===")

            chip_bams = find_bams(f"*{sample}*{target}*{BAM_SUFFIX}")
            if not chip_bams:
                print("  WARNING: No ChIP BAMs found — skipping")
                continue

            ctrl_bams = find_bams(f"*{sample}*{CTRL_LABEL}*{BAM_SUFFIX}")
            if not ctrl_bams:
                print(f"  WARNING: No {CTRL_LABEL} BAMs found — skipping")
                continue

            submit_opt_job(sample, target, chip_bams, ctrl_bams)


if __name__ == "__main__":
    main()
