#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_6b_seacr.py  —  SEACR peak calling for H3K27me3 CUT&TAG samples.

Strategy
--------
For each replicate in a group (K27M, K27MKO):
  A. Find treatment BAM.
  B. If using control BAM (default):
       Find matching IgG control BAM.
       Convert both BAMs to fragment-level bedgraph.
       Run SEACR_1.3.sh with "norm" mode.
     If using a numeric threshold:
       Convert treatment BAM to fragment-level bedgraph.
       Run SEACR_1.3.sh with the numeric threshold and "non" mode.
"""

import os
import glob
import time
import re
import argparse

import config

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SEACR_SCRIPT   = os.path.join(config.CONDA_ENV_BIN, "SEACR_1.3.sh")
work_dir       = config.SEACR_OUTDIR
bedgraph_dir   = config.BEDGRAPHDIR
log_dir        = os.path.join(work_dir, "log")
processed_dir  = config.PROCESSEDBAMDIR

os.makedirs(work_dir,      exist_ok=True)
os.makedirs(bedgraph_dir,  exist_ok=True)
os.makedirs(log_dir,       exist_ok=True)

NUM_CORES    = 2
MAX_MEM      = 16000
MEM_PER_CORE = MAX_MEM // NUM_CORES

BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if config.REMOVE_DUPLICATES
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

PEAK_GROUPS = [
    ("K27M_DMSO",    "K27M_DMSO_K27me3",    "K27M_DMSO_IgG"),
    ("K27M_500nM",   "K27M_500nM_K27me3",   "K27M_500nM_IgG"),
    ("K27MKO_DMSO",  "K27MKO_DMSO_K27me3",  "K27MKO_DMSO_IgG"),
    ("K27MKO_500nM", "K27MKO_500nM_K27me3", "K27MKO_500nM_IgG"),
]

def clean_sample(sample_name):
    return re.sub(r"\.qc\..*$", "", sample_name)

def find_bams(pattern):
    return sorted(glob.glob(os.path.join(processed_dir, f"*{pattern}*{BAM_SUFFIX}")))

def parse_args():
    parser = argparse.ArgumentParser(description="Run SEACR peak calling.")
    parser.add_argument("--threshold", type=float, default=None,
                        help="Numeric threshold between 0 and 1 to replace IgG control bedgraph in SEACR.")
    return parser.parse_args()

def run_seacr_replicate(treat_bam, ctrl_bam=None, threshold=None):
    sample_name = clean_sample(os.path.basename(treat_bam))
    job_name   = f"{sample_name}.seacr"
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    treat_bg      = os.path.join(bedgraph_dir, f"{sample_name}_treat.bedgraph")
    treat_bg_sort = os.path.join(bedgraph_dir, f"{sample_name}_treat.sorted.bedgraph")
    if threshold is not None:
        seacr_prefix  = os.path.join(work_dir, f"{sample_name}_t{threshold}")
        ctrl_block = ""
        seacr_ctrl_arg = str(threshold)
        seacr_norm_arg = "non"
        cleanup_files = f"{treat_bg} {treat_bg_sort}"
    else:
        seacr_prefix  = os.path.join(work_dir, f"{sample_name}")
        ctrl_bg       = os.path.join(bedgraph_dir, f"{sample_name}_ctrl.bedgraph")
        ctrl_bg_sort  = os.path.join(bedgraph_dir, f"{sample_name}_ctrl.sorted.bedgraph")
        ctrl_block = f"""\
bedtools genomecov \\
    -ibam {ctrl_bam} \\
    -bg -pc \\
    -g {config.GNM_SIZES} \\
    > {ctrl_bg}

sort -k1,1 -k2,2n {ctrl_bg}  > {ctrl_bg_sort}
"""
        seacr_ctrl_arg = ctrl_bg_sort
        seacr_norm_arg = "norm"
        cleanup_files = f"{treat_bg} {ctrl_bg} {treat_bg_sort} {ctrl_bg_sort}"

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

echo '=== Step A: BAM to bedgraph ==='

bedtools genomecov \\
    -ibam {treat_bam} \\
    -bg -pc \\
    -g {config.GNM_SIZES} \\
    > {treat_bg}

sort -k1,1 -k2,2n {treat_bg} > {treat_bg_sort}

{ctrl_block}

echo 'Bedgraphs created and sorted.'

echo '=== Step B: SEACR peak calling ==='

bash {SEACR_SCRIPT} \\
    {treat_bg_sort} \\
    {seacr_ctrl_arg} \\
    {seacr_norm_arg} \\
    relaxed \\
    {seacr_prefix}

echo 'SEACR complete. Output files:'
ls -lh {seacr_prefix}.stringent.bed {seacr_prefix}.relaxed.bed 2>/dev/null || \\
    echo 'WARNING: expected output files not found — check SEACR log above.'

echo '=== Step C: Cleanup ==='
rm -f {cleanup_files}

echo '=== Job complete ==='
echo '=== Resource usage summary ==='
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, 'w', encoding='utf-8') as fh:
        fh.write(batch_cmd)

    print(f"Submitting: {job_name}")
    print(f"  Treatment: {os.path.basename(treat_bam)}")
    if threshold is not None:
        print(f"  Control  : Numeric threshold {threshold}")
    else:
        print(f"  Control  : {os.path.basename(ctrl_bam)}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    args = parse_args()

    if args.threshold is not None:
        if not (0 <= args.threshold <= 1):
            print("ERROR: Threshold must be between 0 and 1.")
            return

    for group_label, treat_pattern, ctrl_pattern in PEAK_GROUPS:
        treatment_bams = find_bams(treat_pattern)

        if not treatment_bams:
            print(f"WARNING: No BAMs found for '{treat_pattern}' — skipping {group_label}.")
            continue

        for treat_bam in treatment_bams:
            if args.threshold is not None:
                run_seacr_replicate(treat_bam, threshold=args.threshold)
            else:
                bam_name = os.path.basename(treat_bam)
                ctrl_bam_name = bam_name.replace(treat_pattern, ctrl_pattern)
                ctrl_bam = os.path.join(processed_dir, ctrl_bam_name)

                if not os.path.exists(ctrl_bam):
                    print(f"ERROR: Missing control BAM {ctrl_bam} for {treat_bam}")
                    continue

                run_seacr_replicate(treat_bam, ctrl_bam=ctrl_bam)

if __name__ == "__main__":
    main()
