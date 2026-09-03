#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7a_k27me3_classification.py — Submit a batch job to classify H3K27me3
regions into Retained and Lost using set operations on peak BED files.

Logic (mirrors step_7a_define_nucleation_spreading.R):
  - Retained H3K27me3 : reduce(H3K27me3_1uMcpd peaks)
                         → peaks that persist under EZH2 inhibition
  - Lost H3K27me3     : setdiff(reduce(DMSO peaks), 1uM peaks), regions >= 200 bp
                         → DMSO-only peaks with no overlap to any 1uM peak

No spike-in normalisation or edgeR is used; classification is purely
based on peak presence / absence across the two conditions.

Usage (run from project root on a login node):

  # Auto-discover via csaw peak directory:
  python Scripts/step_7a_k27me3_classification.py

  # Supply explicit BED file paths:
  python Scripts/step_7a_k27me3_classification.py \\
      --bed-dmso Analysis_Data/peaks/csaw/H3K27me3_DMSO.w150.d50.filt3.w2000.d500.filt1.5.lfc1.merge100.tmm.bed \\
      --bed-1um  Analysis_Data/peaks/csaw/H3K27me3_1uMcpd.w150.d50.filt3.w2000.d500.filt1.5.lfc1.merge100.tmm.bed

Outputs: Analysis_Data/k27me3_classification/
  retained_H3K27me3.bed
  lost_H3K27me3.bed
"""

import argparse
import os
import time

import config

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
NUM_CORES    = 2
MAX_MEM_MB   = 16384
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES

DEFAULT_BED_DMSO = (
    "Analysis_Data/peaks/csaw/"
    "H3K27me3_DMSO.w150.d50.filt3.w2000.d500.filt2.lfc1.merge500.tmm.bed"
)
DEFAULT_BED_1UM = (
    "Analysis_Data/peaks/csaw/"
    "H3K27me3_1uMcpd.w150.d50.filt3.w2000.d500.filt2.lfc1.merge500.tmm.bed"
)
DEFAULT_BED_K27ME2 = (
    "Analysis_Data/peaks/csaw/"
    "H3K27me2_1uMcpd.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit bsub job: H3K27me3 Retained vs Lost classification (set-based)."
    )
    parser.add_argument(
        "--peak-caller",
        default="csaw",
        choices=["csaw", "macs", "seacr"],
        help=(
            "Peak caller to use when auto-discovering BED files "
            "(default: csaw). Ignored when --bed-dmso/--bed-1um are supplied."
        ),
    )
    parser.add_argument(
        "--bed-dmso",
        default=DEFAULT_BED_DMSO,
        metavar="PATH",
        help="Explicit path to the H3K27me3 DMSO csaw BED file.",
    )
    parser.add_argument(
        "--bed-1um",
        default=DEFAULT_BED_1UM,
        metavar="PATH",
        help="Explicit path to the H3K27me3 1uMcpd csaw BED file.",
    )
    parser.add_argument(
        "--bed-k27me2",
        default=DEFAULT_BED_K27ME2,
        metavar="PATH",
        help="Explicit path to the K27me2 csaw BED file for intersection with lost K27me3.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    out_dir = os.path.join(config.ANALYSIS_DATA, "k27me3_classification")
    os.makedirs(out_dir, exist_ok=True)

    r_script = os.path.join(config.CODEDIR, "utils", "k27me3_classification.R")

    job_name   = "k27me3_classification"
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

Rscript {r_script} \\
  --peak-caller {args.peak_caller} \\
  --bed-dmso    {args.bed_dmso} \\
  --bed-1um     {args.bed_1um} \\
  --bed-k27me2  {args.bed_k27me2}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"\nSubmitting bsub job: {job_name}")
    print(f"  Logic        : Retained = reduce(1uM peaks); Lost = DMSO setdiff(1uM), >=200 bp")
    print(f"  Peak caller  : {args.peak_caller}")
    print(f"  BED DMSO     : {args.bed_dmso}")
    print(f"  BED 1uM      : {args.bed_1um}")
    print(f"  BED K27me2   : {args.bed_k27me2}")
    print(f"  Memory       : {MAX_MEM_MB} MB ({NUM_CORES} cores × {MEM_PER_CORE} MB)")
    print(f"  Output       : {out_dir}")
    print(f"  Log          : {log_out}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\nDone. Monitor with: bjobs -J {job_name}")


if __name__ == "__main__":
    main()
