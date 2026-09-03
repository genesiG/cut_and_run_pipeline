#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7d_classify_regions.py — Submit a batch job to classify genomic regions
into Retained and Lost based on spike-in normalised log2 fold-change.

Uses the generic helper `Scripts/utils/classify_regions_lfc.R`.

##############################################################################
# USER TOGGLES — edit these to change the run configuration
##############################################################################
"""

import argparse
import os
import sys
import time

import config

# ===========================================================================
# USER TOGGLES
# ===========================================================================
TARGET          = "K27me2"          # Antibody target label (must match SF file prefix)
                                    # e.g. "K27me2", "K27me3", "CBX2", "EZH2"

REGIONS         = (                 # Pre-specified BED regions to classify
    "Analysis_Data/peaks/csaw/"
    "H3K27me2_DMSO.w150.d50.filt2.w2000.d500.filt1.lfc1.merge50.tmm.bed"
    #"H3K27me3_DMSO.w150.d50.filt3.w2000.d500.filt1.5.lfc1.merge50.tmm.bed"
)

LFC_RETAINED    = -1.0              # logFC (1uM vs DMSO) ABOVE which → Retained

LFC_LOST        = -2.0              # logFC (1uM vs DMSO) BELOW which → Lost
                                    # (also requires FDR < FDR_THRESHOLD)

FDR_THRESHOLD   = 0.05               # FDR cutoff applied to Lost regions only
# ===========================================================================

# ---------------------------------------------------------------------------
# HPC resource defaults
# ---------------------------------------------------------------------------
NUM_CORES    = 4
MAX_MEM_MB   = 32768
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Submit bsub job: classify genomic regions into Retained vs Lost "
            "using spike-in normalised log2FC (generic step_7d)."
        )
    )
    parser.add_argument(
        "--target",
        default=TARGET,
        help=(
            "Antibody target label (default: %(default)s). "
            "Used to auto-discover metadata ({target}_processed.txt) and "
            "SF file ({target}_spikein_SF.txt), and to name all output files."
        ),
    )
    parser.add_argument(
        "--regions",
        default=REGIONS,
        metavar="PATH",
        help="Path to pre-specified BED regions to classify (default: %(default)s).",
    )
    parser.add_argument(
        "--sf_file",
        default=None,
        metavar="PATH",
        help=(
            "Path to spike-in SF file. "
            "Defaults to Analysis_Data/normalization/spikein/{TARGET}_spikein_SF.txt."
        ),
    )
    parser.add_argument(
        "--metadata",
        default=None,
        metavar="PATH",
        help=(
            "Path to sample metadata TSV (ID, ANTIBODY, GROUP columns). "
            "Defaults to Metadata/sample_metadata_{TARGET}_processed.txt."
        ),
    )
    parser.add_argument(
        "--lfc_retained",
        type=float,
        default=LFC_RETAINED,
        metavar="FLOAT",
        help="logFC (1uM vs DMSO) above which a region is Retained (default: %(default)s).",
    )
    parser.add_argument(
        "--lfc_lost",
        type=float,
        default=LFC_LOST,
        metavar="FLOAT",
        help="logFC (1uM vs DMSO) below which a region is Lost (default: %(default)s).",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=FDR_THRESHOLD,
        metavar="FLOAT",
        help="FDR threshold applied to Lost regions (default: %(default)s).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    target = args.target
    target_lc = target.lower()

    # -----------------------------------------------------------------------
    # Resolve defaults for optional paths
    # -----------------------------------------------------------------------
    sf_file  = args.sf_file  or os.path.join(
        config.ANALYSIS_DATA, "normalization", "spikein",
        f"{target}_spikein_SF.txt"
    )
    metadata = args.metadata or os.path.join(
        config.WORKDIR, "Metadata",
        f"sample_metadata_{target}_processed.txt"
    )

    # -----------------------------------------------------------------------
    # Output directory / files — all keyed on `target`
    # -----------------------------------------------------------------------
    out_dir    = os.path.join(config.ANALYSIS_DATA, f"{target_lc}_classification")
    os.makedirs(out_dir, exist_ok=True)

    job_name   = f"{target_lc}_classification"
    log_dir    = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    r_script = os.path.join(config.CODEDIR, "utils", "classify_regions_lfc.R")

    # -----------------------------------------------------------------------
    # Batch script
    # -----------------------------------------------------------------------
    batch_cmd = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM_MB}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

# ---- Run configuration (embedded for reproducibility) --------------------
# Target        : {target}
# Regions       : {args.regions}
# SF file       : {sf_file}
# Metadata      : {metadata}
# LFC retained  : > {args.lfc_retained}  (no FDR filter)
# LFC lost      : < {args.lfc_lost}  AND  FDR < {args.fdr}
# -------------------------------------------------------------------------

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Export RETICULATE_PYTHON so R's reticulate uses the correct interpreter
# without triggering the "use_python will be ignored" warning.
export RETICULATE_PYTHON=$(which python3)

ulimit -v unlimited 2>/dev/null || true

Rscript {r_script} \\
  --target        {target} \\
  --metadata      {metadata} \\
  --regions       {args.regions} \\
  --sf_file       {sf_file} \\
  --out_dir       {out_dir} \\
  --lfc_retained  {args.lfc_retained} \\
  --lfc_lost      {args.lfc_lost} \\
  --fdr           {args.fdr}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    # -----------------------------------------------------------------------
    # Submit
    # -----------------------------------------------------------------------
    print(f"\n{'='*60}")
    print(f"  step_7d: classify_regions_lfc")
    print(f"{'='*60}")
    print(f"  Target        : {target}")
    print(f"  Regions       : {args.regions}")
    print(f"  SF file       : {sf_file}")
    print(f"  Metadata      : {metadata}")
    print(f"  LFC retained  : > {args.lfc_retained}  (no FDR filter)")
    print(f"  LFC lost      : < {args.lfc_lost}  AND  FDR < {args.fdr}")
    print(f"  Memory        : {MAX_MEM_MB} MB ({NUM_CORES} cores × {MEM_PER_CORE} MB)")
    print(f"  Output dir    : {out_dir}")
    print(f"  Log           : {log_out}")
    print(f"  Batch script  : {batch_file}")
    print()

    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\nMonitor with: bjobs -J {job_name}")


if __name__ == "__main__":
    main()
