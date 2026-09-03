#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7i_compute_lfc_bins.py
----------------------------
Submits LSF batch job to compute genome-wide spike-in-normalized, shrunken
log2FC (DMSO → 1uM) across tiled genomic bins for a single antibody target.

USER TOGGLES (edit or pass as CLI arguments)
============================================
"""

import os
import argparse
import subprocess
import config

# ===========================================================================
# USER TOGGLES
# ===========================================================================
TARGET    = "K27me2"       # Target label; must match metadata ANTIBODY column
                            # and spikein SF filename prefix
                            # e.g. "K27me2", "CBX2", "GSTCBX7", "EZH2"
BIN_SIZE  = 10000           # Genomic bin width in bp (10000 or 2000)
REGIONS   = "Analysis_Data/peaks/csaw/H3K27me2_DMSO.w150.d50.filt2.w2000.d500.filt1.lfc1.merge100.tmm.bed"
FILTER_MIN_LOGCPM = 0.0 if REGIONS and REGIONS != "None" else 1.0  # Light coverage floor (aveLogCPM > 0.0) for regions
                            #   or aveLogCPM(counts) > 1.0 for tiled bins
OUT_SUBDIR = "delta_vs_delta"
# ===========================================================================


def parse_args():
    p = argparse.ArgumentParser(description="Genome-wide LFC bin/region computation (step 7i)")
    p.add_argument("--target",       default=TARGET,            help="Antibody target label")
    p.add_argument("--bin_size",     default=BIN_SIZE, type=int,help="Genomic bin width (bp)")
    p.add_argument("--filter_min",   default=FILTER_MIN_LOGCPM, type=float,
                   help="aveLogCPM threshold for bin filtering")
    p.add_argument("--regions",      default=REGIONS,           help="Pre-specified BED regions to count over via regionCounts")
    p.add_argument("--out_subdir",   default=OUT_SUBDIR,        help="Output subdirectory under Analysis_Data/")
    return p.parse_args()


def main():
    args = parse_args()

    target     = args.target
    bin_size   = args.bin_size
    filter_min = args.filter_min
    regions    = args.regions
    out_subdir = args.out_subdir

    out_dir  = os.path.join(config.ANALYSIS_DATA, out_subdir)
    log_dir  = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)

    target_safe = target.lower().replace("-", "_")
    tag         = "" if out_subdir == "delta_vs_delta" else "_" + os.path.basename(out_subdir.rstrip("/"))
    job_name    = f"lfc_bins_{target_safe}_{bin_size}bp{tag}"
    log_out     = os.path.join(log_dir, f"{job_name}.log")
    log_err     = os.path.join(log_dir, f"{job_name}.error")
    batch_file  = os.path.join(log_dir, f"{job_name}.batch")

    script_dir  = os.path.dirname(os.path.abspath(__file__))
    r_script    = os.path.join(script_dir, "utils", "compute_lfc_bins.R")

    r_cmd = (
        f"cd {config.WORKDIR} && "
        f"conda run -n chipseq Rscript {r_script} "
        f"--target {target} "
        f"--bin_size {bin_size} "
        f"--filter_min {filter_min} "
        f"--regions {regions} "
        f"--out_dir {out_dir}"
    )

    batch = f"""#!/bin/bash
#BSUB -J "{job_name}"
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -n 2
#BSUB -M 16384
#BSUB -q {config.BSUB_QUEUE}

{r_cmd}
"""

    if os.path.exists(log_out):
        os.remove(log_out)
    if os.path.exists(log_err):
        os.remove(log_err)

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch)

    print("=" * 62)
    print("  step_7i_compute_lfc_bins: submitting bsub job")
    print("=" * 62)
    print(f"  Job name      : {job_name}")
    print(f"  Target        : {target}")
    print(f"  Bin size      : {bin_size:,} bp")
    print(f"  Filter min    : aveLogCPM > {filter_min}")
    print(f"  Output dir    : {out_dir}")
    print(f"  Log           : {log_out}")
    print(f"  Batch script  : {batch_file}")
    print()

    result = subprocess.run(["bsub"], input=batch, capture_output=True, text=True, check=False)
    print(result.stdout.strip())
    if result.stderr.strip():
        print("STDERR:", result.stderr.strip())
    print(f"\nMonitor with: bjobs -J {job_name}")


if __name__ == "__main__":
    main()
