#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7l_correlate_occupancy.py — Submit bsub job that quantifies occupancy / enrichment
correlations between CBX targets, H3K27me2, and EZH2 across specified regions or bins.

Compares Ordinary Least Squares (OLS / lm) vs Generalized Additive Models (GAM / gam),
evaluates multiple linear regressions with Spatial HAC CIs and sr^2 variance decompositions,
and stratifies across EZH2 occupancy tertiles (including EZH2-absent/low regions).
"""

import argparse
import os
import sys
import time

import config

# ===========================================================================
# USER TOGGLES
# ===========================================================================
BIN_SIZE             = 10000
BIN_SIZE_SENSITIVITY = 2000
OUT_SUBDIR           = "occupancy_vs_occupancy/k27me2_dmso_peaks"
LFC_DIR              = "Analysis_Data/delta_vs_delta/k27me2_dmso_peaks"
OCC_TYPE             = "avelogcpm"
# ===========================================================================

NUM_CORES    = 2
MAX_MEM_MB   = 16384
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES

R_SCRIPT = os.path.join(os.path.dirname(__file__), "utils", "correlate_occupancy.R")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit bsub job: correlate target occupancy (step_7l)."
    )
    parser.add_argument("--bin_size",             type=int, default=BIN_SIZE,             help="Primary bin size (default: %(default)s).")
    parser.add_argument("--bin_size_sensitivity", type=int, default=BIN_SIZE_SENSITIVITY, help="Sensitivity bin size (default: %(default)s).")
    parser.add_argument("--out_subdir",           default=OUT_SUBDIR,                     help="Sub-directory under Analysis_Data for outputs (default: %(default)s).")
    parser.add_argument("--lfc_dir",              default=LFC_DIR,                        help="Path containing existing LFC/occupancy RDS files (default: %(default)s).")
    parser.add_argument("--occ_type",             default=OCC_TYPE,                       help="Occupancy metric ('avelogcpm' or 'dmso_logcpm', default: %(default)s).")
    parser.add_argument("--direct",               action="store_true", default=False,     help="Run directly without bsub.")
    return parser.parse_args()


def build_r_cmd(args, out_dir):
    lfc_dir_abs = os.path.abspath(args.lfc_dir)
    parts = [
        f"conda run -n chipseq Rscript {os.path.abspath(R_SCRIPT)}",
        f"  --bin_size             {args.bin_size}",
        f"  --bin_size_sensitivity {args.bin_size_sensitivity}",
        f"  --out_dir              {out_dir}",
        f"  --lfc_dir              {lfc_dir_abs}",
        f"  --occ_type             {args.occ_type}",
    ]
    return " \\\n".join(parts)


def run_direct(args):
    out_dir = os.path.join(config.ANALYSIS_DATA, args.out_subdir)
    os.makedirs(out_dir, exist_ok=True)
    r_cmd = build_r_cmd(args, out_dir)
    print(f"\n==============================================================")
    print("  step_7l_correlate_occupancy: running directly")
    print(f"==============================================================\n")
    ret = os.system(r_cmd)
    if ret != 0:
        sys.exit(f"Rscript failed with exit code {ret}")


def run_bsub(args):
    out_dir = os.path.join(config.ANALYSIS_DATA, args.out_subdir)
    os.makedirs(out_dir, exist_ok=True)
    log_dir = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)

    region_stem = os.path.basename(args.out_subdir.rstrip("/"))
    job_name = f"corr_occupancy_{args.bin_size}bp_{region_stem}"
    log_out = os.path.join(log_dir, f"{job_name}.log")
    log_err = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    r_cmd = build_r_cmd(args, out_dir)

    batch_cmd = f"""#!/bin/bash
#BSUB -J "{job_name}"
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM_MB}
#BSUB -q {config.BSUB_QUEUE}

cd {os.path.abspath(config.HOMEDIR + "/GG_EPICYPHER_CBX2")} && {r_cmd}
"""
    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"\nSubmitting step_7l_correlate_occupancy job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    args = parse_args()
    if args.direct:
        run_direct(args)
    else:
        run_bsub(args)


if __name__ == "__main__":
    main()
