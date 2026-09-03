#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7j_correlate_deltas.py
----------------------------
Submits an LSF batch job to run correlate_deltas.R:
  - Merges per-target shrunken LFC bin tables (from step_7i)
  - Produces core delta-vs-delta scatterplots (CBX2 and CBX7 vs H3K27me2)
  - Runs partial correlation controlling for ΔEZH2 (ppcor)
  - Runs Steiger's Z test comparing CBX2 vs CBX7 correlation strengths
  - Stratifies by ΔEZH2 tertiles and re-checks Spearman correlations
  - Runs sensitivity check at 2kb bin size
"""

import os
import argparse
import subprocess
import config

# ===========================================================================
# USER TOGGLES
# ===========================================================================
BIN_SIZE_PRIMARY   = 10000   # Primary bin size (used for all analyses)
BIN_SIZE_SENSITIVITY = 2000  # Bin size for spatial sensitivity check
OUT_SUBDIR         = "delta_vs_delta"
# ===========================================================================


def parse_args():
    p = argparse.ArgumentParser(description="Delta-vs-delta correlation (step 7j)")
    p.add_argument("--bin_size",           default=BIN_SIZE_PRIMARY,     type=int)
    p.add_argument("--bin_size_sensitivity", default=BIN_SIZE_SENSITIVITY, type=int)
    p.add_argument("--out_subdir",         default=OUT_SUBDIR)
    p.add_argument("--lfc_type",           default="shrunk",             choices=["shrunk", "unshrunk"])
    p.add_argument("--lfc_dir",            default=None)
    return p.parse_args()


def main():
    args = parse_args()

    out_dir    = os.path.join(config.ANALYSIS_DATA, args.out_subdir)
    lfc_dir    = os.path.join(config.ANALYSIS_DATA, args.lfc_dir) if args.lfc_dir else out_dir
    log_dir    = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)

    tag        = "" if args.out_subdir == "delta_vs_delta" else "_" + os.path.basename(args.out_subdir.rstrip("/"))
    lfc_tag    = f"_{args.lfc_type}" if args.lfc_type != "shrunk" else ""
    job_name   = f"delta_vs_delta_{args.bin_size}bp{tag}{lfc_tag}"
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script   = os.path.join(script_dir, "utils", "correlate_deltas.R")

    r_cmd = (
        f"cd {config.WORKDIR} && "
        f"conda run -n chipseq Rscript {r_script} "
        f"--bin_size {args.bin_size} "
        f"--bin_size_sensitivity {args.bin_size_sensitivity} "
        f"--out_dir {out_dir} "
        f"--lfc_dir {lfc_dir} "
        f"--lfc_type {args.lfc_type}"
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
    print("  step_7j_correlate_deltas: submitting bsub job")
    print("=" * 62)
    print(f"  Job name            : {job_name}")
    print(f"  Primary bin size    : {args.bin_size:,} bp")
    print(f"  Sensitivity bin size: {args.bin_size_sensitivity:,} bp")
    print(f"  Output dir          : {out_dir}")
    print(f"  Log                 : {log_out}")
    print(f"  Batch script        : {batch_file}")
    print()

    result = subprocess.run(["bsub"], input=batch, capture_output=True, text=True, check=False)
    print(result.stdout.strip())
    if result.stderr.strip():
        print("STDERR:", result.stderr.strip())
    print(f"\nMonitor with: bjobs -J {job_name}")


if __name__ == "__main__":
    main()
