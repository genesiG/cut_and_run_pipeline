#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7h_correlate_logcpm.py — Submit a bsub job that correlates two antibody
targets (logCPM) over pre-specified genomic regions.

Uses the R helper: Scripts/utils/correlate_targets_logcpm.R

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
REGIONS = (                     # BED file defining the regions to count over
    "Analysis_Data/k27me3_classification/lost_k27me3.bed"
)

TARGET_A = "K27me2"             # First antibody target label (used to
                                #   auto-discover metadata file)

# Second antibody target label (e.g. "CBX2", "GSTCBX7", "EZH2")
TARGET_B = "CBX2"               

GROUP = "1uM"                   # Experimental group to compare
                                #   (must match GROUP column in metadata)

# Antibody values in the ANTIBODY column of each metadata file.
# If the metadata uses a different identifier than the target label, set it here.
# For CBX2:    ANTIBODY_B = "CBX2",    LABEL_B = "CBX2"
# For GSTCBX7: ANTIBODY_B = "GSTCBX7", LABEL_B = "GST-CBX7"
# For EZH2:    ANTIBODY_B = "EZH2",    LABEL_B = "EZH2"
ANTIBODY_A = "K27me2"          # ANTIBODY column value for target A
ANTIBODY_B = "CBX2"            # ANTIBODY column value for target B

# Display labels used on figure axes/titles
LABEL_A = "H3K27me2"
LABEL_B = "CBX2"

OUT_SUBDIR = "correlation"      # Subdirectory under Analysis_Data for outputs
# ===========================================================================

# ---------------------------------------------------------------------------
# HPC resource defaults
# ---------------------------------------------------------------------------
NUM_CORES    = 4
MAX_MEM_MB   = 32768
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES

R_SCRIPT = os.path.join(os.path.dirname(__file__),
                        "utils", "correlate_targets_logcpm.R")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Submit bsub job: correlate two antibody target logCPMs "
            "over pre-specified BED regions (step_7h)."
        )
    )
    parser.add_argument("--regions",    default=REGIONS,    metavar="PATH",
        help="BED regions to count reads over (default: %(default)s).")
    parser.add_argument("--target_a",   default=TARGET_A,   metavar="STR",
        help="First antibody target label (default: %(default)s).")
    parser.add_argument("--target_b",   default=TARGET_B,   metavar="STR",
        help="Second antibody target label (default: %(default)s).")
    parser.add_argument("--group",      default=GROUP,      metavar="STR",
        help="Group column value to select (default: %(default)s).")
    parser.add_argument("--antibody_a", default=ANTIBODY_A, metavar="STR",
        help="ANTIBODY column value for target A metadata (default: %(default)s).")
    parser.add_argument("--antibody_b", default=ANTIBODY_B, metavar="STR",
        help="ANTIBODY column value for target B metadata (default: %(default)s).")
    parser.add_argument("--label_a",    default=LABEL_A,    metavar="STR",
        help="Display label for target A in figures (default: %(default)s).")
    parser.add_argument("--label_b",    default=LABEL_B,    metavar="STR",
        help="Display label for target B in figures (default: %(default)s).")
    parser.add_argument("--metadata_a", default=None,       metavar="PATH",
        help="Metadata for target A (auto-discovered if omitted).")
    parser.add_argument("--metadata_b", default=None,       metavar="PATH",
        help="Metadata for target B (auto-discovered if omitted).")
    parser.add_argument("--out_subdir", default=OUT_SUBDIR, metavar="STR",
        help="Sub-directory under Analysis_Data for outputs (default: %(default)s).")
    parser.add_argument("--out_prefix", default=None,       metavar="STR",
        help="Output file prefix (auto-derived from targets+regions if omitted).")
    parser.add_argument("--direct",     action="store_true", default=False,
        help="Run Rscript directly without bsub (used internally by batch scripts).")
    return parser.parse_args()


def _safe(s: str) -> str:
    """Lower-case, replace non-alphanumeric with underscore."""
    return "".join(c if c.isalnum() else "_" for c in s).lower()


def build_r_cmd(args, out_dir):
    """Construct the Rscript command to pass all args to the R helper."""
    regions_abs = os.path.abspath(args.regions)

    meta_a = (args.metadata_a or
              os.path.join("Metadata",
                           f"sample_metadata_{args.target_a}_processed.txt"))
    meta_b = (args.metadata_b or
              os.path.join("Metadata",
                           f"sample_metadata_{args.target_b}_processed.txt"))

    # Auto-derive prefix if not provided
    regions_stem = os.path.splitext(os.path.basename(args.regions))[0]
    prefix = (args.out_prefix or
              f"{_safe(args.label_a)}_vs_{_safe(args.label_b)}"
              f"_{_safe(args.group)}_over_{_safe(regions_stem)}")

    parts = [
        f"Rscript {os.path.abspath(R_SCRIPT)}",
        f"  --regions    {regions_abs}",
        f"  --target_a   {args.target_a}",
        f"  --target_b   {args.target_b}",
        f"  --group      {args.group}",
        f"  --antibody_a {args.antibody_a}",
        f"  --antibody_b {args.antibody_b}",
        f"  --label_a    '{args.label_a}'",
        f"  --label_b    '{args.label_b}'",
        f"  --metadata_a {meta_a}",
        f"  --metadata_b {meta_b}",
        f"  --out_dir    {out_dir}",
        f"  --out_prefix {prefix}",
    ]
    return " \\\n".join(parts), prefix, meta_a, meta_b


def run_direct(args):
    out_dir = os.path.join(config.ANALYSIS_DATA, args.out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    r_cmd, prefix, meta_a, meta_b = build_r_cmd(args, out_dir)

    sep = "=" * 62
    print(f"\n{sep}")
    print("  step_7h_correlate_logcpm: running Rscript directly")
    print(sep)
    print(f"  Regions       : {args.regions}")
    print(f"  Target A      : {args.target_a} ({args.label_a})")
    print(f"  Target B      : {args.target_b} ({args.label_b})")
    print(f"  Group         : {args.group}")
    print(f"  Output dir    : {out_dir}")
    print()

    ret = os.system(r_cmd)
    if ret != 0:
        sys.exit(f"Rscript failed with exit code {ret}")
    print(f"\n  Done → {out_dir}/{prefix}_scatter.svg")


def run_bsub(args):
    out_dir = os.path.join(config.ANALYSIS_DATA, args.out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    log_dir = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)

    regions_stem = os.path.splitext(os.path.basename(args.regions))[0]
    job_name = (f"corr_{_safe(args.target_a)}_vs_{_safe(args.target_b)}"
                f"_{_safe(args.group)}_over_{_safe(regions_stem)}")

    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    r_cmd, prefix, meta_a, meta_b = build_r_cmd(args, out_dir)

    # Self-invocation with --direct so bsub script runs Rscript on compute node
    self_cmd = (
        f"python3 {os.path.abspath(__file__)}"
        f" --regions    {os.path.abspath(args.regions)}"
        f" --target_a   {args.target_a}"
        f" --target_b   {args.target_b}"
        f" --group      {args.group}"
        f" --antibody_a {args.antibody_a}"
        f" --antibody_b {args.antibody_b}"
        f" --label_a    '{args.label_a}'"
        f" --label_b    '{args.label_b}'"
        f" --metadata_a {meta_a}"
        f" --metadata_b {meta_b}"
        f" --out_subdir {args.out_subdir}"
        f"{' --out_prefix ' + args.out_prefix if args.out_prefix else ''}"
        f" --direct"
    )

    batch_cmd = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM_MB}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

# ---- Run configuration (embedded for reproducibility) ----------------------
# Regions    : {args.regions}
# Target A   : {args.target_a} ({args.label_a})  antibody={args.antibody_a}
# Target B   : {args.target_b} ({args.label_b})  antibody={args.antibody_b}
# Group      : {args.group}
# Metadata A : {meta_a}
# Metadata B : {meta_b}
# Norm       : library-size (CPM) only -- no spike-in SFs
# Output     : {out_dir}
# ---------------------------------------------------------------------------

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

export RETICULATE_PYTHON=$(which python3)

ulimit -v unlimited 2>/dev/null || true

{self_cmd}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    sep = "=" * 62
    print(f"\n{sep}")
    print("  step_7h_correlate_logcpm: submitting bsub job")
    print(sep)
    print(f"  Job name      : {job_name}")
    print(f"  Regions       : {args.regions}")
    print(f"  Target A      : {args.target_a} ({args.label_a})")
    print(f"  Target B      : {args.target_b} ({args.label_b})")
    print(f"  Group         : {args.group}")
    print(f"  Output dir    : {out_dir}")
    print(f"  Log           : {log_out}")
    print(f"  Batch script  : {batch_file}")
    print()

    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\nMonitor with: bjobs -J {job_name}")


def main():
    args = parse_args()

    if not os.path.exists(args.regions):
        sys.exit(f"Regions BED not found: {args.regions}")

    if args.direct:
        run_direct(args)
    else:
        run_bsub(args)


if __name__ == "__main__":
    main()
