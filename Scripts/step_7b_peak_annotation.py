#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7b_peak_annotation.py — Submit a batch job to annotate nucleation/spreading
peak BED files using utils/peak_annotation.R.

The script:
  1. Discovers BED files under config.PEAKANNODIR/<caller>/ matching the
     supplied --glob-patterns (e.g. 'nucleation' 'spreading').
  2. Derives a label for each pattern (capitalised pattern, or --labels override).
  3. Builds and submits a bsub batch script that calls:
       Rscript utils/peak_annotation.R
         --bed_files <...>
         --labels    <...>
         --out_dir   <PEAKANNODIR/<caller>/>
         --workdir   <WORKDIR>
         --codedir   <CODEDIR>

Usage (run from project root on a login node — the heavy R runs inside bsub):
  python Scripts/step_7b_peak_annotation.py \\
      --caller   csaw \\
      --glob-patterns nucleation spreading

  # Override labels (order must match --glob-patterns):
  python Scripts/step_7b_peak_annotation.py \\
      --caller   csaw \\
      --glob-patterns nucleation spreading \\
      --labels   "Nucleation sites" "Spreading sites"

  # Use a custom search directory (overrides PEAKANNODIR/<caller>/ default):
  python Scripts/step_7b_peak_annotation.py \\
      --caller   csaw \\
      --search-dir Analysis_Data/k27me3_classification \\
      --glob-patterns retained lost \\
      --labels   "Retained H3K27me3" "Lost H3K27me3" \\
      --out-dir  Analysis_Data/k27me3_classification
"""

import argparse
import glob
import os
import sys
import time

import config

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
NUM_CORES    = 1
MAX_MEM_MB   = 32768
MEM_PER_CORE = MAX_MEM_MB


def parse_args():
    parser = argparse.ArgumentParser(
        description="Submit bsub job: annotate nucleation/spreading BED files."
    )
    parser.add_argument(
        "--caller",
        required=True,
        choices=["macs", "csaw", "seacr"],
        help="Peak caller used in step_7a (determines BED file search directory).",
    )
    parser.add_argument(
        "--glob-patterns",
        nargs="+",
        required=True,
        metavar="PATTERN",
        help=(
            "One or more filename glob patterns to find BED files inside "
            "config.PEAKANNODIR/<caller>/. "
            "E.g. 'nucleation' 'spreading'."
        ),
    )
    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        metavar="LABEL",
        help=(
            "Human-readable label for each glob pattern (order must match "
            "--glob-patterns). Defaults to pattern capitalised."
        ),
    )
    parser.add_argument(
        "--search-dir",
        default=None,
        metavar="PATH",
        help=(
            "Override the default BED file search directory "
            "(default: config.PEAKANNODIR/<caller>/). "
            "Useful when BED files are in a custom location such as a "
            "classification output directory."
        ),
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        metavar="PATH",
        help=(
            "Output directory for tables and plots. "
            "Defaults to config.PEAKANNODIR/<caller>/."
        ),
    )
    parser.add_argument(
        "--bigwig-files",
        nargs="+",
        default=None,
        metavar="BW",
        help=(
            "Paths to bigWig files. When supplied, EnrichedHeatmap-based TSS "
            "and gene-body heatmaps (plots 01a/01b) are generated using these "
            "signal tracks instead of the ChIPseeker geom_raster fallback."
        ),
    )
    parser.add_argument(
        "--bigwig-labels",
        nargs="+",
        default=None,
        metavar="LABEL",
        help=(
            "Labels for each bigWig file (order must match --bigwig-files). "
            "Defaults to the filename stem of each bigWig."
        ),
    )
    parser.add_argument(
        "--window-bp",
        type=int,
        default=3000,
        metavar="INT",
        help="Half-window in bp for TSS heatmap (default: 3000, i.e. ±3 kb).",
    )
    parser.add_argument(
        "--n-bins",
        type=int,
        default=100,
        metavar="INT",
        help="Bins per half-window for heatmap matrices (default: 100).",
    )
    return parser.parse_args()


def find_bed_file(search_dir: str, pattern: str) -> str:
    """
    Return the single BED file in search_dir whose name contains pattern
    (case-insensitive, anchored by '_sites.bed' or just pattern + '.bed').
    Raises SystemExit if zero or more than one match is found.
    """
    all_beds = glob.glob(os.path.join(search_dir, "*.bed"))
    matches  = [f for f in all_beds
                if pattern.lower() in os.path.basename(f).lower()]

    if not matches:
        print(
            f"ERROR: No BED file matching pattern '{pattern}' found in:\n"
            f"       {search_dir}\n"
            f"  Available files:\n"
            + "\n".join(f"    {os.path.basename(f)}" for f in sorted(all_beds)),
            file=sys.stderr,
        )
        sys.exit(1)

    if len(matches) > 1:
        print(
            f"ERROR: Multiple BED files match pattern '{pattern}' in:\n"
            f"       {search_dir}\n"
            f"  Matches:\n"
            + "\n".join(f"    {os.path.basename(f)}" for f in sorted(matches))
            + "\n  Refine your --glob-patterns to be more specific.",
            file=sys.stderr,
        )
        sys.exit(1)

    return matches[0]


def main():
    args = parse_args()

    caller         = args.caller
    patterns       = args.glob_patterns
    labels         = args.labels or [p.capitalize() for p in patterns]
    search_dir     = args.search_dir or os.path.join(config.PEAKANNODIR, caller)
    out_dir        = args.out_dir or search_dir
    bigwig_files   = args.bigwig_files or []
    bigwig_labels  = args.bigwig_labels or [
        os.path.splitext(os.path.basename(bw))[0] for bw in bigwig_files
    ]
    window_bp      = args.window_bp
    n_bins         = args.n_bins

    if len(labels) != len(patterns):
        print(
            f"ERROR: --labels ({len(labels)}) must match "
            f"--glob-patterns ({len(patterns)}) in count.",
            file=sys.stderr,
        )
        sys.exit(1)

    if bigwig_files and len(bigwig_files) != len(bigwig_labels):
        print(
            f"ERROR: --bigwig-labels ({len(bigwig_labels)}) must match "
            f"--bigwig-files ({len(bigwig_files)}) in count.",
            file=sys.stderr,
        )
        sys.exit(1)

    # --- Locate BED files ---------------------------------------------------
    print(f"\nSearching for BED files in: {search_dir}")
    bed_files = [find_bed_file(search_dir, p) for p in patterns]

    print("\nResolved:")
    for lbl, bed in zip(labels, bed_files):
        print(f"  [{lbl}]  {os.path.basename(bed)}")

    # --- Build bsub script --------------------------------------------------
    r_script = os.path.join(config.CODEDIR, "utils", "peak_annotation.R")

    bed_arg   = " \\\n    ".join(bed_files)
    label_arg = " \\\n    ".join(f'"{lbl}"' for lbl in labels)

    # Optional bigWig block appended only when files are supplied
    bigwig_block = ""
    if bigwig_files:
        bw_arg  = " \\\n    ".join(bigwig_files)
        bwl_arg = " \\\n    ".join(f'"{lbl}"' for lbl in bigwig_labels)
        bigwig_block = f""" \\
  --bigwig_files \\
    {bw_arg} \\
  --bigwig_labels \\
    {bwl_arg} \\
  --window_bp {window_bp} \\
  --n_bins    {n_bins}"""

    job_name   = f"peak_annotation.{caller}"
    log_dir    = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    batch_cmd = f"""
#!/bin/bash
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
  --bed_files \\
    {bed_arg} \\
  --labels \\
    {label_arg} \\
  --out_dir   {out_dir} \\
  --workdir   {config.WORKDIR} \\
  --codedir   {config.CODEDIR}{bigwig_block}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    # --- Submit job ---------------------------------------------------------
    print(f"\nSubmitting bsub job: {job_name}")
    print(f"  Memory : {MAX_MEM_MB} MB ({NUM_CORES} cores × {MEM_PER_CORE} MB)")
    print(f"  Output : {out_dir}")
    print(f"  Log    : {log_out}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\nDone. Monitor with: bjobs -J {job_name}")


if __name__ == "__main__":
    main()
