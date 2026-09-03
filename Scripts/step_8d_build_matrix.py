#!/usr/bin/env python3

"""
step_8d_build_matrix.py — Phase 1: build normalizedMatrix objects

Submits one bsub batch job per window size. Each job calls:
  MC_CORES=N Rscript Scripts/utils/build_matrix.R ...

The R script uses mclapply (fork-based, true Linux parallelism) to import
all bigWig files concurrently, one per core, and saves:
  {prefix}_norm_mats.rds         — normalizedMatrix list (one per sample)
  {prefix}_targets.rds           — original GRanges peak regions
  {prefix}_targets_center.rds    — 1-bp centred GRanges
  {prefix}_cluster_matrix.tsv.gz — H3K27me2 DMSO + 1uM columns for step_8e

A JSON config is also written to {out_dir}/log/step_8d_config.json;
steps 8e and 8f read it to discover paths without duplicating parameters.

Usage (from project root):
  python3 Scripts/step_8d_build_matrix.py
"""

### Import modules
import json
import os
import sys
import time

import config
###


# =============================================================================
# ===  USER CONFIGURATION — edit ONLY this section  ===========================
# =============================================================================

# --- BigWig directory and filenames ---
# Separate BW_DIR from filenames so relocating the bigwig folder requires
# changing only one variable.
BW_DIR = config.BIGWIGDIR

BW_FILENAMES = [
    "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw",
    "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw",
    "RP_036_H3K27me2_DMSO_R1.spikein.rpkm.bw",
    "RP_060_H3K27me2_1uMcpd_R1.spikein.rpkm.bw",
    "RP_038_H3K27me3_DMSO_R1.spikein.rpkm.bw",
    "RP_062_H3K27me3_1uMcpd_R1.spikein.rpkm.bw",
    "RP_042_EZH2_DMSO_R1.spikein.rpkm.bw",
    "RP_066_EZH2_1uMcpd_R1.spikein.rpkm.bw",
]

# Internal sample labels (keys into norm_mats, must be R-safe identifiers)
BW_LABELS = [
    "CBX2_DMSO",    "CBX2_1uM",
    "H3K27me2_DMSO","H3K27me2_1uM",
    "H3K27me3_DMSO","H3K27me3_1uM",
    "EZH2_DMSO",    "EZH2_1uM",
]

# Full absolute paths — derived automatically, do not edit
BW_FILES = [os.path.join(BW_DIR, f) for f in BW_FILENAMES]

# --- Genomic regions ---
REGIONS_BED = os.path.join(
    config.ANALYSIS_DATA, "csaw", "peaks",
    "AntiCBX2RP_1uMcpd.w150.d50.filt2.w2000.d500.filt1.5.lfc1.merge100.tmm.bed"
)

# --- Window sizes (one batch job submitted per entry) ---
# Each entry is a half-window in bp (signal is extracted +/- WINDOW_BP from
# each peak centre). Edit the list to restrict to a single window if needed.
WINDOW_CONFIGS = [10_000, 5_000, 2_000]  # 10 kb, 5 kb, 2 kb

# --- Binning ---
N_BINS = 40  # bins per half-window; bin_size = WINDOW_BP / N_BINS

# --- Output ---
BASE_PREFIX   = "cbx2_retained"
OUT_DIR_BASE  = os.path.join(config.HEATMAPDIR, "cbx2_retained")

# --- LSF resources ---
NUM_CORES    = 4       # must match len(BW_FILENAMES) for full parallelism
MAX_MEM_MB   = 65_536  # 64 GB total
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES
QUEUE        = "normal"
WALL_TIME    = "04:00"

# =============================================================================


def submit_window(window_bp: int):
    """
    Write the JSON config and batch script, then submit via bsub,
    for a single window size.
    """
    half_kb    = window_bp // 1_000
    prefix     = f"{BASE_PREFIX}_{half_kb}kb"
    out_dir    = os.path.join(OUT_DIR_BASE, f"{half_kb}kb")
    log_dir    = os.path.join(out_dir, "log")
    job_name   = f"buildmat_{half_kb}kb"
    log_out    = os.path.join(log_dir, f"{prefix}.log")
    log_err    = os.path.join(log_dir, f"{prefix}.error")
    batch_file = os.path.join(log_dir, f"{prefix}_build_matrix.batch")
    config_json = os.path.join(log_dir, "step_8d_config.json")

    # Create directories on the login node (lightweight, no compute required)
    os.makedirs(log_dir, exist_ok=True)

    # --- Validate inputs before submitting ---
    errors = []
    for bw in BW_FILES:
        if not os.path.exists(bw):
            errors.append(f"  bigWig not found: {bw}")
    if not os.path.exists(REGIONS_BED):
        errors.append(f"  Regions BED not found: {REGIONS_BED}")
    if errors:
        print(f"[{prefix}] ERROR: missing inputs:")
        for e in errors:
            print(e)
        sys.exit(1)

    # --- Write JSON config (read by steps 8e and 8f) ---
    cfg_data = {
        "window_bp":           window_bp,
        "n_bins":              N_BINS,
        "out_dir":             out_dir,
        "prefix":              prefix,
        "bw_dir":              BW_DIR,
        "bw_filenames":        BW_FILENAMES,
        "bw_labels":           BW_LABELS,
        "bw_files":            BW_FILES,
        "regions_bed":         REGIONS_BED,
    }
    with open(config_json, "w", encoding="utf-8") as fh:
        json.dump(cfg_data, fh, indent=2)
    print(f"  Config written: {config_json}")

    # --- Build batch script ---
    r_script   = os.path.join(config.CODEDIR, "utils", "build_matrix.R")
    bw_files_r = " ".join(BW_FILES)
    bw_labels_r = " ".join(BW_LABELS)

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM_MB}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"
#BSUB -W {WALL_TIME}

# ── Memory rationale ──────────────────────────────────────────────────────────
# {NUM_CORES} R child processes run in parallel (mclapply / fork).
# Each child imports one bigWig subset (~90 MB) and builds a normalizedMatrix
# Total request: {MAX_MEM_MB // 1024} GB ({MEM_PER_CORE} MB/slot × {NUM_CORES} cores)
# ─────────────────────────────────────────────────────────────────────────────

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# MC_CORES tells build_matrix.R how many fork workers to spawn.
# Each worker imports one bigWig file — genuine Linux fork() parallelism.
export MC_CORES={NUM_CORES}

Rscript {r_script} \\
  --regions_bed {REGIONS_BED} \\
  --bigwig_files {bw_files_r} \\
  --bigwig_labels {bw_labels_r} \\
  --window_bp {window_bp} \\
  --n_bins {N_BINS} \\
  --out_dir {out_dir} \\
  --out_prefix {prefix}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"  Submitting {job_name}  (window={half_kb}kb, cores={NUM_CORES})")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """Submit one Phase-1 job per window configuration."""
    print("=== step_8d_build_matrix.py ===")
    print(f"BigWig directory  : {BW_DIR}")
    print(f"Regions BED       : {REGIONS_BED}")
    print(f"Window configs    : {[f'{w//1000}kb' for w in WINDOW_CONFIGS]}")
    print(f"Output base dir   : {OUT_DIR_BASE}")
    print()

    for window_bp in WINDOW_CONFIGS:
        submit_window(window_bp)


if __name__ == "__main__":
    main()
