#!/usr/bin/env python3

"""
step_8f_plot_heatmaps.py — Phase 3: plot EnrichedHeatmaps

Submits one bsub batch job per window size. Each job calls:
  Rscript Scripts/utils/plot_enriched_heatmaps_clustered.R <config.json>

The R script reads all parameters from a JSON file written by this script,
which merges the step_8d matrix config with the plotting parameters defined
below. This avoids shell-quoting issues with multi-line display labels.

Reads:
  {out_dir}/log/step_8d_config.json        (written by step_8d)
  {out_dir}/{prefix}_norm_mats.rds         (written by step_8d Phase 1)
  {out_dir}/{prefix}_targets.rds           (written by step_8d Phase 1)
  {out_dir}/cluster_assignments_k{K}.csv   (written by step_8e Phase 2)

Writes:
  {out_dir}/log/step_8f_plot_config.json              (combined config)
  {out_dir}/{prefix}_kmeans_k{K}_heatmap.svg          (main heatmap)
  {out_dir}/regions/{prefix}_k{K}_{cluster_label}.bed (per-cluster BEDs)

Usage (from project root):
  python3 Scripts/step_8f_plot_heatmaps.py
"""

### Import modules
import json
import os
import time

import config
###


# =============================================================================
# ===  USER CONFIGURATION — edit ONLY this section  ===========================
# =============================================================================

# --- Sample display configuration ---
# SAMPLE_LABELS must match the BW_LABELS defined in step_8d.
# SAMPLE_DISPLAY_LABELS can contain "\\n" for line breaks in the plot title.
SAMPLE_LABELS = [
    "CBX2_DMSO",    "CBX2_1uM",
    "H3K27me2_DMSO","H3K27me2_1uM",
    "H3K27me3_DMSO","H3K27me3_1uM",
    "EZH2_DMSO",    "EZH2_1uM",
]

SAMPLE_DISPLAY_LABELS = [
    "CBX2\nDMSO",      "CBX2\n1\u00b5M",
    "H3K27me2\nDMSO",  "H3K27me2\n1\u00b5M",
    "H3K27me3\nDMSO",  "H3K27me3\n1\u00b5M",
    "EZH2\nDMSO",      "EZH2\n1\u00b5M",
]

# --- Colors (from step_8c color_list) ---
# white → saturated hue ramp per sample.
SAMPLE_HUES = [
    "#58135e", "#58135e",   # CBX2 DMSO / 1uM
    "#8c96c6", "#8c96c6",   # H3K27me2 DMSO / 1uM
    "#8c6bb1", "#8c6bb1",   # H3K27me3 DMSO / 1uM
    "#39107b", "#39107b",   # EZH2 DMSO / 1uM
]

# --- Color scale limits ---
# Set to None for data-driven scaling (cap = max column mean).
Z_MIN = [0,  0,  0,  0,  0,   0,  0,  0]
Z_MAX = [15, 15, 52, 52, 190, 190, 35, 35]

# --- Cluster labels and colors ---
# Applied in cluster-ID order (1 → first label, 2 → second label).
# Inspect the output heatmap and swap if the assignment is inverted.
CLUSTER_LABELS = ["a", "b"]

# Manual cluster colors for summary profile plots and legend.
# Set to None to automatically use RColorBrewer Set1 palette.
CLUSTER_COLORS = ["#E41A1C", "#377EB8"]

# Linewidth for the average summary profile curves above the heatmap
PROFILE_LINEWIDTH = 2.5

# --- Font sizes (in pt) ---
FONT_COL_TITLE    = 16  # Sample title above each heatmap
FONT_ROW_TITLE    = 16   # Cluster title on the left
FONT_AXIS_NAME    = 13   # Heatmap x-axis labels (-10kb, center, +10kb)
FONT_ANNO_AXIS    = 13   # Average summary profile y-axis numbers
FONT_LEGEND_TITLE = 13   # Legend titles
FONT_LEGEND_TEXT  = 13   # Legend item text

# --- Within-cluster row ordering ---
# Rows within each cluster are sorted by decreasing mean signal of this sample.
SORT_BY_SAMPLE = "H3K27me2_1uM"

# --- Number of clusters to use for the heatmap ---
# Set to an integer to force a specific k (must match a k run by step_8e).
# Set to None to automatically read from {out_dir}/best_k.txt
# (written by step_8e based on the gap statistic).
N_CLUSTERS = None

# --- Layout ---
HEATMAP_WIDTH_CM = 4.0
COLUMN_GAP_CM    = 1.0
ANNO_HEIGHT_CM   = 2.5
SPLIT_WITH_GAPS  = True   # draw gap + border between cluster panels
OUT_FORMAT       = "svg"  # "svg" or "pdf"

# --- Windows to process (must match step_8d) ---
WINDOW_CONFIGS = [10_000, 5_000, 2_000]

# --- Output base (must match step_8d) ---
BASE_PREFIX  = "cbx2_retained"
OUT_DIR_BASE = os.path.join(config.HEATMAPDIR, "cbx2_retained")

# --- LSF resources ---
# Phase 3 is single-threaded R plotting. Only 2 cores are requested
# (1 active + 1 for R GC/raster threads). Memory is sized for
# use_raster=TRUE rendering of 8 heatmap panels × ~10k regions.
NUM_CORES    = 2
MAX_MEM_MB   = 16_384
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES
QUEUE        = "normal"
WALL_TIME    = "01:00"

# =============================================================================


def submit_window(window_bp: int):
    """
    Read the step_8d config JSON, merge plotting parameters, write a
    step_8f plot config, and submit the bsub job.
    """
    half_kb        = window_bp // 1_000
    prefix         = f"{BASE_PREFIX}_{half_kb}kb"
    out_dir        = os.path.join(OUT_DIR_BASE, f"{half_kb}kb")
    log_dir        = os.path.join(out_dir, "log")
    step8d_json    = os.path.join(log_dir, "step_8d_config.json")
    plot_cfg_json  = os.path.join(log_dir, "step_8f_plot_config.json")
    job_name       = f"plothm_{half_kb}kb"
    log_out        = os.path.join(log_dir, f"step_8f_{half_kb}kb.log")
    log_err        = os.path.join(log_dir, f"step_8f_{half_kb}kb.error")
    batch_file     = os.path.join(log_dir, f"step_8f_{half_kb}kb_plot.batch")

    if not os.path.exists(step8d_json):
        print(f"  [{half_kb}kb] WARNING: step_8d config not found: {step8d_json}")
        print("           Run step_8d_build_matrix.py first.")
        return

    # Load matrix-build config
    with open(step8d_json, encoding="utf-8") as fh:
        step8d_cfg = json.load(fh)

    # Resolve number of clusters
    if N_CLUSTERS is None:
        best_k_path = os.path.join(out_dir, "best_k.txt")
        if not os.path.exists(best_k_path):
            print(f"  [{half_kb}kb] WARNING: best_k.txt not found: {best_k_path}")
            print("           Run step_8e_clustering.py first (or set N_CLUSTERS explicitly).")
            return
        with open(best_k_path, encoding="utf-8") as fh:
            n_clusters = int(fh.read().strip())
        print(f"  [{half_kb}kb] best_k.txt → k={n_clusters}")
    else:
        n_clusters = int(N_CLUSTERS)

    # Validate that cluster assignment CSV exists
    asgn_csv = os.path.join(out_dir,
                            f"cluster_assignments_k{n_clusters}.csv")
    if not os.path.exists(asgn_csv):
        print(f"  [{half_kb}kb] WARNING: cluster assignments not found: {asgn_csv}")
        print("           Run step_8e_clustering.py first.")
        return

    # Pad cluster labels if k > len(CLUSTER_LABELS)
    # (happens when N_CLUSTERS=None and gap statistic picks a k > 2)
    cluster_labels_eff = list(CLUSTER_LABELS)
    while len(cluster_labels_eff) < n_clusters:
        cluster_labels_eff.append(f"Cluster {len(cluster_labels_eff) + 1}")
    if len(cluster_labels_eff) > n_clusters:
        cluster_labels_eff = cluster_labels_eff[:n_clusters]
    if cluster_labels_eff != list(CLUSTER_LABELS[:n_clusters]):
        print(f"  NOTE: CLUSTER_LABELS auto-adjusted to k={n_clusters}: "
              f"{cluster_labels_eff}")

    # Build the combined plot config (matrix paths + plot parameters)
    plot_cfg = {
        # Paths (from step_8d config)
        "norm_mats_rds":         os.path.join(out_dir,
                                              f"{prefix}_norm_mats.rds"),
        "targets_rds":           os.path.join(out_dir,
                                              f"{prefix}_targets.rds"),
        "cluster_assignments_csv": asgn_csv,
        "regions_bed":           step8d_cfg.get("regions_bed", ""),
        "window_bp":             step8d_cfg["window_bp"],
        "out_dir":               out_dir,
        "out_prefix":            prefix,
        # Plot parameters (easy to modify above)
        "sample_labels":         SAMPLE_LABELS,
        "sample_display_labels": SAMPLE_DISPLAY_LABELS,
        "sample_hues":           SAMPLE_HUES,
        "cluster_labels":        cluster_labels_eff,
        "cluster_colors":        CLUSTER_COLORS[:n_clusters] if CLUSTER_COLORS else None,
        "profile_linewidth":     PROFILE_LINEWIDTH,
        "font_col_title":        FONT_COL_TITLE,
        "font_row_title":        FONT_ROW_TITLE,
        "font_axis_name":        FONT_AXIS_NAME,
        "font_anno_axis":        FONT_ANNO_AXIS,
        "font_legend_title":     FONT_LEGEND_TITLE,
        "font_legend_text":      FONT_LEGEND_TEXT,
        "sort_by_sample":        SORT_BY_SAMPLE,
        "z_min":                 Z_MIN,
        "z_max":                 Z_MAX,
        "n_clusters":            n_clusters,
        "out_format":            OUT_FORMAT,
        "heatmap_width_cm":      HEATMAP_WIDTH_CM,
        "column_gap_cm":         COLUMN_GAP_CM,
        "anno_height_cm":        ANNO_HEIGHT_CM,
        "split_with_gaps":       SPLIT_WITH_GAPS,
    }

    with open(plot_cfg_json, "w", encoding="utf-8") as fh:
        json.dump(plot_cfg, fh, indent=2)
    print(f"  Plot config written: {plot_cfg_json}")

    # Build batch script
    r_script = os.path.join(config.CODEDIR, "utils",
                            "plot_enriched_heatmaps_clustered.R")

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

# ── Resource rationale ─────────────────────────────────────────────────────
# Plotting in R (EnrichedHeatmap + ComplexHeatmap) is single-threaded.
# 2 cores: 1 active R process + 1 for OS/GC overhead.
# 32 GB: raster rendering of 8 panels × ~10k regions in float32.
# ───────────────────────────────────────────────────────────────────────────

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

Rscript {r_script} {plot_cfg_json}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"  Submitting {job_name}  (window={half_kb}kb, k={n_clusters}, "
          f"cores={NUM_CORES})")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """Submit one Phase-3 plotting job per window configuration."""
    print("=== step_8f_plot_heatmaps.py ===")
    print(f"Window configs  : {[f'{w//1000}kb' for w in WINDOW_CONFIGS]}")
    print(f"k clusters      : {N_CLUSTERS}")
    print(f"Cluster labels  : {CLUSTER_LABELS}")
    print(f"Sort by sample  : {SORT_BY_SAMPLE}")
    print(f"Output format   : {OUT_FORMAT}")
    print()

    for window_bp in WINDOW_CONFIGS:
        submit_window(window_bp)


if __name__ == "__main__":
    main()
