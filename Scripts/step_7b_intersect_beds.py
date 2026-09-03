#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_7b_intersect_beds.py — Run `bedtools intersect` between two BED files,
write the overlapping records from the -a file to a named output, and
produce publication-ready SVG figure panels summarising the intersect.

Figure panels generated
-----------------------
  Panel A  Euler / proportional Venn diagram of region counts
  Panel B  Stacked bar chart of -a region overlap status
  Panel C  Violin + strip plot: region-width distributions (overlap vs unique)
  Panel D  Per-chromosome overlap proportion bar chart

All panels are exported both individually and as a single multi-panel figure.

Equivalent bedtools call:
    bedtools intersect -a <BED_A> -b <BED_B> [options] > <output>

################################################################################
# USER TOGGLES — edit these to change the run configuration
################################################################################
"""

import argparse
import gc
import os
import subprocess
import sys
import time

# ---------------------------------------------------------------------------
# Matplotlib must use a non-interactive backend before any other mpl import
# ---------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

try:
    from matplotlib_venn import venn2
    HAS_VENN = True
except ImportError:
    HAS_VENN = False

import config

# ===========================================================================
# USER TOGGLES
# ===========================================================================
BED_B = (                           # -a file: records from THIS file are returned
    "Analysis_Data/k27me2_classification/retained_k27me2.bed"
)

BED_A = (                           # -b file: used for overlap testing only
    "Analysis_Data/k27me3_classification/lost_k27me3.bed"
)
LABEL_A = "Lost H3K27me3"       # Display label for -a dataset (used in figures & filenames)
LABEL_B = "Retained H3K27me2"           # Display label for -b dataset (used in figures & filenames)

MIN_OVERLAP = 0.05

OUT_SUBDIR = "intersections"        # Subdirectory under Analysis_Data for BED + figure output

# bedtools intersect flags (extend as needed)
BEDTOOLS_EXTRA_FLAGS = f"-f {MIN_OVERLAP} -r"
        # e.g. "-v" (non-overlapping), "-f 0.5" (min overlap)
        # Leave empty for default (any overlap, -a records returned)
# ===========================================================================

# ---------------------------------------------------------------------------
# HPC resource defaults (only used when --bsub is passed)
# ---------------------------------------------------------------------------
NUM_CORES    = 1
MAX_MEM_MB   = 8192
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES

# ---------------------------------------------------------------------------
# Publication figure style
# ---------------------------------------------------------------------------
# Colour palette (accessible, Nature-style)
COL_OVERLAP   = "#39107b"   # overlap regions  (blue)
COL_UNIQUE_A  = "#8c6bb1"   # -a only regions  (red)
COL_UNIQUE_B  = "#8c96c6"   # -b only regions  (gold)
COL_NEUTRAL   = "#Bdbdbd"   # reference bars / strips

FONT_FAMILY   = "sans-serif"
FONT_SANS     = ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]
BASE_FONTSIZE = 11           # pt — Nature / Cell figure standard
PANEL_FONTSIZE = 12          # pt — panel labels (A, B, C, D)
DPI           = 600


def _apply_pub_style():
    """Apply a clean, publication-ready rcParams style."""
    plt.rcParams.update({
        "font.family":          FONT_FAMILY,
        "font.sans-serif":      FONT_SANS,
        "font.size":            BASE_FONTSIZE,
        "axes.titlesize":       BASE_FONTSIZE,
        "axes.labelsize":       BASE_FONTSIZE,
        "xtick.labelsize":      BASE_FONTSIZE - 1,
        "ytick.labelsize":      BASE_FONTSIZE - 1,
        "legend.fontsize":      BASE_FONTSIZE - 1,
        "figure.dpi":           DPI,
        "svg.fonttype":         "none",   # keep text as text in SVG
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "axes.linewidth":       0.6,
        "xtick.major.width":    0.6,
        "ytick.major.width":    0.6,
        "xtick.major.size":     3,
        "ytick.major.size":     3,
        "lines.linewidth":      0.8,
        "patch.linewidth":      0.6,
        "savefig.bbox":         "tight",
        "savefig.pad_inches":   0.05,
    })


def _panel_label(ax, text, x=-0.18, y=1.05):
    """Add a bold upper-case panel label (A, B, C …) to an axis."""
    ax.text(x, y, text, transform=ax.transAxes,
            fontsize=PANEL_FONTSIZE, fontweight="bold",
            va="top", ha="left")


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run bedtools intersect between two BED files and generate "
            "publication-ready SVG figure panels."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bed_a",      default=BED_A,   metavar="PATH",
                        help="-a BED file: records from THIS file are written to the output.")
    parser.add_argument("--bed_b",      default=BED_B,   metavar="PATH",
                        help="-b BED file: used for overlap testing only.")
    parser.add_argument("--label_a",    default=LABEL_A, metavar="STR",
                        help="Display label for the -a dataset.")
    parser.add_argument("--label_b",    default=LABEL_B, metavar="STR",
                        help="Display label for the -b dataset.")
    parser.add_argument("--out_subdir", default=OUT_SUBDIR, metavar="STR",
                        help="Sub-directory under Analysis_Data for all outputs.")
    parser.add_argument("--extra_flags", default=BEDTOOLS_EXTRA_FLAGS, metavar="FLAGS",
                        help="Extra bedtools intersect flags (e.g. '-f 0.5 -r').")
    parser.add_argument("--bsub", action="store_true", default=False,
                        help="Submit as a bsub job (default behavior unless --direct is passed).")
    parser.add_argument("--direct", action="store_true", default=False,
                        help="Run directly without submitting a bsub job (used internally by bsub scripts).")
    parser.add_argument("--out", default=None, metavar="PATH",
                        help="Explicit output BED path (auto-derived if omitted).")
    parser.add_argument("--no_figures", action="store_true", default=False,
                        help="Skip figure generation (output BED only).")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _safe_label(label: str) -> str:
    """Convert a display label to a filesystem-safe lower-case string."""
    return label.replace(" ", "_").replace("/", "-").lower()


def build_output_path(args):
    out_dir = os.path.join(config.ANALYSIS_DATA, args.out_subdir)
    os.makedirs(out_dir, exist_ok=True)
    if args.out:
        return args.out
    fn = f"{_safe_label(args.label_a)}_in_{_safe_label(args.label_b)}.bed"
    return os.path.join(out_dir, fn)


def build_figure_dir(args):
    fig_dir = os.path.join(config.ANALYSIS_DATA, args.out_subdir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    return fig_dir


# ---------------------------------------------------------------------------
# BED loading
# ---------------------------------------------------------------------------

def _ncols(path):
    """Return the number of columns in the first non-comment line."""
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                return len(line.rstrip().split("\t"))
    return 3


def load_bed(path):
    """Load a BED file into a DataFrame. Handles 3–6+ column files."""
    ncols = _ncols(path)
    base  = ["chrom", "start", "end", "name", "score", "strand"]
    names = base[:min(ncols, 6)]
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     names=names, usecols=range(len(names)),
                     dtype=str)
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)
    df["width"] = df["end"] - df["start"]
    return df


def get_overlap_ids(out_path):
    """Return the set of (chrom, start, end) tuples present in the output BED."""
    df = load_bed(out_path)
    return set(zip(df["chrom"], df["start"], df["end"]))


# ---------------------------------------------------------------------------
# Figure panels
# ---------------------------------------------------------------------------

def panel_a_euler(ax, n_a_only, n_overlap, n_b_only, label_a, label_b):
    """
    Panel A — Euler / proportional Venn diagram.
    Uses matplotlib_venn if available; falls back to a simple bar-of-three.
    """
    if HAS_VENN:
        v = venn2(subsets=(n_a_only, n_b_only, n_overlap),
                  set_labels=("", ""),
                  ax=ax,
                  set_colors=(COL_UNIQUE_A, COL_UNIQUE_B),
                  alpha=0.75)
        # Override patch colours
        if v.get_patch_by_id("10"):
            v.get_patch_by_id("10").set_color(COL_UNIQUE_A)
        if v.get_patch_by_id("01"):
            v.get_patch_by_id("01").set_color(COL_UNIQUE_B)
        if v.get_patch_by_id("11"):
            v.get_patch_by_id("11").set_color(COL_OVERLAP)

        # Format subset labels with commas
        for lbl_id, n in [("10", n_a_only), ("01", n_b_only), ("11", n_overlap)]:
            lbl = v.get_label_by_id(lbl_id)
            if lbl:
                lbl.set_text(f"{n:,}")
                lbl.set_fontsize(BASE_FONTSIZE)

        # Custom legend
        patches = [
            mpatches.Patch(color=COL_UNIQUE_A, label=f"{label_a} only ({n_a_only:,})"),
            mpatches.Patch(color=COL_OVERLAP,  label=f"Overlap ({n_overlap:,})"),
            mpatches.Patch(color=COL_UNIQUE_B, label=f"{label_b} only ({n_b_only:,})"),
        ]
        ax.legend(handles=patches, loc="upper left", bbox_to_anchor=(1.05, 1.0), frameon=False,
                  fontsize=BASE_FONTSIZE - 1, handlelength=1)
        ax.set_axis_off()

    else:
        # Fallback: horizontal stacked bar
        total = n_a_only + n_overlap + n_b_only
        bars  = [n_a_only / total, n_overlap / total, n_b_only / total]
        colors = [COL_UNIQUE_A, COL_OVERLAP, COL_UNIQUE_B]
        labels = [f"{label_a} only ({n_a_only:,})",
                  f"Overlap ({n_overlap:,})",
                  f"{label_b} only ({n_b_only:,})"]
        left = 0
        for bar, col, lbl in zip(bars, colors, labels):
            ax.barh(0, bar, left=left, color=col, height=0.5, label=lbl)
            left += bar
        ax.set_xlim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("Fraction of total regions")
        ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1.0), frameon=False, fontsize=BASE_FONTSIZE - 1)

    ax.set_title("Region overlap", pad=4)
    _panel_label(ax, "A")


def panel_b_bar(ax, n_a, n_overlap, n_a_only, label_a, label_b):
    """
    Panel B — Stacked bar showing -a region composition.
    """
    pct_overlap = 100.0 * n_overlap / n_a if n_a > 0 else 0
    pct_unique  = 100.0 * n_a_only  / n_a if n_a > 0 else 0

    ax.bar(0, pct_overlap, color=COL_OVERLAP,  label=f"Overlap with {label_b} ({n_overlap:,})")
    ax.bar(0, pct_unique,  bottom=pct_overlap, color=COL_UNIQUE_A, label=f"{label_a} only ({n_a_only:,})")

    # Annotate
    ax.text(0, pct_overlap / 2, f"{pct_overlap:.1f}%",
            ha="center", va="center", fontsize=BASE_FONTSIZE - 1,
            color="white", fontweight="bold")
    ax.text(0, pct_overlap + pct_unique / 2, f"{pct_unique:.1f}%",
            ha="center", va="center", fontsize=BASE_FONTSIZE - 1,
            color="white", fontweight="bold")

    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(0, 105)
    ax.set_xticks([0])
    ax.set_xticklabels([label_a], rotation=0, ha="center")
    ax.set_ylabel("% of regions")
    ax.set_title(f"Overlap with {label_b}", pad=4)
    ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1.0), frameon=False, fontsize=BASE_FONTSIZE - 1)
    _panel_label(ax, "B")


def panel_c_width(ax, df_a, overlap_keys, label_a):
    """
    Panel C — Violin + strip plot of region width for overlapping vs unique -a regions.
    """
    in_ovlp = df_a[df_a.apply(
        lambda r: (r["chrom"], r["start"], r["end"]) in overlap_keys, axis=1
    )]["width"].values
    unique  = df_a[df_a.apply(
        lambda r: (r["chrom"], r["start"], r["end"]) not in overlap_keys, axis=1
    )]["width"].values

    groups  = {"Overlap": in_ovlp, f"{label_a}\nonly": unique}
    colors  = [COL_OVERLAP, COL_UNIQUE_A]
    pos     = [1, 2]

    for p, (_, vals), col in zip(pos, groups.items(), colors):
        if len(vals) == 0:
            continue
        parts = ax.violinplot(vals, positions=[p], widths=0.6,
                              showmedians=True, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(col)
            pc.set_alpha(0.6)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(1.0)

        # Jitter strip
        rng   = np.random.default_rng(42)
        jitter = rng.uniform(-0.15, 0.15, min(len(vals), 300))
        sample = vals[:300] if len(vals) > 300 else vals
        ax.scatter(p + jitter, sample, s=2, color=col, alpha=0.4,
                   linewidths=0, zorder=3)

    ax.set_xticks(pos)
    ax.set_xticklabels(list(groups.keys()))
    ax.set_ylabel("Region width (bp)")
    ax.set_yscale("log")
    ax.set_title("Region width distribution", pad=4)
    _panel_label(ax, "C")


def panel_d_chrom(ax, df_a, overlap_keys, label_a, label_b):
    """
    Panel D — Per-chromosome overlap proportion bar chart.
    """
    df_a = df_a.copy()
    df_a["overlaps"] = df_a.apply(
        lambda r: (r["chrom"], r["start"], r["end"]) in overlap_keys, axis=1
    )

    # Keep only canonical chromosomes (chr1–22, chrX, chrY)
    canon = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]
    df_a  = df_a[df_a["chrom"].isin(canon)]

    grp   = df_a.groupby("chrom")["overlaps"].agg(["sum", "count"])
    grp.columns = ["n_overlap", "n_total"]
    grp["pct"]  = 100.0 * grp["n_overlap"] / grp["n_total"]

    # Sort by canonical chromosome order
    grp = grp.reindex([c for c in canon if c in grp.index])
    grp = grp.dropna()

    x   = np.arange(len(grp))
    ax.bar(x, grp["pct"], color=COL_OVERLAP, width=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(grp.index,
                       rotation=45, ha="right", fontsize=BASE_FONTSIZE - 2)
    ax.set_ylim(0, 110)
    ax.set_ylabel(f"% of {label_a}\noverlapping {label_b}")
    ax.set_title("Per-chromosome overlap", pad=4)
    ax.axhline(100, color="grey", linewidth=0.5, linestyle="--", zorder=0)
    _panel_label(ax, "D")


# ---------------------------------------------------------------------------
# Figure assembly
# ---------------------------------------------------------------------------

def generate_figures(args, out_path, fig_dir):
    """Run all four panels and save individual + composite SVG figures."""
    _apply_pub_style()

    # Load data
    df_a = load_bed(args.bed_a)
    df_b = load_bed(args.bed_b)

    n_a       = len(df_a)
    n_b       = len(df_b)
    overlap_keys = get_overlap_ids(out_path)
    n_overlap = len(overlap_keys)
    n_a_only  = n_a - n_overlap

    # n_overlap (above) is the number of -a regions that overlap -b.
    # To compute n_b_only correctly, we need the number of -b regions that overlap -a.
    extra = f" {args.extra_flags.strip()}" if args.extra_flags.strip() else ""
    cmd_b = f"bedtools intersect -a {args.bed_b} -b {args.bed_a}{extra} -u | wc -l"
    try:
        n_b_overlap = int(subprocess.check_output(cmd_b, shell=True, text=True).strip())
    except Exception:
        n_b_overlap = n_overlap # Fallback

    n_b_only  = n_b - n_b_overlap

    label_a = args.label_a
    label_b = args.label_b
    stem    = f"{_safe_label(label_a)}_in_{_safe_label(label_b)}"

    # Programmatically calculate legend width to adjust figure width & spacing
    legend_strs = [
        f"{label_a} only ({n_a_only:,})",
        f"Overlap ({n_overlap:,})",
        f"{label_b} only ({n_b_only:,})",
        f"Overlap with {label_b} ({n_overlap:,})",
        f"Overlap with {label_a} ({n_overlap:,})",
    ]
    max_legend_chars = max(len(s) for s in legend_strs)
    # At 7pt font (~0.055 in/char) + handle width and padding (~0.4 in)
    legend_width_in = max(1.5, max_legend_chars * 0.055 + 0.4)

    # ---- Individual panels -------------------------------------------------
    panels = {}

    # --- Panel A ---
    fig_a, ax_a = plt.subplots(figsize=(2.6 + legend_width_in, 2.4))
    panel_a_euler(ax_a, n_a_only, n_overlap, n_b_only, label_a, label_b)
    path_a = os.path.join(fig_dir, f"{stem}_panelA_euler.svg")
    fig_a.savefig(path_a, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig_a)
    gc.collect()
    panels["A"] = path_a

    # --- Panel B ---
    fig_b, ax_b = plt.subplots(figsize=(1.8 + legend_width_in, 2.8))
    panel_b_bar(ax_b, n_a, n_overlap, n_a_only, label_a, label_b)
    path_b = os.path.join(fig_dir, f"{stem}_panelB_bar.svg")
    fig_b.savefig(path_b, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig_b)
    gc.collect()
    panels["B"] = path_b

    # --- Panel C ---
    fig_c, ax_c = plt.subplots(figsize=(2.8, 2.8))
    panel_c_width(ax_c, df_a, overlap_keys, label_a)
    fig_c.tight_layout()
    path_c = os.path.join(fig_dir, f"{stem}_panelC_width.svg")
    fig_c.savefig(path_c, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig_c)
    gc.collect()
    panels["C"] = path_c

    # --- Panel D ---
    fig_d, ax_d = plt.subplots(figsize=(4.5, 2.8))
    panel_d_chrom(ax_d, df_a, overlap_keys, label_a, label_b)
    fig_d.tight_layout()
    path_d = os.path.join(fig_dir, f"{stem}_panelD_chrom.svg")
    fig_d.savefig(path_d, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig_d)
    gc.collect()
    panels["D"] = path_d

    # ---- Composite multi-panel figure -------------------------------------
    w_col0 = 3.0
    w_col1 = 3.0
    gap_in = legend_width_in + 0.3
    right_margin_in = legend_width_in + 0.2
    total_fig_width = w_col0 + gap_in + w_col1 + right_margin_in
    wspace_val = gap_in / 3.0

    fig, axes = plt.subplots(
        2, 2,
        figsize=(total_fig_width, 6.0),
        gridspec_kw={"width_ratios": [1.1, 1.0], "hspace": 0.5, "wspace": wspace_val},
    )

    panel_a_euler(axes[0, 0], n_a_only, n_overlap, n_b_only, label_a, label_b)
    panel_b_bar(axes[0, 1],   n_a,      n_overlap, n_a_only,         label_a, label_b)
    panel_c_width(axes[1, 0], df_a,     overlap_keys,                label_a)
    panel_d_chrom(axes[1, 1], df_a,     overlap_keys,                label_a, label_b)

    path_all = os.path.join(fig_dir, f"{stem}_all_panels.svg")
    fig.savefig(path_all, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    gc.collect()
    panels["all"] = path_all

    return panels, {"n_a": n_a, "n_b": n_b,
                    "n_overlap": n_overlap, "n_a_only": n_a_only}


# ---------------------------------------------------------------------------
# Execution modes
# ---------------------------------------------------------------------------

def _count_lines(path):
    try:
        with open(path) as fh:
            return sum(1 for ln in fh if not ln.startswith("#"))
    except OSError:
        return "?"


def run_direct(cmd, out_path, args):
    print(f"\n{'='*62}")
    print(f"  step_7b_intersect_beds: bedtools intersect")
    print(f"{'='*62}")
    print(f"  -a  (query)   : {args.bed_a}")
    print(f"  -b  (filter)  : {args.bed_b}")
    print(f"  Extra flags   : {args.extra_flags or '(none)'}")
    print(f"  Output        : {out_path}")
    print()

    ret = os.system(cmd)
    if ret != 0:
        sys.exit(f"bedtools intersect failed with exit code {ret}")

    n_a   = _count_lines(args.bed_a)
    n_b   = _count_lines(args.bed_b)
    n_out = _count_lines(out_path)

    print(f"\n  Summary:")
    print(f"    Input -a regions : {n_a:,}" if isinstance(n_a, int) else f"    Input -a regions : {n_a}")
    print(f"    Input -b regions : {n_b:,}" if isinstance(n_b, int) else f"    Input -b regions : {n_b}")
    print(f"    Output regions   : {n_out:,}" if isinstance(n_out, int) else f"    Output regions   : {n_out}")
    if isinstance(n_a, int) and n_a > 0 and isinstance(n_out, int):
        print(f"    Overlap rate     : {100.0 * n_out / n_a:.1f}% of -a regions")

    if not args.no_figures:
        print("\n  Generating figure panels …")
        fig_dir = build_figure_dir(args)
        panels, stats = generate_figures(args, out_path, fig_dir)
        print(f"  Figures saved to  : {fig_dir}")
        for key, path in panels.items():
            print(f"    Panel {key:3s} → {os.path.basename(path)}")

    print(f"\n  Done → {out_path}")


def run_bsub(cmd, out_path, args):
    out_dir   = os.path.dirname(out_path)
    os.makedirs(out_dir, exist_ok=True)
    fig_dir   = build_figure_dir(args)

    job_name   = f"intersect_{_safe_label(args.label_a)}_in_{_safe_label(args.label_b)}"
    log_dir    = os.path.join(out_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    fig_flag = "--no_figures" if args.no_figures else ""
    self_cmd = (
        f"python3 {os.path.abspath(__file__)}"
        f" --bed_a {args.bed_a}"
        f" --bed_b {args.bed_b}"
        f" --label_a '{args.label_a}'"
        f" --label_b '{args.label_b}'"
        f" --out_subdir {args.out_subdir}"
        f" --out {out_path}"
        f" --direct"
        f" {fig_flag}"
        f" {'--extra_flags ' + repr(args.extra_flags) if args.extra_flags.strip() else ''}"
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

# ---- Run configuration (embedded for reproducibility) --------------------
# -a (query)  : {args.bed_a}
# -b (filter) : {args.bed_b}
# Label A     : {args.label_a}
# Label B     : {args.label_b}
# Extra flags : {args.extra_flags or '(none)'}
# Output      : {out_path}
# -------------------------------------------------------------------------

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

{self_cmd}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"\n{'='*62}")
    print(f"  step_7b_intersect_beds: submitting bsub job")
    print(f"{'='*62}")
    print(f"  Job name      : {job_name}")
    print(f"  -a (query)    : {args.bed_a}")
    print(f"  -b (filter)   : {args.bed_b}")
    print(f"  Extra flags   : {args.extra_flags or '(none)'}")
    print(f"  Output        : {out_path}")
    print(f"  Figures dir   : {fig_dir}")
    print(f"  Log           : {log_out}")
    print(f"  Batch script  : {batch_file}")
    print()

    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    print(f"\nMonitor with: bjobs -J {job_name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    for attr, label in [("bed_a", "-a"), ("bed_b", "-b")]:
        path = getattr(args, attr)
        if not os.path.exists(path):
            sys.exit(f"BED file not found ({label}): {path}")

    out_path = build_output_path(args)
    extra    = f" {args.extra_flags.strip()}" if args.extra_flags.strip() else ""
    cmd = (
        f"bedtools intersect"
        f" -a {args.bed_a}"
        f" -b {args.bed_b}"
        f"{extra}"
        f" > {out_path}"
    )

    if args.direct:
        run_direct(cmd, out_path, args)
    else:
        run_bsub(cmd, out_path, args)


if __name__ == "__main__":
    main()
