#!/usr/bin/env python3

"""
qc_9_peak_stats.py - Summarize peak calling results (MACS3 / csaw)

Steps:
1. Load sample metadata via ldsample.load_samples().
2. For each sample, locate its peak file (narrowPeak, broadPeak, or csaw BED).
3. Compute per-sample peak statistics:
   (a) Total peak count
   (b) Peak width distribution (min, median, max, IQR)
   (c) Per-chromosome peak counts
4. Generate per-antibody-target replicate overlap plots:
   - UpSet plot (≥3 replicates, upsetplot package)
   - Venn diagram (exactly 2 replicates, matplotlib-venn)
5. Export:
   (a) TSV: sample, n_peaks, median_width, pct_peaks_chr1-22
   (b) Peak width histograms (seaborn, overlaid)
   (c) UpSet / Venn plots per antibody target

No bsub required — lightweight, runs locally or as a single lightweight job.
"""

### Import modules
import os
import glob
import warnings
from collections import defaultdict

import pandas as pd
import numpy as np

import config
from utils import ldsample
###

# Working directory and log directory
work_dir = os.path.join(config.QCDIR1, "peak_stats")
log_dir  = os.path.join(work_dir, "log")

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir,  exist_ok=True)

# Control label
CTRL_LABEL = "IgG" if config.EXPERIMENT in ("cutandrun", "cutandtag") else "Input"


def find_peak_file(sample_name: str, antibody: str) -> str | None:
    """
    Locate a peak file for a given sample.

    Search order:
      1. MACS3 narrowPeak / broadPeak in config.PEAKDIR matching antibody.
      2. csaw BED files in config.PEAKDIR/csaw/ matching antibody.
    """
    for ext in ("broadPeak", "narrowPeak", "bed"):
        for pat in (f"*{antibody}*{ext}", f"*{sample_name}*{ext}"):
            hits = sorted(glob.glob(os.path.join(config.PEAKDIR, pat)))
            if hits:
                return hits[0]

    csaw_dir = os.path.join(config.PEAKDIR, "csaw")
    if os.path.isdir(csaw_dir):
        for pat in (f"*{antibody}*.bed", f"*{sample_name}*.bed"):
            hits = sorted(glob.glob(os.path.join(csaw_dir, pat)))
            if hits:
                return hits[0]
    return None


def load_bed(peak_file: str) -> pd.DataFrame | None:
    """
    Load a BED/narrowPeak/broadPeak file into a DataFrame.
    Returns None if file is missing or empty.
    """
    if not os.path.exists(peak_file):
        return None
    try:
        df = pd.read_csv(peak_file, sep="\t", header=None, comment="#",
                         usecols=[0, 1, 2])
        df.columns = ["chr", "start", "end"]
        df["width"] = df["end"] - df["start"]
        return df[df["width"] > 0]
    except Exception as e:
        warnings.warn(f"Could not parse {peak_file}: {e}")
        return None


def compute_peak_stats(sample_name: str, peaks: pd.DataFrame) -> dict:
    """Compute summary statistics for a single sample's peaks."""
    widths = peaks["width"].values
    chr_counts = peaks.groupby("chr")["chr"].count()

    total = len(peaks)
    # Percentage of peaks on canonical chromosomes chr1-22
    autosomes = [f"chr{i}" for i in range(1, 23)]
    n_autosome = chr_counts.reindex(autosomes, fill_value=0).sum()
    pct_chr1_22 = round(100 * n_autosome / total, 2) if total > 0 else 0

    return {
        "sample":          sample_name,
        "n_peaks":         total,
        "min_width":       int(np.min(widths)),
        "median_width":    float(np.median(widths)),
        "max_width":       int(np.max(widths)),
        "IQR_width":       float(np.percentile(widths, 75) - np.percentile(widths, 25)),
        "pct_peaks_chr1_22": pct_chr1_22,
    }


def plot_width_histograms(sample_peaks: dict):
    """
    Overlaid seaborn histogram of peak widths for all samples.
    One PDF per antibody target for readability.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Group by antibody
    ab_groups = defaultdict(dict)
    for sn, (peaks, ab) in sample_peaks.items():
        ab_groups[ab][sn] = peaks

    for ab, sp_dict in ab_groups.items():
        if not sp_dict:
            continue
        fig, ax = plt.subplots(figsize=(8, 5))
        for sn, peaks in sp_dict.items():
            widths = peaks["width"].values
            widths = widths[widths < 10000]  # cap outliers for visualisation
            sns.kdeplot(widths, ax=ax, label=sn, alpha=0.7, linewidth=1.2)

        ax.set_xlabel("Peak width (bp)")
        ax.set_ylabel("Density")
        ax.set_title(f"{config.PROJECT_NAME} — Peak Width Distribution: {ab}")
        ax.legend(fontsize=6, ncol=2)
        plt.tight_layout()
        out_path = os.path.join(work_dir, f"peak_width_hist_{ab}.pdf")
        fig.savefig(out_path, bbox_inches="tight")
        plt.close()
        print(f"  Width histogram written: {out_path}")


def run_upset_or_venn(target: str, sample_peak_files: dict):
    """
    Compute pairwise replicate peak overlaps via bedtools intersect.
    - 2 samples   → Venn diagram (matplotlib-venn)
    - ≥3 samples  → UpSet plot (upsetplot)

    sample_peak_files : {sample_name: peak_file_path}
    """
    import subprocess
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    samples = list(sample_peak_files.keys())
    if len(samples) < 2:
        print(f"  Skipping overlap for {target}: only {len(samples)} sample(s).")
        return

    # Build set memberships via bedtools
    membership = {}
    for sn, pf in sample_peak_files.items():
        try:
            tmp_out = os.path.join(log_dir, f"{sn}_peaks_sorted.bed")
            subprocess.run(
                f"sort -k1,1 -k2,2n {pf} > {tmp_out}",
                shell=True, check=True
            )
            membership[sn] = tmp_out
        except subprocess.CalledProcessError as e:
            warnings.warn(f"bedtools sort failed for {sn}: {e}")

    if len(samples) == 2:
        try:
            from matplotlib_venn import venn2
        except ImportError:
            warnings.warn("matplotlib-venn not installed; skipping Venn diagram.")
            return

        s1, s2 = samples
        f1, f2 = membership[s1], membership[s2]

        # Count unique and overlapping peaks
        def count_lines(bedfile):
            with open(bedfile) as fh:
                return sum(1 for l in fh if l.strip() and not l.startswith("#"))

        try:
            overlap_cmd = (f"bedtools intersect -a {f1} -b {f2} -u | wc -l")
            n_overlap = int(subprocess.check_output(overlap_cmd, shell=True).decode().strip())
            n1 = count_lines(f1)
            n2 = count_lines(f2)
            n1_only = n1 - n_overlap
            n2_only = n2 - n_overlap

            fig, ax = plt.subplots(figsize=(5, 5))
            venn2(subsets=(n1_only, n2_only, n_overlap),
                  set_labels=(s1, s2), ax=ax)
            ax.set_title(f"{target} — Peak Overlap")
            out_path = os.path.join(work_dir, f"peak_overlap_{target}_venn.pdf")
            fig.savefig(out_path, bbox_inches="tight")
            plt.close()
            print(f"  Venn diagram written: {out_path}")
        except Exception as e:
            warnings.warn(f"Venn diagram failed for {target}: {e}")

    else:
        try:
            from upsetplot import UpSet, from_memberships
        except ImportError:
            warnings.warn("upsetplot not installed; skipping UpSet plot.")
            return

        try:
            # Determine which peaks from s1 overlap each other sample
            # Build binary membership per peak in sample 1 (anchor)
            anchor = samples[0]
            anchor_file = membership[anchor]

            with open(anchor_file) as fh:
                anchor_peaks = [l.strip() for l in fh if l.strip()]

            # For each other sample, find overlapping peaks
            memberships_list = []
            for peak_line in anchor_peaks:
                mem = [anchor]
                for sn in samples[1:]:
                    if sn not in membership:
                        continue
                    try:
                        res = subprocess.run(
                            f"echo '{peak_line}' | bedtools intersect -a stdin -b {membership[sn]} -u | wc -l",
                            shell=True, capture_output=True, text=True
                        )
                        if res.stdout.strip() == "1":
                            mem.append(sn)
                    except Exception:
                        pass
                memberships_list.append(frozenset(mem))

            data = from_memberships(memberships_list)
            fig = plt.figure(figsize=(10, 5))
            upset = UpSet(data, subset_size="count", show_counts=True)
            upset.plot(fig)
            plt.suptitle(f"{target} — Peak Overlap (UpSet)")
            out_path = os.path.join(work_dir, f"peak_overlap_{target}_upset.pdf")
            fig.savefig(out_path, bbox_inches="tight")
            plt.close()
            print(f"  UpSet plot written: {out_path}")
        except Exception as e:
            warnings.warn(f"UpSet plot failed for {target}: {e}")


def main():
    """
    Main function: loads samples, computes peak statistics,
    exports summary TSV and plots.
    """
    import matplotlib
    matplotlib.use("Agg")

    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        all_samples = list(ldsample.SAMPLES_CTL.keys())
    else:
        ldsample.load_samples("samples.txt")
        all_samples = list(ldsample.SAMPLES.keys())

    if not all_samples:
        print("WARNING: No samples loaded.")
        return

    stats_rows      = []
    sample_peaks    = {}   # {sample_name: (DataFrame, antibody)}
    target_peak_files = defaultdict(dict)  # {antibody: {sample_name: peak_file}}

    for sn in all_samples:
        # Identify antibody
        antibody = "unknown"
        for pattern, label in config.ANTIBODY_PATTERNS.items():
            if pattern.lower() in sn.lower():
                antibody = label
                break

        # Skip control samples for peak stats (they typically lack peak files)
        if antibody.lower() in (CTRL_LABEL.lower(), "igg", "input"):
            continue

        peak_file = find_peak_file(sn, antibody)
        if peak_file is None:
            print(f"WARNING: No peak file for {sn} (antibody={antibody}) — skipping.")
            continue

        peaks = load_bed(peak_file)
        if peaks is None or len(peaks) == 0:
            print(f"WARNING: Empty or unreadable peak file for {sn} — skipping.")
            continue

        print(f"  {sn}: {len(peaks)} peaks from {os.path.basename(peak_file)}")
        stats_rows.append(compute_peak_stats(sn, peaks))
        sample_peaks[sn] = (peaks, antibody)
        target_peak_files[antibody][sn] = peak_file

    # --- Export TSV ---
    if not stats_rows:
        print("WARNING: No peak statistics computed. Check peak file locations.")
        return

    stats_df = pd.DataFrame(stats_rows)
    tsv_path = os.path.join(work_dir, "peak_stats_summary.tsv")
    stats_df.to_csv(tsv_path, sep="\t", index=False)
    print(f"\nSummary TSV written to: {tsv_path}")
    print(stats_df.to_string(index=False))

    # --- Width histograms ---
    print("\nGenerating peak width histograms...")
    plot_width_histograms(sample_peaks)

    # --- Overlap plots per target ---
    print("\nGenerating replicate overlap plots...")
    for target, spf in target_peak_files.items():
        print(f"  Target: {target} ({len(spf)} samples)")
        run_upset_or_venn(target, spf)

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
