#!/usr/bin/env python3

"""
step_8e_clustering.py — Phase 2: K-means clustering over genomic regions

Two modes, selected automatically:

  Mode 1 (submit, default):  python3 Scripts/step_8e_clustering.py
    Reads step_8d_config.json for each window size and submits one bsub job
    per window. Each job runs this script in cluster mode.

  Mode 2 (cluster):  invoked inside the batch job as
    python3 Scripts/step_8e_clustering.py --cluster \\
      --config_json <path/to/step_8d_config.json>

    For each k in K_CLUST (or k_min..k_max):
      - Fits KMeans (k-means++ init, n_init=K_NSTART)
      - Computes silhouette score (subsampled for speed)
      - Records within-cluster inertia (for elbow plot)
      - Computes gap statistic (Monte Carlo reference distribution)
      - Saves cluster_assignments_k{k}.csv (1-indexed)

    QC outputs (written to {out_dir}/qc/):
      elbow.png           — inertia vs k
      silhouette.png      — avg silhouette vs k
      gap_statistic.png   — gap ± 1 SE vs k with optimal-k marker
      silhouette_scores.csv
      gap_statistic.csv

    For the optimal k (max gap statistic):
      - Saves cluster_assignments_k{best_k}.csv (already done above)
      - Saves per-cluster BED files to {out_dir}/regions/
      - Writes best_k.txt

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Notes on K-means for high-dimensional genomic data
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Dimensionality
  The cluster matrix is (n_regions × 2·N_BINS) — e.g. 200 columns for a
  ±5 kb window at 100-bin resolution. This is genuinely high-dimensional.
  The curse-of-dimensionality makes Euclidean distance less discriminative
  as dimensions grow, but in practice ChIP/CUT&RUN data are correlated
  across adjacent bins, so the effective dimensionality is much lower.
  PCA/UMAP pre-reduction is an alternative but is not implemented here.

k-means++ initialisation (init='k-means++')
  STRONGLY RECOMMENDED for high-dimensional data:
  • Guarantees O(log k) approximation to the optimal solution.
  • Dramatically reduces the variance in inertia across restarts.
  • With n_init=50 restarts, sklearn tries 50 independent k-means++ seeds
    and keeps the best (lowest inertia) solution — this is the global best
    practice and what is used here.
  • Plain random init ('random') would require far more restarts to reach
    a comparable quality solution.

n_init = 50 (K_NSTART)
  Genomic data often has regions at the boundary of clusters.  50 restarts
  with k-means++ is robust for datasets up to ~50k regions.

max_iter = 500
  Default is 300; increased to 500 because high-dimensional data can
  require more Lloyd iterations to converge fully.

Silhouette score
  Measures how well each region fits its assigned cluster vs the nearest
  other cluster (range −1 to +1; higher = better separation).
  Computed on a subsample (SILHOUETTE_SUBSAMPLE_MAX) due to O(N²) cost.

Elbow method (inertia)
  Total within-cluster sum of squares vs k.  The "elbow" indicates
  diminishing returns from adding more clusters.  Useful as a sanity
  check but subjective.

Gap statistic (Tibshirani et al. 2001)
  Compares observed inertia to the inertia expected under a null reference
  distribution (uniform random data within the data bounding box).
  Optimal k = smallest k s.t. Gap(k) ≥ Gap(k+1) - SE(k+1).
  More principled than the elbow; recommended as the primary criterion.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Reads:   {out_dir}/{prefix}_cluster_matrix.tsv.gz
         {out_dir}/log/step_8d_config.json  (for BED file path)
Writes:  {out_dir}/cluster_assignments_k{k}.csv      (for all k in K_CLUST)
         {out_dir}/best_k.txt
         {out_dir}/qc/elbow.png
         {out_dir}/qc/silhouette.png
         {out_dir}/qc/gap_statistic.png
         {out_dir}/qc/silhouette_scores.csv
         {out_dir}/qc/gap_statistic.csv
         {out_dir}/regions/{prefix}_k{best_k}_cluster{c}.bed

Usage (from project root):
  python3 Scripts/step_8e_clustering.py
"""

### Import modules
import argparse
import csv
import json
import os
import sys
import time

import config
###


# =============================================================================
# ===  USER CONFIGURATION — edit ONLY this section  ===========================
# =============================================================================

# --- Samples to use for clustering ---
# Define which binned signal profiles to extract from norm_mats.rds to perform
# k-means clustering. Changing these does NOT require re-running step_8d.
CLUSTERING_SAMPLES = [
    "CBX2_DMSO",
    "CBX2_1uM",
    "H3K27me2_DMSO",
    "H3K27me2_1uM",
    "H3K27me3_DMSO"
    ]

# Number of clusters to evaluate.
#
#   Single integer  → cluster only for that k.  Example: K_CLUST = 2
#   Tuple (min, max) → sweep from min to max inclusive.  Example: K_CLUST = (2, 6)
#
# When a range is given, the gap statistic selects the optimal k automatically
# and per-cluster BED files are exported for that k only.
K_CLUST = (2)

# K-means parameters
# init='k-means++': smart seeding — strongly recommended for high-dimensional
# genomic data (see docstring above for rationale).
K_INIT       = "k-means++"  # 'k-means++' or 'random'
K_NSTART     = 100   # n_init: independent restarts; best result is kept
RANDOM_STATE = 42
MAX_ITER     = 500  # increased from default 300 for high-dimensional data

# Silhouette subsampling: max regions to use for the O(N²) distance matrix.
# Set to None to use all regions (very slow for >10k regions).
SILHOUETTE_SUBSAMPLE_MAX = 5_000

# Gap statistic: number of Monte Carlo reference datasets.
# Higher = more accurate SE estimate; each adds ~1 KMeans run.
GAP_N_REF = 20

# --- Windows to process (must match step_8d) ---
WINDOW_CONFIGS = [10_000, 5_000, 2_000]  # half-window in bp

# --- Output base (must match step_8d) ---
BASE_PREFIX  = "cbx2_retained"
OUT_DIR_BASE = os.path.join(config.HEATMAPDIR, "cbx2_retained")

# --- LSF resources ---
# Cores allow BLAS multi-threading (OMP_NUM_THREADS) during KMeans.
# Gap statistic Monte Carlo runs are parallelised with joblib.
NUM_CORES    = 8
MAX_MEM_MB   = 49_152  # 48 GB
MEM_PER_CORE = MAX_MEM_MB // NUM_CORES
QUEUE        = "normal"
WALL_TIME    = "02:00"

# =============================================================================


def _parse_k_clust(k_clust):
    """Return a sorted list of k values from K_CLUST (int or tuple)."""
    if isinstance(k_clust, int):
        return [k_clust]
    lo, hi = int(k_clust[0]), int(k_clust[1])
    if lo < 2:
        raise ValueError("K_CLUST minimum must be >= 2")
    return list(range(lo, hi + 1))


# =============================================================================
# === MODE 2: CLUSTERING (runs inside the batch job) ==========================
# =============================================================================

def run_clustering(cfg_path: str):
    """
    Load the cluster matrix, run k-means++ for each k, compute QC metrics,
    and save outputs.  Called inside the bsub batch job.
    """
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    from joblib import Parallel, delayed
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend for HPC
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    print("=== step_8e_clustering.py — cluster mode ===")

    # --- Load config and matrix ---
    with open(cfg_path, encoding="utf-8") as fh:
        cfg = json.load(fh)

    out_dir      = cfg["out_dir"]
    prefix       = cfg["prefix"]
    regions_bed  = cfg.get("regions_bed", "")

    qc_dir      = os.path.join(out_dir, "qc")
    regions_dir = os.path.join(out_dir, "regions")
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(regions_dir, exist_ok=True)

    matrix_path = os.path.join(out_dir, f"{prefix}_cluster_matrix.tsv.gz")
    rds_path    = os.path.join(out_dir, f"{prefix}_norm_mats.rds")

    if os.path.exists(rds_path):
        print(f"Extracting cluster matrix for samples: {CLUSTERING_SAMPLES}")
        extract_script = os.path.join(config.CODEDIR, "utils", "extract_cluster_matrix.R")
        samples_arg    = " ".join(CLUSTERING_SAMPLES)
        cmd = f"Rscript {extract_script} {rds_path} {matrix_path} {samples_arg}"
        if os.system(cmd) != 0:
            sys.exit("ERROR: failed to extract cluster matrix from norm_mats.rds")
    elif not os.path.exists(matrix_path):
        sys.exit(
            f"ERROR: neither {rds_path} nor {matrix_path} found.\n"
            "       Run step_8d first."
        )

    print(f"Loading cluster matrix: {matrix_path}")
    X = np.loadtxt(matrix_path, dtype=np.float32)
    n_regions, n_features = X.shape
    print(f"  Matrix shape: {X.shape}  ({n_regions} regions × {n_features} features)")

    # --- Determine k values ---
    k_values = _parse_k_clust(K_CLUST)
    print(f"\nK_CLUST = {K_CLUST}  →  evaluating k = {k_values}")
    print(f"  init     : {K_INIT}")
    print(f"  n_init   : {K_NSTART}")
    print(f"  max_iter : {MAX_ITER}")
    print(f"  gap refs : {GAP_N_REF} Monte Carlo datasets")

    # --- Silhouette subsample ---
    rng = np.random.default_rng(RANDOM_STATE)
    if SILHOUETTE_SUBSAMPLE_MAX and n_regions > SILHOUETTE_SUBSAMPLE_MAX:
        sil_idx = rng.choice(n_regions, size=SILHOUETTE_SUBSAMPLE_MAX, replace=False)
        print(f"  Silhouette subsample: {SILHOUETTE_SUBSAMPLE_MAX} / {n_regions} regions")
    else:
        sil_idx = np.arange(n_regions)

    # =========================================================================
    # Fit k-means for each k; collect inertia, labels, silhouette
    # =========================================================================
    print("\n--- Fitting k-means ---")
    inertia_list  = []
    sil_list      = []
    labels_by_k   = {}

    for k in k_values:
        km = KMeans(
            n_clusters   = k,
            init         = K_INIT,
            n_init       = K_NSTART,
            max_iter     = MAX_ITER,
            random_state = RANDOM_STATE,
            algorithm    = "lloyd",  # explicit; default may change across sklearn versions
        )
        labels = km.fit_predict(X)
        labels_by_k[k] = labels
        inertia_list.append(km.inertia_)

        # Silhouette on subsample
        X_sub  = X[sil_idx]
        lb_sub = labels[sil_idx]
        if len(np.unique(lb_sub)) < 2:
            sil = float("nan")
        else:
            sil = float(silhouette_score(X_sub, lb_sub, metric="euclidean"))
        sil_list.append(sil)

        # Cluster size report
        sizes = " | ".join(
            f"C{c+1}: {int((labels == c).sum())}"
            for c in range(k)
        )
        print(f"  k={k}  inertia={km.inertia_:.1f}  silhouette={sil:.4f}  [{sizes}]")

        # Save cluster assignment CSV (1-indexed, always)
        out_csv = os.path.join(out_dir, f"cluster_assignments_k{k}.csv")
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            for lbl in labels:
                w.writerow([int(lbl) + 1])

    # =========================================================================
    # Gap statistic (Tibshirani et al. 2001)
    # Parallelised over reference datasets using joblib.
    # =========================================================================
    print(f"\n--- Gap statistic ({GAP_N_REF} reference datasets, {NUM_CORES} jobs) ---")

    # Bounding box of the observed data (for uniform reference generation)
    col_min = X.min(axis=0)
    col_max = X.max(axis=0)

    def _ref_inertia_for_k(k, seed):
        """Fit k-means on one uniform random reference dataset."""
        rng_ref = np.random.default_rng(seed)
        X_ref = rng_ref.uniform(low=col_min, high=col_max, size=X.shape).astype(np.float32)
        km_ref = KMeans(
            n_clusters   = k,
            init         = K_INIT,
            n_init       = max(1, K_NSTART // 5),   # fewer restarts for ref runs
            max_iter     = MAX_ITER,
            random_state = seed,
            algorithm    = "lloyd",
        )
        km_ref.fit(X_ref)
        return np.log(km_ref.inertia_)

    gap_values = []
    gap_se     = []

    for k, obs_inertia in zip(k_values, inertia_list):
        log_obs = np.log(obs_inertia)
        seeds   = [RANDOM_STATE + i * 97 + k * 13 for i in range(GAP_N_REF)]
        log_refs = Parallel(n_jobs=NUM_CORES)(
            delayed(_ref_inertia_for_k)(k, s) for s in seeds
        )
        log_refs = np.array(log_refs)
        gap  = float(log_refs.mean() - log_obs)
        se   = float(log_refs.std(ddof=1) * np.sqrt(1 + 1 / GAP_N_REF))
        gap_values.append(gap)
        gap_se.append(se)
        print(f"  k={k}  gap={gap:.4f}  SE={se:.4f}")

    # Optimal k: Tibshirani 1-SE rule —
    # smallest k s.t. gap(k) >= gap(k+1) - se(k+1)
    best_k = k_values[0]
    for i in range(len(k_values) - 1):
        if gap_values[i] >= gap_values[i + 1] - gap_se[i + 1]:
            best_k = k_values[i]
            break
    else:
        best_k = k_values[int(np.argmax(gap_values))]

    print(f"\nOptimal k = {best_k}  (gap statistic, 1-SE rule)")
    with open(os.path.join(out_dir, "best_k.txt"), "w", encoding="utf-8") as fh:
        fh.write(str(best_k) + "\n")

    # =========================================================================
    # Save QC CSVs
    # =========================================================================
    sil_csv = os.path.join(qc_dir, "silhouette_scores.csv")
    with open(sil_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["k", "avg_silhouette"])
        for k, s in zip(k_values, sil_list):
            w.writerow([k, f"{s:.6f}"])

    gap_csv = os.path.join(qc_dir, "gap_statistic.csv")
    with open(gap_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["k", "gap", "se", "gap_minus_se", "is_optimal"])
        for k, g, se in zip(k_values, gap_values, gap_se):
            w.writerow([k, f"{g:.6f}", f"{se:.6f}", f"{g - se:.6f}",
                        "TRUE" if k == best_k else "FALSE"])

    # =========================================================================
    # QC plots
    # =========================================================================
    _STYLE = {
        "figure.facecolor": "white",
        "axes.spines.top":  False,
        "axes.spines.right": False,
        "font.family": "sans-serif",
    }
    plt.rcParams.update(_STYLE)

    ks = np.array(k_values)
    k_ticks = k_values  # integer x-axis

    # --- 1. Elbow plot ---
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(ks, inertia_list, "o-", color="#2c7bb6", linewidth=1.8, markersize=6)
    ax.axvline(best_k, color="#d7191c", linestyle="--", linewidth=1.2,
               label=f"optimal k={best_k}")
    ax.set_xlabel("Number of clusters (k)", fontsize=11)
    ax.set_ylabel("Within-cluster inertia", fontsize=11)
    ax.set_title("Elbow method", fontsize=12, fontweight="bold")
    ax.xaxis.set_major_locator(ticker.FixedLocator(k_ticks))
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(qc_dir, "elbow.png"), dpi=150)
    plt.close(fig)

    # --- 2. Silhouette plot ---
    valid_sil = [s if not (isinstance(s, float) and s != s) else 0.0
                 for s in sil_list]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(ks, valid_sil, "o-", color="#1a9641", linewidth=1.8, markersize=6)
    ax.axvline(best_k, color="#d7191c", linestyle="--", linewidth=1.2,
               label=f"optimal k={best_k}")
    ax.set_xlabel("Number of clusters (k)", fontsize=11)
    ax.set_ylabel("Average silhouette score", fontsize=11)
    ax.set_title("Silhouette score", fontsize=12, fontweight="bold")
    ax.xaxis.set_major_locator(ticker.FixedLocator(k_ticks))
    ax.legend(fontsize=9)
    ax.set_ylim(bottom=max(0, min(valid_sil) - 0.05))
    fig.tight_layout()
    fig.savefig(os.path.join(qc_dir, "silhouette.png"), dpi=150)
    plt.close(fig)

    # --- 3. Gap statistic plot ---
    gap_arr = np.array(gap_values)
    se_arr  = np.array(gap_se)
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.errorbar(ks, gap_arr, yerr=se_arr,
                fmt="o-", color="#7b3294", linewidth=1.8, markersize=6,
                capsize=4, elinewidth=1.2, label="Gap ± 1 SE")
    ax.axvline(best_k, color="#d7191c", linestyle="--", linewidth=1.2,
               label=f"optimal k={best_k}")
    ax.set_xlabel("Number of clusters (k)", fontsize=11)
    ax.set_ylabel("Gap statistic", fontsize=11)
    ax.set_title("Gap statistic (Tibshirani 2001)", fontsize=12, fontweight="bold")
    ax.xaxis.set_major_locator(ticker.FixedLocator(k_ticks))
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(qc_dir, "gap_statistic.png"), dpi=150)
    plt.close(fig)

    print(f"\nQC plots saved to: {qc_dir}/")

    # =========================================================================
    # Per-cluster BED files for the optimal k
    # =========================================================================
    if not regions_bed or not os.path.exists(regions_bed):
        print(f"WARNING: regions_bed not found ({regions_bed}); "
              "skipping BED export.")
    else:
        print(f"\nExporting per-cluster BED files (k={best_k}) ...")
        # Load the original BED (skip comment lines)
        bed_rows = []
        with open(regions_bed, encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith("#"):
                    bed_rows.append(line.rstrip("\n"))

        best_labels = labels_by_k[best_k]
        if len(bed_rows) != n_regions:
            print(f"  WARNING: BED has {len(bed_rows)} lines but matrix has "
                  f"{n_regions} rows — skipping BED export.")
        else:
            for c in range(best_k):
                cluster_bed = os.path.join(
                    regions_dir,
                    f"{prefix}_k{best_k}_cluster{c + 1}.bed"
                )
                cluster_rows = [bed_rows[i] for i in range(n_regions)
                                if best_labels[i] == c]
                with open(cluster_bed, "w", encoding="utf-8") as fh:
                    fh.write("\n".join(cluster_rows) + "\n")
                print(f"  Cluster {c+1}: {len(cluster_rows)} regions → "
                      f"{os.path.basename(cluster_bed)}")

    print("\n=== Phase 2 complete ===")


# =============================================================================
# === MODE 1: SUBMIT (runs on the login/interactive node) =====================
# =============================================================================

def submit_window(window_bp: int):
    """
    Locate the step_8d config JSON for this window, write a bsub script,
    and submit it.
    """
    half_kb     = window_bp // 1_000
    out_dir     = os.path.join(OUT_DIR_BASE, f"{half_kb}kb")
    log_dir     = os.path.join(out_dir, "log")
    config_json = os.path.join(log_dir, "step_8d_config.json")

    if not os.path.exists(config_json):
        print(f"  [{half_kb}kb] WARNING: step_8d config not found: {config_json}")
        print("           Run step_8d_build_matrix.py first.")
        return

    job_name   = f"cbx2_cluster_{half_kb}kb"
    log_out    = os.path.join(log_dir, f"step_8e_{half_kb}kb.log")
    log_err    = os.path.join(log_dir, f"step_8e_{half_kb}kb.error")
    batch_file = os.path.join(log_dir, f"step_8e_{half_kb}kb_clustering.batch")
    py_script  = os.path.join(config.CODEDIR, "step_8e_clustering.py")

    k_vals = _parse_k_clust(K_CLUST)
    k_summary = f"k={k_vals[0]}" if len(k_vals) == 1 else f"k={k_vals[0]}..{k_vals[-1]}"

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

# ── Resource rationale ──────────────────────────────────────────────────────
# OMP_NUM_THREADS: BLAS threads inside each KMeans run.
# NUM_CORES={NUM_CORES} also used by joblib for gap statistic parallelisation
# ({GAP_N_REF} Monte Carlo reference datasets in parallel).
# Memory: gap statistic spawns {GAP_N_REF} parallel KMeans processes, each
# holding a copy of X (~{round((10_000 * 200 * 4) / 1e6, 1)} MB for 10k×200 float32).
# ────────────────────────────────────────────────────────────────────────────

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# OMP_NUM_THREADS lets BLAS (used internally by sklearn) use multiple threads.
export OMP_NUM_THREADS={NUM_CORES}

python3 {py_script} --cluster --config_json {config_json}

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print(f"  Submitting {job_name}  (window={half_kb}kb, {k_summary}, "
          f"cores={NUM_CORES})")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="Run clustering (called automatically inside the batch job).",
    )
    parser.add_argument(
        "--config_json",
        type=str,
        default="",
        help="Path to step_8d JSON config (required for --cluster mode).",
    )
    args = parser.parse_args()

    if args.cluster:
        if not args.config_json:
            sys.exit("ERROR: --config_json is required in --cluster mode.")
        run_clustering(args.config_json)
    else:
        k_vals = _parse_k_clust(K_CLUST)
        k_summary = f"k={k_vals[0]}" if len(k_vals) == 1 else f"k={k_vals[0]}..{k_vals[-1]}"
        print("=== step_8e_clustering.py — submit mode ===")
        print(f"Window configs : {[f'{w//1000}kb' for w in WINDOW_CONFIGS]}")
        print(f"K_CLUST        : {K_CLUST}  →  {k_summary}")
        print(f"init           : {K_INIT}")
        print(f"n_init         : {K_NSTART}")
        print(f"Gap refs       : {GAP_N_REF}")
        print()
        for window_bp in WINDOW_CONFIGS:
            submit_window(window_bp)


if __name__ == "__main__":
    main()
