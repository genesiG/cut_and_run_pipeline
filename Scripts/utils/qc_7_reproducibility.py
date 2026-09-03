#!/usr/bin/env python3

"""
qc_7_reproducibility.py - Inter-replicate reproducibility via deepTools

Steps:
1. Load sample metadata and separate BAMs into IgG controls and IP targets.
2. For each IP antibody target (config.ANTIBODY_PATTERNS), submit one bsub job:
   a. Runs multiBamSummary bins (--binSize 5000 --distanceBetweenBins 0).
      IgG BAMs from all groups are included alongside the target BAMs.
   b. Runs plotCorrelation (Spearman, heatmap with dendrogram) → SVG.
   c. Runs plotPCA (samples coloured by group) → SVG.
   d. Exports raw count NPZ/tab file.
3. For each cell-line / sample GROUP (config.GROUP_PATTERNS), submit an
   additional bsub job with the same analysis (all IP targets + IgG for that
   group).  Outputs are prefixed with "group_<label>_".

IgG samples are NOT processed as a separate target/group; they are appended
to every IP target/group job so all plots share a common control reference.

Inputs:  processed BAMs in config.PROCESSEDBAMDIR
         sample metadata TSV (config.METADATA/sample_metadata.txt)
Outputs: config.REPRODUCIBILITY_DIR
         *.svg  (PCA biplots and correlation heatmaps)
"""

### Import modules
import os
import time

import config
from utils import ldsample
###

# Working directory and log directory
work_dir = config.REPRODUCIBILITY_DIR
log_dir  = os.path.join(work_dir, "log")

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir,  exist_ok=True)

NUM_CORES    = 16
MAX_MEM      = 64000
MEM_PER_CORE = MAX_MEM / NUM_CORES

# BAM suffix — for spike-in species use the human-only split BAM
BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if getattr(config, 'REMOVE_DUPLICATES', False)
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

# Metadata file for group annotation in PCA
META_FILE = os.path.join(config.METADATA, "sample_metadata.txt")

# Control label
CTRL_LABEL = "IgG" if config.EXPERIMENT in ("cutandrun", "cutandtag") else "Input"


def run_reproducibility(target: str, bam_paths: list, sample_labels: list,
                        prefix: str = ""):
    """
    Submits a bsub job that runs multiBamSummary + plotCorrelation + plotPCA
    for a group of BAMs.  IgG samples are always included in every job.

    Parameters
    ----------
    target        : label used for plot titles (e.g. 'K27me3', 'BT54')
    bam_paths     : list of BAM file paths (IP + IgG)
    sample_labels : list of clean sample labels (same order as bam_paths)
    prefix        : filename prefix to avoid collisions (e.g. 'group_')
    """
    target_clean = target.replace("_", "")
    safe_target  = target.replace(" ", "_").replace("/", "_")
    file_stem    = f"{prefix}{safe_target}"
    job_name     = f"{file_stem}.reproducibility"
    bam_str      = " ".join(bam_paths)
    label_str    = " ".join(sample_labels)

    # Output files — SVG for heatmap and PCA
    npz_file     = os.path.join(work_dir, f"{file_stem}_multiBamSummary.npz")
    tab_file     = os.path.join(work_dir, f"{file_stem}_readCounts.tab")
    corr_heatmap = os.path.join(work_dir, f"{file_stem}_correlation_heatmap.svg")
    corr_matrix  = os.path.join(work_dir, f"{file_stem}_correlation_matrix.tab")
    pca_plot     = os.path.join(work_dir, f"{file_stem}_pca_biplot.svg")

    log_out_file = os.path.join(log_dir, f"{file_stem}.log")
    log_err_file = os.path.join(log_dir, f"{file_stem}.error")
    batch_file   = os.path.join(log_dir, f"{file_stem}.batch")

    # Derive PCA colours from sample_metadata.txt (embedded Python snippet)
    python_colors_snippet = f"""\
import os, sys
meta_file = "{META_FILE}"
labels    = {sample_labels!r}

try:
    with open(meta_file) as fh:
        lines = fh.readlines()
    header = lines[0].strip().split("\\t")
    rows   = [dict(zip(header, l.strip().split("\\t"))) for l in lines[1:] if l.strip()]
    meta   = {{r["ID"]: r for r in rows}}

    groups = [meta.get(lb, {{}}).get("GROUP", "unknown") for lb in labels]
    unique_grps = list(dict.fromkeys(groups))

    palette = ["#7B5EA7AA", "#4E9BB5AA", "#E88A29AA",
               "#8CC97AAA", "#D65C5CAA", "#999999AA"]
    grp2col = {{g: palette[i % len(palette)] for i, g in enumerate(unique_grps)}}
    colors  = " ".join(grp2col[g] for g in groups)
    print(colors)
except Exception as e:
    sys.stderr.write(f"WARNING: Could not parse group colors: {{e}}\\n")
    print("")  # fall back to deepTools defaults
"""

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n {NUM_CORES}
#BSUB -M {MAX_MEM}
#BSUB -R "rusage[mem={MEM_PER_CORE}] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

echo "Running reproducibility QC for: {target_clean} (prefix={prefix})"
echo "Samples ({len(bam_paths)}): {label_str}"

### 1. multiBamSummary bins — 5000 bp bins for higher resolution
multiBamSummary bins \\
    --bamfiles {bam_str} \\
    --labels {label_str} \\
    --binSize 5000 \\
    --distanceBetweenBins 0 \\
    --numberOfProcessors {NUM_CORES} \\
    --outFileName {npz_file} \\
    --outRawCounts {tab_file}

### 2. plotCorrelation — Spearman heatmap with dendrogram (SVG)
plotCorrelation \\
    --corData {npz_file} \\
    --corMethod spearman \\
    --whatToPlot heatmap \\
    --skipZeros \\
    --removeOutliers \\
    --plotTitle "{target_clean} — Spearman Correlation (IP + IgG)" \\
    --colorMap RdYlBu \\
    --plotNumbers \\
    --outFileCorMatrix {corr_matrix} \\
    --plotFile {corr_heatmap}

### 3. plotPCA — derive group colors from metadata (SVG)
PCA_COLORS=$(python3 - << 'PYEOF'
{python_colors_snippet}
PYEOF
)

if [ -n "$PCA_COLORS" ]; then
    COLOR_ARG="--colors $PCA_COLORS"
else
    COLOR_ARG=""
fi

plotPCA \\
    --corData {npz_file} \\
    --labels {label_str} \\
    --plotTitle "{target_clean} — PCA (5000 bp bins, IP + IgG)" \\
    --plotFile {pca_plot} \\
    --plotWidth 7 \\
    --plotHeight 8.5 \\
    $COLOR_ARG

echo "Outputs written to: {work_dir}"

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"Submitting job: {job_name}  ({len(bam_paths)} BAMs: "
          f"{len(bam_paths) - sum(1 for lb in sample_labels if CTRL_LABEL in lb)} IP "
          f"+ {sum(1 for lb in sample_labels if CTRL_LABEL in lb)} {CTRL_LABEL})")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """
    Groups processed IP BAMs by antibody target (ANTIBODY_PATTERNS) and by
    sample group (GROUP_PATTERNS), then submits one reproducibility job per
    grouping.  IgG BAMs (all groups combined) are appended to every job.
    """
    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        all_samples = list(ldsample.SAMPLES_CTL.keys())
    else:
        ldsample.load_samples("samples.txt")
        all_samples = list(ldsample.SAMPLES.keys())

    if not all_samples:
        print("WARNING: No samples loaded.")
        return

    # Resolve BAM paths; skip missing files
    def resolve_bam(sn):
        bam = os.path.join(config.PROCESSEDBAMDIR, f"{sn}.{BAM_SUFFIX}")
        if os.path.exists(bam):
            return bam
        print(f"WARNING: BAM not found for {sn} — skipping.")
        return None

    # Separate IgG from IP samples
    igg_bams:   list = []
    igg_labels: list = []

    # Antibody-target grouping (existing behaviour)
    antibody_samples: dict = {}
    # Cell-line / group grouping (new)
    group_samples: dict = {}

    for sn in all_samples:
        bam = resolve_bam(sn)
        if bam is None:
            continue

        is_ctrl = CTRL_LABEL.lower() in sn.lower()

        if is_ctrl:
            igg_bams.append(bam)
            igg_labels.append(sn)

        # --- Antibody grouping ---
        matched_ab = False
        for pattern, label in config.ANTIBODY_PATTERNS.items():
            if pattern.lower() in sn.lower():
                antibody_samples.setdefault(label, []).append((sn, bam))
                matched_ab = True
                break
        if not matched_ab and not is_ctrl:
            antibody_samples.setdefault("unknown", []).append((sn, bam))

        # --- Cell-line / group grouping ---
        matched_grp = False
        for pattern, label in config.GROUP_PATTERNS.items():
            if pattern.lower() in sn.lower():
                group_samples.setdefault(label, []).append((sn, bam))
                matched_grp = True
                break
        if not matched_grp:
            group_samples.setdefault("other", []).append((sn, bam))

    if not igg_bams:
        print(f"WARNING: No {CTRL_LABEL} BAMs found — jobs will run without control.")

    # ── 1. Jobs grouped by ANTIBODY TARGET ──────────────────────────────────
    print("\n=== Submitting jobs grouped by ANTIBODY TARGET ===")
    for target, sample_bam_pairs in antibody_samples.items():
        if target in config.ANTIBODY_PATTERNS.values() and \
                CTRL_LABEL.lower() in target.lower():
            # Skip IgG as a standalone antibody target group
            continue

        ip_labels  = [sb[0] for sb in sample_bam_pairs]
        ip_bams    = [sb[1] for sb in sample_bam_pairs]

        all_bams   = ip_bams   + igg_bams
        all_labels = ip_labels + igg_labels

        if len(all_bams) < 2:
            print(f"WARNING: Fewer than 2 BAMs for antibody '{target}' — skipping.")
            continue

        run_reproducibility(target, all_bams, all_labels, prefix="")

    # ── 2. Jobs grouped by SAMPLE GROUP (GROUP_PATTERNS) ────────────────────
    print("\n=== Submitting jobs grouped by SAMPLE GROUP ===")
    for grp_label, sample_bam_pairs in group_samples.items():
        # For the group view, include ALL samples (IP + IgG) from that group
        # instead of appending all IgGs from every group.
        grp_bams   = [sb[1] for sb in sample_bam_pairs]
        grp_labels = [sb[0] for sb in sample_bam_pairs]

        if len(grp_bams) < 2:
            print(f"WARNING: Fewer than 2 BAMs for group '{grp_label}' — skipping.")
            continue

        run_reproducibility(grp_label, grp_bams, grp_labels, prefix="group_")


if __name__ == "__main__":
    main()
