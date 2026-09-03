#!/usr/bin/env python3

"""
qc_8_frip.py - Compute FRiP (Fraction of Reads in Peaks) for each sample

Steps:
1. Load sample metadata via ldsample.load_samples().
2. For each sample:
   a. Locate the processed BAM in config.PROCESSEDBAMDIR.
   b. Locate the corresponding peak file in config.PEAKDIR (MACS3 or csaw).
   c. Submit a bsub job that:
      - Counts total reads via samtools view -c.
      - Counts reads overlapping peaks via bedtools intersect -u.
      - Appends a line to a shared results file.
3. After all per-sample jobs finish (handled externally), an aggregation job
   reads the results file, exports the TSV and bar chart.

Inputs:
   - Processed BAMs: config.PROCESSEDBAMDIR/{sample}.{BAM_SUFFIX}
   - Peak files: config.PEAKDIR/*.narrowPeak | *.broadPeak | *.bed
     Alternatively, csaw BED files: config.PEAKDIR/csaw/*.bed

Outputs (in config.ANALYSIS_DATA/qc/frip/):
   - frip_results.tsv        : sample, total_reads, reads_in_peaks, FRiP
   - frip_barplot.pdf        : bar chart coloured by antibody target
"""

### Import modules
import os
import time
import glob

import config
from utils import ldsample
###

# Working directory and log directory
work_dir = os.path.join(config.QCDIR1, "frip")
log_dir  = os.path.join(work_dir, "log")

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir,  exist_ok=True)

NUM_CORES    = 4
MAX_MEM      = 16000
MEM_PER_CORE = MAX_MEM / NUM_CORES

# BAM suffix — mirrors step_3_process_bam.py logic
BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if config.REMOVE_DUPLICATES
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

# Shared TSV that each per-sample job appends to
RESULTS_TSV = os.path.join(work_dir, "frip_results.tsv")

# Control label
CTRL_LABEL = "IgG" if config.EXPERIMENT in ("cutandrun", "cutandtag") else "Input"


def find_peak_file(sample_name: str, antibody: str) -> str | None:
    """
    Locate a peak BED/narrowPeak/broadPeak file for a given sample.

    Search order:
      1. MACS3 narrowPeak/broadPeak files in config.PEAKDIR matching antibody.
      2. csaw BED files in config.PEAKDIR/csaw/ matching antibody.

    Returns the first match, or None if not found.
    """
    # Collect all antibody patterns that match this sample name
    matched_ab = antibody  # use the passed antibody label directly

    # 1. MACS3 peaks (broadPeak or narrowPeak)
    for ext in ("broadPeak", "narrowPeak", "bed"):
        patterns = [
            os.path.join(config.PEAKDIR, f"*{matched_ab}*{ext}"),
            os.path.join(config.PEAKDIR, f"*{sample_name}*{ext}"),
        ]
        for pat in patterns:
            hits = sorted(glob.glob(pat))
            if hits:
                return hits[0]

    # 2. csaw peaks
    csaw_dir = os.path.join(config.PEAKDIR, "csaw")
    if os.path.isdir(csaw_dir):
        for pat in (f"*{matched_ab}*.bed", f"*{sample_name}*.bed"):
            hits = sorted(glob.glob(os.path.join(csaw_dir, pat)))
            if hits:
                return hits[0]

    return None


def run_frip_sample(sample_name: str, antibody: str):
    """
    Creates and submits a batch script to compute FRiP for a single sample.

    Steps inside batch:
      1. samtools view -c  → total reads
      2. bedtools intersect -u → reads overlapping peaks
      3. Append tab-delimited line to RESULTS_TSV
    """

    bam_path  = os.path.join(config.PROCESSEDBAMDIR, f"{sample_name}.{BAM_SUFFIX}")
    peak_file = find_peak_file(sample_name, antibody)

    if not os.path.exists(bam_path):
        print(f"WARNING: BAM not found for {sample_name} at {bam_path} — skipping.")
        return

    if peak_file is None:
        print(f"WARNING: No peak file found for {sample_name} (antibody={antibody}) — skipping.")
        return

    job_name     = f"{sample_name}.frip"
    log_out_file = os.path.join(log_dir, f"{sample_name}.log")
    log_err_file = os.path.join(log_dir, f"{sample_name}.error")
    batch_file   = os.path.join(log_dir, f"{sample_name}.batch")

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

echo "Computing FRiP for: {sample_name}"
echo "  BAM : {bam_path}"
echo "  Peaks: {peak_file}"

TOTAL_READS=$(samtools view -c {bam_path})
echo "  Total reads: $TOTAL_READS"

READS_IN_PEAKS=$(bedtools intersect -u \\
    -a {bam_path} \\
    -b {peak_file} \\
    | samtools view -c)
echo "  Reads in peaks: $READS_IN_PEAKS"

FRIP=$(python3 -c "print(round($READS_IN_PEAKS / $TOTAL_READS, 6) if $TOTAL_READS > 0 else 0)")
echo "  FRiP: $FRIP"

# Append to shared results file (use flock for concurrency safety)
(
  flock -x 200
  echo -e "{sample_name}\\t{antibody}\\t$TOTAL_READS\\t$READS_IN_PEAKS\\t$FRIP" >> {RESULTS_TSV}
) 200>{RESULTS_TSV}.lock

echo "Done."

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"Submitting job: {job_name}  peak_file={os.path.basename(peak_file)}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def run_frip_aggregation(all_samples_antibodies: list):
    """
    Submits a lightweight aggregation job that reads RESULTS_TSV and
    produces a bar chart coloured by antibody target (with IgG alongside IP).

    all_samples_antibodies : list of (sample_name, antibody) tuples
    """

    job_name     = f"{config.PROJECT_NAME}.frip_aggregate"
    log_out_file = os.path.join(log_dir, "frip_aggregate.log")
    log_err_file = os.path.join(log_dir, "frip_aggregate.error")
    batch_file   = os.path.join(log_dir, "frip_aggregate.batch")

    python_plot = f"""\
import os, sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

results_tsv = "{RESULTS_TSV}"
work_dir    = "{work_dir}"

if not os.path.exists(results_tsv):
    sys.exit(f"Results file not found: {{results_tsv}}")

df = pd.read_csv(results_tsv, sep="\\t",
                 names=["sample", "antibody", "total_reads", "reads_in_peaks", "FRiP"])
df = df.drop_duplicates(subset="sample", keep="last")
df = df.sort_values(["antibody", "sample"])

# Colour by antibody target
unique_ab = df["antibody"].unique()
palette   = plt.cm.get_cmap("tab20", len(unique_ab))
ab_colors = {{ab: palette(i) for i, ab in enumerate(unique_ab)}}
colors    = [ab_colors[ab] for ab in df["antibody"]]

fig, ax = plt.subplots(figsize=(max(8, 0.4 * len(df)), 5))
bars = ax.bar(df["sample"], df["FRiP"], color=colors, edgecolor="white", linewidth=0.5)

# Legend
from matplotlib.patches import Patch
handles = [Patch(facecolor=ab_colors[ab], label=ab) for ab in unique_ab]
ax.legend(handles=handles, title="Antibody", bbox_to_anchor=(1.01, 1),
          loc="upper left", fontsize=8)

ax.set_title("{config.PROJECT_NAME} — FRiP Score per Sample")
ax.set_xlabel("")
ax.set_ylabel("FRiP (Fraction of Reads in Peaks)")
ax.set_ylim(0, min(1.0, df["FRiP"].max() * 1.2 + 0.02))
plt.xticks(rotation=45, ha="right", fontsize=7)
plt.tight_layout()
fig.savefig(os.path.join(work_dir, "frip_barplot.pdf"), bbox_inches="tight")
plt.close()
print("Bar chart written.")
print(df.to_string(index=False))
"""

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out_file}
#BSUB -eo {log_err_file}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 2
#BSUB -M 8000
#BSUB -R "rusage[mem=4000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Write header if file is new
if [ ! -f "{RESULTS_TSV}" ]; then
    echo -e "sample\\tantibody\\ttotal_reads\\treads_in_peaks\\tFRiP" > {RESULTS_TSV}
fi

python3 - << 'EOF'
{python_plot}
EOF

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"\nSubmitting aggregation job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """
    Loads sample metadata, matches samples to their antibody target
    and peak files, and submits per-sample FRiP jobs plus one aggregation job.
    """
    # Initialise results TSV header
    if not os.path.exists(RESULTS_TSV):
        with open(RESULTS_TSV, 'w', encoding='utf-8') as fh:
            fh.write("sample\tantibody\ttotal_reads\treads_in_peaks\tFRiP\n")

    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        all_samples = list(ldsample.SAMPLES_CTL.keys())
    else:
        ldsample.load_samples("samples.txt")
        all_samples = list(ldsample.SAMPLES.keys())

    if not all_samples:
        print("WARNING: No samples loaded.")
        return

    samples_antibodies = []
    for sn in all_samples:
        # Identify antibody from sample name using ANTIBODY_PATTERNS
        antibody = "unknown"
        for pattern, label in config.ANTIBODY_PATTERNS.items():
            if pattern.lower() in sn.lower():
                antibody = label
                break
        samples_antibodies.append((sn, antibody))
        run_frip_sample(sn, antibody)

    run_frip_aggregation(samples_antibodies)


if __name__ == "__main__":
    main()
