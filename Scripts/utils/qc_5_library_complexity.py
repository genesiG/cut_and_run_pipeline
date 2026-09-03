#!/usr/bin/env python3

"""
qc_5_library_complexity.py - Estimate library complexity and duplication metrics

Steps:
1. Load sample metadata via ldsample.load_samples().
2. For each sample:
   a. Submit a bsub job to run Picard EstimateLibraryComplexity on the processed BAM.
3. After Picard jobs are done, run a separate aggregation job that:
   a. Parses Picard output to compute NRF, PBC1.
   b. Exports a TSV table and two bar charts.

Outputs (in config.ANALYSIS_DATA/qc/library_complexity/):
   - picard/        : Picard EstimateLibraryComplexity output files (one per sample)
        - library_complexity_summary.tsv
        - library_complexity_pct_dup_barplot.pdf
        - library_complexity_nrf_barplot.pdf
"""

### Import modules
import os
import time

import config
from utils import ldsample
###

# Working directory and log directory
work_dir = os.path.join(config.QCDIR1, "library_complexity")
log_dir  = os.path.join(work_dir, "log")
picard_dir = os.path.join(work_dir, "picard")

# Create directories if needed
os.makedirs(work_dir,  exist_ok=True)
os.makedirs(log_dir,   exist_ok=True)
os.makedirs(picard_dir, exist_ok=True)

NUM_CORES    = 4
MAX_MEM      = 32000
MEM_PER_CORE = MAX_MEM / NUM_CORES

# BAM suffix — mirrors step_3_process_bam.py logic
BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if config.REMOVE_DUPLICATES
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

# Markdup BAM always has the markdup suffix (for log parsing)
MARKDUP_SUFFIX = f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"

# Control label
CTRL_LABEL = "IgG" if config.EXPERIMENT == "cutandrun" else "Input"


def run_picard_complexity(sample_name):
    """
    Creates and submits a batch script to run Picard EstimateLibraryComplexity
    on the processed BAM for a single sample.

    Input:  {config.PROCESSEDBAMDIR}/{sample_name}.{BAM_SUFFIX}
    Output: {picard_dir}/{sample_name}.complexity_metrics.txt
    """

    job_name   = f"{sample_name}.picard_complexity"
    bam_path   = os.path.join(config.ORIGINALBAMDIR, f"{sample_name}.bam")
    metrics_out = os.path.join(picard_dir, f"{sample_name}.complexity_metrics.txt")

    log_out_file = os.path.join(log_dir, f"{sample_name}.picard.log")
    log_err_file = os.path.join(log_dir, f"{sample_name}.picard.error")
    batch_file   = os.path.join(log_dir, f"{sample_name}.picard.batch")

    if not os.path.exists(bam_path):
        print(f"WARNING: BAM not found for {sample_name} at {bam_path} — skipping.")
        return

    if os.path.exists(metrics_out) and os.path.getsize(metrics_out) > 0:
        print(f"Picard metrics already exist for {sample_name} — skipping.")
        return

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

export BASH_ENV=""
export JAVA_HOME=""

echo "Running Picard EstimateLibraryComplexity for: {sample_name}"

picard EstimateLibraryComplexity \\
    I={bam_path} \\
    O={metrics_out} \\
    VERBOSITY=WARNING \\
    QUIET=true

echo "Picard output written to: {metrics_out}"

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"Submitting job: {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def run_aggregation(sample_names):
    """
    Submits a bsub aggregation job that:
    1. Parses Picard EstimateLibraryComplexity output files.
    2. Parses sambamba markdup logs for duplication rates.
    3. Computes NRF (Non-Redundant Fraction) and PBC1.
    4. Exports summary TSV and bar charts (matplotlib / seaborn).

    NRF = unique read positions / total mapped reads
        ≈ (mapped_reads - duplicate_reads) / mapped_reads

    PBC1 = N1 / N_distinct, where N1 = positions covered by exactly 1 read
         (taken from Picard ESTIMATED_LIBRARY_SIZE and PERCENT_DUPLICATION)
    """

    job_name     = f"{config.PROJECT_NAME}.complexity_aggregate"
    log_out_file = os.path.join(log_dir, "aggregate.log")
    log_err_file = os.path.join(log_dir, "aggregate.error")
    batch_file   = os.path.join(log_dir, "aggregate.batch")

    # Build the sample list as a Python literal embedded in the batch script
    samples_repr = repr(sample_names)
    flagstat_dir = os.path.join(config.QCDIR1, "alignment")

    python_script = f"""\
import os, glob, re
import pandas as pd
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

picard_dir    = "{picard_dir}"
flagstat_dir  = "{flagstat_dir}"
work_dir      = "{work_dir}"
markdup_dir   = "{config.PROCESSEDBAMDIR}"
mapq          = {config.MAPQ}
sample_names  = {samples_repr}

rows = []
for sn in sample_names:
    row = {{"sample": sn, "estimated_library_size": None,
            "pct_duplicates": None, "NRF": None, "PBC1": None}}

    # --- Picard output ---
    picard_out = os.path.join(picard_dir, f"{{sn}}.complexity_metrics.txt")
    if os.path.exists(picard_out):
        with open(picard_out) as fh:
            lines = fh.readlines()
        # Find the METRICS header line
        for i, l in enumerate(lines):
            if l.startswith("ESTIMATED_LIBRARY_SIZE") or "ESTIMATED_LIBRARY_SIZE" in l:
                try:
                    vals = lines[i+1].strip().split("\\t")
                    hdr  = l.strip().split("\\t")
                    d    = dict(zip(hdr, vals))
                    row["estimated_library_size"] = int(float(d.get("ESTIMATED_LIBRARY_SIZE", 0) or 0))
                    row["pct_duplicates"]         = float(d.get("PERCENT_DUPLICATION", 0) or 0) * 100
                except Exception:
                    pass
                break
    else:
        print(f"WARNING: Picard output not found for {{sn}}")

    # --- Flagstat for NRF ---
    flagstat_path = os.path.join(flagstat_dir, f"{{sn}}.flagstat.txt")
    if os.path.exists(flagstat_path):
        with open(flagstat_path) as fh:
            content = fh.read()
        total_m  = re.search(r"^(\\d+).*in total", content, re.M)
        mapped_m = re.search(r"^(\\d+).*mapped \\(", content, re.M)
        if total_m and mapped_m:
            total  = int(total_m.group(1))
            mapped = int(mapped_m.group(1))
            dup_frac = (row["pct_duplicates"] or 0) / 100.0
            unique_reads = mapped * (1 - dup_frac)
            row["NRF"] = round(unique_reads / mapped, 4) if mapped > 0 else None
            row["PBC1"] = row["NRF"]  # approximation when single-position info unavailable
    else:
        print(f"WARNING: flagstat not found for {{sn}}")

    rows.append(row)

df = pd.DataFrame(rows)

# Reorder df so IgG/Input is last
is_ctrl = df['sample'].str.contains("IgG|Input", case=False, na=False)
df = pd.concat([df[~is_ctrl], df[is_ctrl]])

tsv_path = os.path.join(picard_dir, "library_complexity_summary.tsv")
df.to_csv(tsv_path, sep="\\t", index=False)
print(f"Summary TSV written to: {{tsv_path}}")
print(df.to_string(index=False))

# --- Bar chart: % duplicates ---
fig_dup = px.bar(df, x="sample", y="pct_duplicates", 
                 title="Library Complexity — % Duplicate Reads",
                 labels={{"sample": "", "pct_duplicates": "% Duplicate reads"}})
fig_dup.update_traces(marker_color="#7B5EA7")
fig_dup.write_html(os.path.join(picard_dir, "library_complexity_pct_dup_barplot.html"))

# --- Bar chart: NRF ---
fig_nrf = px.bar(df, x="sample", y="NRF", 
                 title="Library Complexity — NRF (Non-Redundant Fraction)",
                 labels={{"sample": "", "NRF": "NRF"}})
fig_nrf.update_traces(marker_color="#4E9BB5")
fig_nrf.add_hline(y=0.9, line_dash="dash", line_color="red", annotation_text="NRF ≥ 0.9 (ideal)")
fig_nrf.add_hline(y=0.8, line_dash="dash", line_color="orange", annotation_text="NRF ≥ 0.8 (acceptable)")
fig_nrf.update_yaxes(range=[0, 1.05])
fig_nrf.write_html(os.path.join(picard_dir, "library_complexity_nrf_barplot.html"))

print("Interactive HTML bar charts written.")
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
#BSUB -w "ended(*.picard_complexity)"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

python3 - << 'EOF'
{python_script}
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
    Loads sample metadata and submits:
      1. One Picard EstimateLibraryComplexity bsub job per sample.
      2. A single aggregation job to parse and plot results.
    """
    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        for sample_name, paired_sample in ldsample.SAMPLES_CTL.items():
            if not paired_sample:
                print(f"WARNING: No paired sample found for {sample_name}")
            run_picard_complexity(sample_name)
    else:
        ldsample.load_samples("samples.txt")
        for sample_name in ldsample.SAMPLES.keys():
            run_picard_complexity(sample_name)

    all_samples = (list(ldsample.SAMPLES_CTL.keys())
                   if config.IS_PAIRED_END
                   else list(ldsample.SAMPLES.keys()))

    run_aggregation(all_samples)


if __name__ == "__main__":
    main()
