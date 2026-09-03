#!/usr/bin/env python3

"""
step_4c_spikein_normalization.py
---------------------------------

Computes spike-in normalization scale factors following the EpiCypher
CUTANA CUT&RUN protocol:

    spike_in_fraction = spike_in_reads / total_uniquely_aligned_reads

    SF = 1 / spike_in_fraction

where:

  spike_in_reads
      Read count from  *.qc.mapq{MAPQ}.spike_in.bam
      These are MAPQ-filtered E. coli reads, split from the combined-genome
      BAM before any human-genome processing.  No deduplication is applied
      to spike-in reads — this matches the EpiCypher protocol Step 2, which
      requires only unique alignment filtering.

  total_uniquely_aligned_reads
      Read count from  *.qc.sort.[markdup|rmdup].mapq{MAPQ}.final.bam
      These are the fully processed human BAMs: MAPQ-filtered, chrM removed,
      ENCODE DAC blacklist removed, and duplicate-handled by Picard.
      Using this as the denominator matches the EpiCypher protocol Step 1
      definition: "filter out multi-mapping reads, reads assigned to ENCODE
      DAC exclusion list regions, and duplicate reads (as desired) to
      determine the total number of uniquely aligned reads."

Scale factors are then min-normalised so that the sample with the highest
spike-in fraction (smallest raw SF) anchors at SF = 1.0, and all other
samples scale proportionally upward.  A bamCoverage-compatible SF column
is also written (inverse of the normalised SF, for use with --scaleFactor).

Output: {run_name}_spikein_SF.txt  with columns:
    ID  raw_SF  SF  bamCov_SF
"""

import os
import re
import glob
import subprocess
import time

import config

# ---------------------------------------------------------------------------
# Directory setup
# ---------------------------------------------------------------------------
scaling_dir = config.SCALINGDIR
work_dir    = os.path.join(scaling_dir, "spikein")
log_dir     = os.path.join(work_dir, "log")
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

processed_dir = config.PROCESSEDBAMDIR
MAPQ          = config.MAPQ

# ---------------------------------------------------------------------------
# BAM suffix definitions
# ---------------------------------------------------------------------------
# Spike-in BAM: MAPQ-filtered E. coli reads split in step_3 Step C,
# before any human-genome filtering or deduplication.
SPIKE_SUFFIX = f"qc.mapq{MAPQ}.spike_in.bam"

# Human BAM: fully processed by step_3 Steps D + E.
# This is the correct denominator per the EpiCypher protocol.
if config.REMOVE_DUPLICATES:
    TARGET_SUFFIX = f"qc.sort.rmdup.mapq{MAPQ}.final.bam"
else:
    TARGET_SUFFIX = f"qc.sort.markdup.mapq{MAPQ}.final.bam"


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def run_cmd(cmd):
    """Run a shell command and return stdout as text."""
    return subprocess.check_output(cmd, shell=True).decode("utf-8").strip()


def count_mapped_reads(bam_file):
    """
    Return the total number of mapped reads in a BAM file using
    samtools idxstats (sums the mapped-read column across all chromosomes).
    """
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM not found: {bam_file}")

    out = run_cmd(f"samtools idxstats {bam_file}").split("\n")
    total = 0
    for line in out:
        fields = line.split("\t")
        if len(fields) == 4:
            total += int(fields[2])
    return total


def normalize_sf_file(sf_file):
    """
    Min-normalise raw SF values so the sample with the highest spike-in
    fraction (lowest raw SF) anchors at SF = 1.0.  Also writes a
    bamCov_SF column (= raw_SF / min_raw_SF) for use with bamCoverage
    --scaleFactor.
    """
    samples, raw_sfs = [], []

    with open(sf_file, "r", encoding="utf-8") as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                samples.append(parts[0])
                raw_sfs.append(float(parts[1]))

    if not raw_sfs:
        return

    min_sf       = min(raw_sfs)
    norm_sfs     = [min_sf / sf for sf in raw_sfs]
    bamcov_sfs   = [sf / min_sf for sf in raw_sfs]

    with open(sf_file, "w", encoding="utf-8") as fh:
        fh.write("ID\traw_SF\tSF\tbamCov_SF\n")
        for sample, raw, norm, bcov in zip(samples, raw_sfs, norm_sfs, bamcov_sfs):
            fh.write(f"{sample}\t{raw:.6f}\t{norm:.6f}\t{bcov:.6f}\n")

    print(f"  Normalised SF values (min raw SF = {min_sf:.6f})")


# ---------------------------------------------------------------------------
# Core normalization logic
# ---------------------------------------------------------------------------

def compute_spikein_normalization(target, run_name=None):
    """
    Compute spike-in normalization factors for all samples matching `target`.

    Parameters
    ----------
    target   : str
        Regex pattern (Python re / R-compatible '|' syntax) matched against
        BAM filenames to select samples for this run.
    run_name : str, optional
        Label for the output file.  Defaults to `target` when not provided.
        Passing the NORMALIZATION_RUNS key keeps filenames clean when
        `target` is a complex regex.
    """
    label       = run_name if run_name else target
    output_file = os.path.join(work_dir, f"{label}_spikein_SF.txt")

    with open(output_file, "w", encoding="utf-8") as fh:
        fh.write("ID\tSF\n")

    results_count = 0
    errors        = []

    # Match human final BAMs against the target regex
    all_target_bams = glob.glob(os.path.join(processed_dir, f"*{TARGET_SUFFIX}"))
    target_bams     = [
        f for f in all_target_bams
        if re.search(target, os.path.basename(f), re.IGNORECASE)
    ]

    if not target_bams:
        print(f"WARNING: No human BAMs found for target '{target}' in {processed_dir}")
        return

    print(f"\nProcessing target: {target}")
    print(f"Found {len(target_bams)} sample(s)")

    for target_bam in target_bams:
        # Derive spike-in BAM path from the human BAM path.
        # The sample base name (everything before the suffix) is shared.
        sample_base = os.path.basename(target_bam)
        sample_stem = sample_base.replace(f".{TARGET_SUFFIX}", "")
        spike_bam   = os.path.join(processed_dir, f"{sample_stem}.{SPIKE_SUFFIX}")

        print(f"\n  Sample: {sample_stem}")

        try:
            if not os.path.exists(spike_bam):
                raise FileNotFoundError(f"Spike-in BAM not found: {spike_bam}")

            spike_reads = count_mapped_reads(spike_bam)
            target_reads = count_mapped_reads(target_bam)
            total_reads = spike_reads + target_reads

            if target_reads == 0:
                raise ValueError("Human (final) BAM contains 0 mapped reads")
            if spike_reads == 0:
                raise ValueError("Spike-in BAM contains 0 mapped reads")

            spike_fraction = spike_reads / total_reads
            SF             = 1.0 / spike_fraction

            print(f"    Spike-in reads (numerator):          {spike_reads:>12,}")
            print(f"    Human reads   (denominator):         {target_reads:>12,}")
            print(f"    Total reads:                         {total_reads:>12,}")
            print(f"    Spike-in fraction:                   {spike_fraction:.6f}")
            print(f"    Raw SF (1 / spike fraction):         {SF:.6f}")

            with open(output_file, "a", encoding="utf-8") as fh:
                fh.write(f"{sample_base}\t{SF}\n")

            results_count += 1

        except Exception as exc:
            print(f"  ERROR: {exc}")
            errors.append((sample_stem, str(exc)))

    print(f"\n{'='*50}")
    print(f"Target '{target}': {results_count} sample(s) processed successfully")
    if errors:
        print(f"  Errors ({len(errors)}):")
        for sample, msg in errors:
            print(f"    {sample}: {msg}")
    print(f"Output: {output_file}")

    if results_count > 0:
        normalize_sf_file(output_file)

    print(f"{'='*50}\n")


# ---------------------------------------------------------------------------
# Batch submission
# ---------------------------------------------------------------------------

def submit_spikein_job():
    """
    Submit a single batch job that runs compute_spikein_normalization()
    sequentially for every entry in config.NORMALIZATION_RUNS.
    """
    runs = [
        (run_name, params["targets"])
        for run_name, params in config.NORMALIZATION_RUNS.items()
        if params.get("targets")
    ]

    if not runs:
        print("ERROR: config.NORMALIZATION_RUNS is empty or has no 'targets' keys")
        return

    job_name   = "spikein_normalization_all"
    log_out    = os.path.join(log_dir, f"{job_name}.log")
    log_err    = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")

    script_dir = os.path.dirname(os.path.abspath(__file__))

    python_script = f"""
import sys
sys.path.insert(0, '{script_dir}')
from step_4c_spikein_normalization import compute_spikein_normalization

runs = {runs!r}
for run_name, target in runs:
    compute_spikein_normalization(target, run_name=run_name)

print("\\nAll runs completed.")
"""

    batch_cmd = f"""#!/bin/bash
#BSUB -J "{job_name}"
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 4
#BSUB -M 8000
#BSUB -R "rusage[mem=2000] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

python3 << 'PYTHON_EOF'
{python_script}
PYTHON_EOF

echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write(batch_cmd)

    print("Submitting spike-in normalization job")
    print(f"Runs: {', '.join(name for name, _ in runs)}")
    print(f"Spike-in BAM suffix:  {SPIKE_SUFFIX}")
    print(f"Human BAM suffix:     {TARGET_SUFFIX}")

    os.system(f"bsub < {batch_file}")
    time.sleep(1)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    submit_spikein_job()
