#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_5c_callpeaks_csaw.py - Absolute peak calling with csaw (differential binding vs control)

Steps:
  1. Load SAMPLES and ANTIBODY_TARGETS from config
  2. For each (sample, target) pair, collect ALL matching target BAMs and ALL
     matching control BAMs (Input for ChIP-seq, IgG for CUT&RUN) using glob
  3. Build a bsub batch script that calls callpeaks_csaw.R passing
     all BAM paths as space-separated --chip_bams / --ctrl_bams arguments
  4. Submit one job per sample/target pair; replicates are passed in full
"""

### Import modules
import os
import time
import glob
import re
import config
###

### Configuration
BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if config.REMOVE_DUPLICATES
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

# Control antibody label depends on experiment type
CTRL_LABEL = "IgG" if config.EXPERIMENT in ("cutandrun", "cutandtag") else "Input"

# Output directories (created by the R script; define here for log path)
CSAWDIR  = getattr(config, "CSAW_OUTDIR", os.path.join(config.PEAKDIR, "csaw"))
LOG_DIR  = os.path.join(CSAWDIR, "log")

os.makedirs(LOG_DIR,  exist_ok=True)
###


def find_bams(pattern: str) -> list[str]:
    """Return sorted list of BAM paths matching pattern in PROCESSEDBAMDIR."""
    return sorted(glob.glob(os.path.join(config.PROCESSEDBAMDIR, pattern)))


def run_csaw_peakcalling(sample: str, target: str,
                         chip_bams: list[str], ctrl_bams: list[str]) -> None:
    """
    Build and submit a bsub batch script that runs callpeaks_csaw.R
    for a sample/target pair, passing ALL replicate BAMs for both target and
    control so that csaw can model them jointly.
    """
    sample_clean = re.sub(r"[_\s]+$", "", sample)
    sample_id  = f"{sample_clean}_{target}"
    job_name   = f"{sample_id}_csaw_callpeaks"

    log_out    = os.path.join(LOG_DIR, f"{job_name}.log")
    log_err    = os.path.join(LOG_DIR, f"{job_name}.error")
    batch_file = os.path.join(LOG_DIR, f"{job_name}.batch")

    r_script   = os.path.join(config.CODEDIR, "utils", "callpeaks_csaw.R")

    # Space-separated lists of BAMs passed verbatim to Rscript
    chip_bams_str = " ".join(chip_bams)
    ctrl_bams_str = " ".join(ctrl_bams)

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 1
#BSUB -M 16384
#BSUB -R "rusage[mem=16384] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq
export RETICULATE_PYTHON=$(which python3)

# Remove the virtual-address-space cap so R can allocate large contiguous
# vectors
ulimit -v unlimited 2>/dev/null || true

Rscript {r_script} \\
  --chip_bams {chip_bams_str} \\
  --ctrl_bams {ctrl_bams_str} \\
  --sample_id {sample_id} \\
  --workdir   {config.WORKDIR} \\
  --codedir   {config.CODEDIR}

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    with open(batch_file, "w", encoding="utf-8") as f:
        f.write(batch_cmd)

    print(f"Submitting {job_name}")
    print(f"  Target BAMs ({len(chip_bams)}):")
    for b in chip_bams:
        print(f"    {os.path.basename(b)}")
    print(f"  Ctrl BAMs ({len(ctrl_bams)}):")
    for b in ctrl_bams:
        print(f"    {os.path.basename(b)}")

    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    samples = config.SAMPLES
    targets = config.ANTIBODY_TARGETS

    for sample in samples:
        for target in targets:
            print(f"\n=== Processing: {sample} / {target} ===")

            # --- Find all Target BAMs for this sample/target ---
            chip_bams = find_bams(f"*{sample}*{target}*{BAM_SUFFIX}")
            if not chip_bams:
                print(f"  WARNING: No BAMs found for {sample}/{target} — skipping")
                continue
            print(f"  Found {len(chip_bams)} target BAM(s)")

            # --- Find all matched control BAMs for this sample ---
            ctrl_bams = find_bams(f"*{sample}*{CTRL_LABEL}*{BAM_SUFFIX}")
            if not ctrl_bams:
                ctrl_bams = find_bams(f"*{CTRL_LABEL}*{sample}*{BAM_SUFFIX}")
            if not ctrl_bams:
                ctrl_bams = find_bams(f"*{CTRL_LABEL}*{target}*{BAM_SUFFIX}")
            if not ctrl_bams:
                print(f"  WARNING: No {CTRL_LABEL} BAMs found for {sample} — skipping")
                continue
            print(f"  Found {len(ctrl_bams)} control BAM(s)")

            run_csaw_peakcalling(sample, target, chip_bams, ctrl_bams)


if __name__ == "__main__":
    main()
