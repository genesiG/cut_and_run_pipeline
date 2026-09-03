#!/usr/bin/env python3

"""
qc_6_fingerprint.py - Run deepTools plotFingerprint for signal-to-noise assessment

Steps:
1. Load sample metadata via ldsample.load_samples().
2. Group processed BAMs by GROUP_PATTERNS (e.g. K27M_DMSO, K27MKO_EZH2i …).
3. For each group:
   a. Collect all IP (non-IgG) BAMs from that group.
   b. Collect all IgG BAMs from the same group as controls.
   c. Submit a bsub job running plotFingerprint with the grouped BAMs.
4. Output per group:
   (a) Fingerprint PNG plot
   (b) Raw tab-separated fingerprint data file

Inputs:  processed BAMs in config.PROCESSEDBAMDIR
Outputs: config.FINGERPRINT_DIR
"""

### Import modules
import os
import time
from typing import Union

import config
from utils import ldsample
###

# Working directory and log directory
work_dir = config.FINGERPRINT_DIR
log_dir  = os.path.join(work_dir, "log")

# Create directories if needed
os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir,  exist_ok=True)

NUM_CORES    = 8
MAX_MEM      = 32000
MEM_PER_CORE = MAX_MEM / NUM_CORES

BAM_SUFFIX = (
    f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
    if config.REMOVE_DUPLICATES
    else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
)

# Control label
CTRL_LABEL = "IgG" if config.EXPERIMENT in ("cutandrun", "cutandtag") else "Input"

def find_bam(sample_name: str) -> Union[str, None]:
    """Return the BAM path for a sample, or None if not found."""
    bam = os.path.join(config.PROCESSEDBAMDIR, f"{sample_name}.{BAM_SUFFIX}")
    return bam if os.path.exists(bam) else None



def run_fingerprint(group: str, chip_bams: list, ctrl_bams: list):
    """
    Creates and submits a bsub batch script to run deepTools plotFingerprint.

    Parameters
    ----------
    group     : group label (e.g. 'K27M_DMSO')
    chip_bams : list of IP BAM paths for this group
    ctrl_bams : list of IgG/Input BAM paths for this group
    """
    job_name = f"{group}.fingerprint"

    # Combine all BAMs; labels = basenames without suffix (IPs first, then controls)
    all_bams   = chip_bams + ctrl_bams
    all_labels = [os.path.basename(b).replace(f".{BAM_SUFFIX}", "")
                  for b in all_bams]

    bam_str   = " ".join(all_bams)
    label_str = " ".join(all_labels)



    out_png  = os.path.join(work_dir, f"{group}_fingerprint.png")
    out_data = os.path.join(work_dir, f"{group}_fingerprint.tab")

    log_out_file = os.path.join(log_dir, f"{group}.log")
    log_err_file = os.path.join(log_dir, f"{group}.error")
    batch_file   = os.path.join(log_dir, f"{group}.batch")

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

echo "Running plotFingerprint for group: {group}"
echo "BAMs: {bam_str}"

plotFingerprint \\
    --bamfiles {bam_str} \\
    --labels {label_str} \\
    --skipZeros \\
    --numberOfSamples 500000 \\
    --binSize 500 \\
    --numberOfProcessors {NUM_CORES} \\
    --plotFile {out_png} \\
    --outRawCounts {out_data} \\
    --plotTitle "{group} — Fingerprint plot"

echo "Fingerprint plot written to: {out_png}"
echo "Raw data written to: {out_data}"

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"

"""

    with open(batch_file, 'w', encoding='utf-8') as batch_fh:
        batch_fh.write(batch_cmd)

    print(f"Submitting job: {job_name}  ({len(chip_bams)} IP + {len(ctrl_bams)} ctrl BAMs)")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)


def main():
    """
    Groups processed BAMs by GROUP_PATTERNS and submits one plotFingerprint
    job per group.  Each job includes IP samples from that group + IgG controls
    from the same group.
    """
    if config.IS_PAIRED_END:
        ldsample.load_samples("paired_samples.txt")
        all_samples = list(ldsample.SAMPLES_CTL.keys())
    else:
        ldsample.load_samples("samples.txt")
        all_samples = list(ldsample.SAMPLES.keys())

    if not all_samples:
        print("WARNING: No samples loaded. Check samples.txt / paired_samples.txt.")
        return

    # Group samples by GROUP_PATTERNS
    group_samples: dict = {}
    for sn in all_samples:
        matched = False
        for pattern, label in config.GROUP_PATTERNS.items():
            if pattern in sn:
                group_samples.setdefault(label, []).append(sn)
                matched = True
                break
        if not matched:
            group_samples.setdefault("unknown", []).append(sn)

    # For each group: separate IP from IgG, then submit one job
    for group, samples in group_samples.items():
        ctrl_bams = []
        chip_bams = []

        for sn in samples:
            bam = find_bam(sn)
            if bam is None:
                print(f"  WARNING: BAM not found for {sn} — skipping.")
                continue
            if CTRL_LABEL.lower() in sn.lower():
                ctrl_bams.append(bam)
            else:
                chip_bams.append(bam)

        if not chip_bams:
            print(f"  WARNING: No IP BAMs found for group '{group}' — skipping.")
            continue
        if not ctrl_bams:
            print(f"  WARNING: No {CTRL_LABEL} BAMs found for group '{group}' — "
                  f"running without control.")

        run_fingerprint(group, chip_bams, ctrl_bams)


if __name__ == "__main__":
    main()
