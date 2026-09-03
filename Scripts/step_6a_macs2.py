#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_X_peak_calling.py - Perform peak calling with MACS2
Steps:
1. Load sample metadata via ldsample
2. For each sample and its control, build MACS2 batch command
3. Submit jobs to cluster using bsub
"""

### Import modules
import os
import time
import glob
import config
from utils import ldsample
###

### Configuration from config.py
# MACS2 parameters (set defaults if not in config)
MIN_LENGTH = getattr(config, 'MACS_MIN_LENGTH', None)
MAX_GAP = getattr(config, 'MACS_MAX_GAP', None)
Q_VALUE = getattr(config, 'MACS_Q_VALUE', 0.05)
FE_CUTOFF = getattr(config, 'MACS_FE_CUTOFF', 1.0)
CALL_SUMMITS = getattr(config, 'MACS_CALL_SUMMITS', False)
CALL_BROAD_PEAKS = getattr(config, 'MACS_CALL_BROAD_PEAKS', False)

# File paths
WORK_DIR = config.PEAKDIR
LOG_DIR = os.path.join(WORK_DIR, "log")
SUBDIRS = {
    "gappedPeaks": os.path.join(WORK_DIR, "gappedPeaks"),
    "XLS": os.path.join(WORK_DIR, "XLS"),
    "summits": os.path.join(WORK_DIR, "summits")
}

# Create output directories
os.makedirs(WORK_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
for name, path in SUBDIRS.items():
    os.makedirs(path, exist_ok=True)

def run_macs3_callpeak(group_name, treatment_param, control_param):
    """Submit MACS2 job for a group of replicates pooled together"""
    job_name = f"{group_name}.peak_calling"
    job_name = (
        job_name.replace(".peak_calling", ".broad.peak_calling")
        if CALL_BROAD_PEAKS else job_name)
    job_name = (
        job_name.replace(".peak_calling", f".maxgap{MAX_GAP}.peak_calling")
        if MAX_GAP else job_name)
    job_name = (
        job_name.replace(".peak_calling", f".qval{Q_VALUE}.peak_calling")
        if Q_VALUE != 0.05 else job_name)
    job_name = (
        job_name.replace(".peak_calling", f".cutoff{FE_CUTOFF}.peak_calling")
        if FE_CUTOFF != 1.0 else job_name)

    # Input parameter (omit -c if no controls)
    if config.EXPERIMENT == "atacseq" or not control_param:
        input_param = f"-t {treatment_param}"
    else:
        input_param = f"-t {treatment_param} \\\n  -c {control_param}"

    bam_param = "BAMPE" if config.IS_PAIRED_END else "BAM"
    keepdup = "all" if config.REMOVE_DUPLICATES else "1"

    additional_params = ""
    additional_params += f"--min-length {MIN_LENGTH} " if MIN_LENGTH else ""
    additional_params += f"--max-gap {MAX_GAP} " if MAX_GAP else ""
    additional_params += f"--fe-cutoff {FE_CUTOFF} " if FE_CUTOFF > 1 else ""
    additional_params += "--call-summits " if CALL_SUMMITS else ""
    additional_params += "--broad " if CALL_BROAD_PEAKS else ""

    out_name = f"{group_name}.qval{Q_VALUE}"
    out_name += f".minlen{MIN_LENGTH}" if MIN_LENGTH else ""
    out_name += f".maxgap{MAX_GAP}" if MAX_GAP else ""
    out_name += f".fe{FE_CUTOFF}" if FE_CUTOFF > 1 else ""
    out_name += ".broad" if CALL_BROAD_PEAKS else ""

    # Log files
    log_out = os.path.join(LOG_DIR, f"{job_name}.log")
    log_err = os.path.join(LOG_DIR, f"{job_name}.error")
    batch_file = os.path.join(LOG_DIR, f"{job_name}.batch")

    # Build MACS2 command
    macs_cmd = f"""macs3 callpeak \\
  {input_param} \\
  -f {bam_param} \\
  -g {config.GENOME_SIZE_FOR_MACS} \\
  --outdir {WORK_DIR} \\
  --keep-dup {keepdup} \\
  -n {out_name} \\
  --qvalue {Q_VALUE} \\
  {additional_params} \\"""

    # Build batch script
    batch_content = f"""#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 2  # number of cores to request
#BSUB -M 6000
#BSUB -R "rusage [mem=3000] span[hosts=1]" # Reserve the different CPU cores under the same node

# Load environment
source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq

# Run MACS2
{macs_cmd}

# Organize output files
mv {WORK_DIR}/{out_name}_peaks.gappedPeak {SUBDIRS['gappedPeaks']}
mv {WORK_DIR}/{out_name}_peaks.xls {SUBDIRS['XLS']}
mv {WORK_DIR}/{out_name}_summits.bed {SUBDIRS['summits']}

# Print resource usage at the end of the job
echo "=== Resource usage summary ==="
bjobs -l $LSB_JOBID | grep -A10 "Resource usage"
"""

    # Write batch file
    with open(batch_file, 'w') as f:
        f.write(batch_content)

    # Submit job
    print(f"Submitting {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    """
    Main function to run peak calling for all replicate groups
    """
    # Load sample metadata
    # sample_file = "paired_samples.txt" if config.IS_PAIRED_END else "samples.txt"
    # ldsample.load_samples(sample_file)
    # sample_ctl = ldsample.SAMPLES_CTL if hasattr(ldsample, 'SAMPLES_CTL') else {}

    # if not sample_ctl:
    #     print("ERROR: No sample-control mappings found. Check your metadata files.")
    #     return
    suffix = (
        f"qc.sort.rmdup.mapq{config.MAPQ}.final.bam"
        if getattr(config, 'REMOVE_DUPLICATES', False)
        else f"qc.sort.markdup.mapq{config.MAPQ}.final.bam"
    )

    # Use replicate groups from config for pattern matching with glob
    replicate_groups = config.REPLICATES
    target_antibodies = config.ANTIBODY_TARGETS

    for group in replicate_groups:
        print(f"\nProcessing replicate group: {group}")

        for target in target_antibodies:
            print(f"  Target antibody: {target}")
            # Use glob for pattern matching in the processed BAM directorys
            treatments = (
                glob.glob(os.path.join(config.PROCESSEDBAMDIR,
                                       f"*{group}*{target}*{suffix}"))
            )
            controls = (
                glob.glob(os.path.join(config.PROCESSEDBAMDIR,
                                       f"*IgG*{target}*{suffix}"))
                if config.EXPERIMENT == "cutandrun"
                else (
                    glob.glob(os.path.join(config.PROCESSEDBAMDIR,
                                       f"*{group}_Input*{suffix}"))
                    if config.EXPERIMENT == "chipseq" else []
                )
            )

            treatment_files = []
            control_files = []

            for sample_id in treatments:
                treatment_files.append(os.path.join(config.PROCESSEDBAMDIR,
                                                    f"{sample_id}"))
            if controls:
                [control_files.append(os.path.join(config.PROCESSEDBAMDIR,
                                                  f"{control_id}")) for control_id in controls]

            if not treatments:
                print(f"WARNING: No treatment samples matched group '{group}'")
                continue

            treatment_param = " ".join(treatment_files)
            control_param = " ".join(control_files) if control_files else ""

            group_name = f"{group}_{target}"
            run_macs3_callpeak(group_name, treatment_param, control_param)


if __name__ == "__main__":
    main()
