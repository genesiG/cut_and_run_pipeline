#!/usr/bin/env python3

"""
step_10a_cbx2_expression_deeptools.py
Phase 3 of CBX2 Expression Analysis

Submits bsub jobs to run deepTools computeMatrix and plotProfile over
the BED files generated in Phase 1.
"""

import os
import time
import argparse
import config

NUM_CORES = 6
MAX_MEM_MB = 32768
MEM_PER_CORE = 8192

# Bigwigs to plot
BW_FILES = [
    os.path.join(config.BIGWIGDIR, "RP_031_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_055_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
    os.path.join(config.BIGWIGDIR, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw")
]
BW_LABELS = ["DMSO_Rep1", "DMSO_Rep2", "1uM_Rep1", "1uM_Rep2"]

def parse_args():
    parser = argparse.ArgumentParser(description="Run deepTools for CBX2 Expression Analysis")
    parser.add_argument("--out_dir", default="Analysis_Data/cbx2_expression", help="Output directory")
    return parser.parse_args()

def main():
    args = parse_args()
    out_dir = os.path.abspath(args.out_dir)
    bed_dir = os.path.join(out_dir, "bed_files")
    dt_dir = os.path.join(out_dir, "deeptools")
    plots_dir = os.path.join(out_dir, "plots")
    
    os.makedirs(dt_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    log_dir = os.path.join(dt_dir, "log")
    os.makedirs(log_dir, exist_ok=True)
    
    conditions = ["control", "ezh2i", "cbx2ko"]
    
    # We submit one large batch script that does everything to keep it simple
    job_name = "cbx2_deeptools"
    log_out = os.path.join(log_dir, f"{job_name}.log")
    log_err = os.path.join(log_dir, f"{job_name}.error")
    batch_file = os.path.join(log_dir, f"{job_name}.batch")
    
    batch_lines = [
        "#!/bin/bash",
        f"#BSUB -P {config.PROJECT_NAME}",
        f"#BSUB -J {job_name}",
        f"#BSUB -oo {log_out}",
        f"#BSUB -eo {log_err}",
        f"#BSUB -q {config.BSUB_QUEUE}",
        f"#BSUB -n {NUM_CORES}",
        f"#BSUB -M {MAX_MEM_MB}",
        f'#BSUB -R "rusage [mem={MEM_PER_CORE}] span[hosts=1]"',
        "",
        f"source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh",
        "conda activate chipseq",
        ""
    ]
    
    for cond in conditions:
        # Promoters BEDs
        p_high = os.path.join(bed_dir, f"promoters_{cond}_High.bed")
        p_med  = os.path.join(bed_dir, f"promoters_{cond}_Medium.bed")
        p_low  = os.path.join(bed_dir, f"promoters_{cond}_Low.bed")
        
        # Gene Body BEDs
        g_high = os.path.join(bed_dir, f"genebodies_{cond}_High.bed")
        g_med  = os.path.join(bed_dir, f"genebodies_{cond}_Medium.bed")
        g_low  = os.path.join(bed_dir, f"genebodies_{cond}_Low.bed")
        
        p_mat = os.path.join(dt_dir, f"promoters_{cond}.gz")
        g_mat = os.path.join(dt_dir, f"genebodies_{cond}.gz")
        
        c_plots_dir = os.path.join(plots_dir, cond)
        os.makedirs(c_plots_dir, exist_ok=True)
        
        p_svg = os.path.join(c_plots_dir, f"profile_promoters_{cond}.svg")
        g_svg = os.path.join(c_plots_dir, f"profile_genebodies_{cond}.svg")
        
        # computeMatrix Promoters (reference-point center)
        batch_lines.append(
            f"computeMatrix reference-point --referencePoint center -b 2000 -a 2000 "
            f"-R {p_high} {p_med} {p_low} "
            f"-S {' '.join(BW_FILES)} "
            f"-o {p_mat} "
            f"-p {NUM_CORES} "
            f"--missingDataAsZero "
            f"--samplesLabel {' '.join(BW_LABELS)}"
        )
        # plotProfile Promoters
        # Colors: High=goldenrod (#DAA520), Medium=gray (#808080), Low=purple (#800080)
        batch_lines.append(
            f"plotProfile -m {p_mat} -out {p_svg} "
            f"--plotTitle 'CBX2 Promoters vs {cond.upper()} Expression' "
            f"--regionsLabel 'High Expression' 'Medium Expression' 'Low Expression' "
            f"--colors '#DAA520' '#808080' '#800080'"
        )
        
        # computeMatrix Gene Bodies (scale-regions)
        batch_lines.append(
            f"computeMatrix scale-regions --regionBodyLength 5000 -b 2000 -a 2000 "
            f"-R {g_high} {g_med} {g_low} "
            f"-S {' '.join(BW_FILES)} "
            f"-o {g_mat} "
            f"-p {NUM_CORES} "
            f"--missingDataAsZero "
            f"--samplesLabel {' '.join(BW_LABELS)}"
        )
        # plotProfile Gene Bodies
        batch_lines.append(
            f"plotProfile -m {g_mat} -out {g_svg} "
            f"--plotTitle 'CBX2 Gene Bodies vs {cond.upper()} Expression' "
            f"--regionsLabel 'High Expression' 'Medium Expression' 'Low Expression' "
            f"--colors '#DAA520' '#808080' '#800080'"
        )
        
    with open(batch_file, "w", encoding="utf-8") as fh:
        fh.write("\n".join(batch_lines) + "\n")
        
    print(f"Submitting {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)
    
if __name__ == "__main__":
    main()
