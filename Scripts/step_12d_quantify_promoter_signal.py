#!/usr/bin/env python3

"""
step_12d_quantify_promoter_signal.py
Analysis 1.2: Quantify promoter signal (TSS +/- 2kb) for:
CBX2, H3K27me3, H3K27me2, H2AK119ub, and IgG
using deepTools multiBigwigSummary.
"""

import os
import sys
import subprocess

WORKDIR = "~/GG_EPICYPHER_CBX2"
BW_DIR  = os.path.join(WORKDIR, "Analysis_Data", "bigwig")
OUT_DIR = os.path.join(WORKDIR, "Analysis_Data", "perturbation", "signal")
BED_DIR = os.path.join(WORKDIR, "Analysis_Data", "perturbation", "bed")
LOG_DIR = os.path.join(OUT_DIR, "log")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

PROMOTER_BED = os.path.join(BED_DIR, "promoters_2kb_all.bed")

TRACKS = {
    "CBX2_DMSO_Rep1": os.path.join(BW_DIR, "RP_031_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    "CBX2_DMSO_Rep2": os.path.join(BW_DIR, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    "CBX2_EZH2i_Rep1": os.path.join(BW_DIR, "RP_055_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
    "CBX2_EZH2i_Rep2": os.path.join(BW_DIR, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
    
    "H3K27me3_DMSO_Rep1": os.path.join(BW_DIR, "RP_030_H3K27me3_DMSO_R1.spikein.rpkm.bw"),
    "H3K27me3_DMSO_Rep2": os.path.join(BW_DIR, "RP_038_H3K27me3_DMSO_R1.spikein.rpkm.bw"),
    "H3K27me3_EZH2i_Rep1": os.path.join(BW_DIR, "RP_054_H3K27me3_1uMcpd_R1.spikein.rpkm.bw"),
    "H3K27me3_EZH2i_Rep2": os.path.join(BW_DIR, "RP_062_H3K27me3_1uMcpd_R1.spikein.rpkm.bw"),
    
    "H3K27me2_DMSO_Rep1": os.path.join(BW_DIR, "RP_028_H3K27me2_DMSO_R1.spikein.rpkm.bw"),
    "H3K27me2_DMSO_Rep2": os.path.join(BW_DIR, "RP_036_H3K27me2_DMSO_R1.spikein.rpkm.bw"),
    "H3K27me2_EZH2i_Rep1": os.path.join(BW_DIR, "RP_052_H3K27me2_1uMcpd_R1.spikein.rpkm.bw"),
    "H3K27me2_EZH2i_Rep2": os.path.join(BW_DIR, "RP_060_H3K27me2_1uMcpd_R1.spikein.rpkm.bw"),
    
    "H2AK119ub_DMSO": os.path.join(BW_DIR, "RP_041_H2AK119ub_DMSO_R1.spikein.rpkm.bw"),
    "H2AK119ub_EZH2i": os.path.join(BW_DIR, "RP_065_H2AK119ub_1uMcpd_R1.spikein.rpkm.bw"),
    
    "IgG_DMSO_Rep1": os.path.join(BW_DIR, "RP_025_RbIgG_DMSO_R1.cpm.bw"),
    "IgG_DMSO_Rep2": os.path.join(BW_DIR, "RP_033_RbIgG_DMSO_R1.cpm.bw"),
    "IgG_EZH2i_Rep1": os.path.join(BW_DIR, "RP_049_RbIgG_1uMcpd_R1.cpm.bw"),
    "IgG_EZH2i_Rep2": os.path.join(BW_DIR, "RP_057_RbIgG_1uMcpd_R1.cpm.bw")
}

def main():
    print("=== Step 12d: Quantify Promoter CUT&RUN Signal ===")
    
    if not os.path.exists(PROMOTER_BED):
        sys.exit(f"Error: Promoter BED not found: {PROMOTER_BED}")
        
    bw_list = []
    labels = []
    for lbl, bw in TRACKS.items():
        if not os.path.exists(bw):
            sys.exit(f"Error: BigWig not found: {bw}")
        bw_list.append(bw)
        labels.append(lbl)
        
    out_npz = os.path.join(OUT_DIR, "promoter_signal_matrix.npz")
    out_tab = os.path.join(OUT_DIR, "promoter_signal_summary.tab")
    
    cmd = [
        "multiBigwigSummary", "BED-file",
        "-b"
    ] + bw_list + [
        "--BED", PROMOTER_BED,
        "-o", out_npz,
        "--outRawCounts", out_tab,
        "-p", "8"
    ]
    
    print(f"Running multiBigwigSummary for {len(bw_list)} tracks...")
    print(" ".join(cmd))
    
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    print(f"Successfully generated summary table: {out_tab}")
    print("=== Step 12d complete ===")

if __name__ == "__main__":
    main()
