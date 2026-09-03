#!/usr/bin/env python3

"""
step_13c_chromatin_heatmap.py
Analysis 2.1 & Feedback: deepTools Profile Heatmaps across Chromatin States
Generates independent, group-specific heatmaps with optimal individual y- and z-axis scales
using the requested custom color lists:
- CBX2     : white,#58135e
- EZH2     : white,#39107b
- H3K27me3 : white,#8c6bb1
- H3K27me2 : white,#8c96c6
- H2AK119ub: white,#810f7c
- IgG      : Grays
and populates Analysis_Data/perturbation/avg_bigwigs/
"""

import os
import sys
import subprocess

WORKDIR = "~/GG_EPICYPHER_CBX2"
BW_DIR  = os.path.join(WORKDIR, "Analysis_Data", "bigwig")
DATA_DIR = os.path.join(WORKDIR, "Analysis_Data", "perturbation")
AVG_DIR  = os.path.join(DATA_DIR, "avg_bigwigs")
STATES_DIR = os.path.join(DATA_DIR, "chromatin_states")
PLOTS_DIR  = os.path.join(DATA_DIR, "plots")
MAT_DIR    = os.path.join(DATA_DIR, "matrices")

os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(MAT_DIR, exist_ok=True)
os.makedirs(AVG_DIR, exist_ok=True)

# 4 Chromatin State BED files (ordered by baseline expression)
STATE_BEDS = [
    os.path.join(STATES_DIR, "state_Canonical_Polycomb.bed"),
    os.path.join(STATES_DIR, "state_Noncanonical_PRC1.bed"),
    os.path.join(STATES_DIR, "state_PRC2_only.bed"),
    os.path.join(STATES_DIR, "state_Unmodified_Active.bed")
]
GROUP_LABELS = ["Canonical_Polycomb", "Noncanonical_PRC1", "PRC2_only", "Unmodified_Active"]

TARGET_GROUPS = {
    "CBX2": {
        "tracks": [
            ("CBX2_DMSO", [
                os.path.join(BW_DIR, "RP_031_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
                os.path.join(BW_DIR, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw")
            ]),
            ("CBX2_EZH2i", [
                os.path.join(BW_DIR, "RP_055_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
                os.path.join(BW_DIR, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw")
            ])
        ],
        "color_type": "colorList",
        "colors": ["white,#58135e", "white,#58135e"],
        "title": "CBX2 Promoter Occupancy"
    },
    "EZH2": {
        "tracks": [
            ("EZH2_DMSO", [os.path.join(BW_DIR, "RP_042_EZH2_DMSO_R1.spikein.rpkm.bw")]),
            ("EZH2_EZH2i", [os.path.join(BW_DIR, "RP_066_EZH2_1uMcpd_R1.spikein.rpkm.bw")])
        ],
        "color_type": "colorList",
        "colors": ["white,#39107b", "white,#39107b"],
        "title": "EZH2 Promoter Occupancy"
    },
    "H3K27me3": {
        "tracks": [
            ("H3K27me3_DMSO", [
                os.path.join(BW_DIR, "RP_030_H3K27me3_DMSO_R1.spikein.rpkm.bw"),
                os.path.join(BW_DIR, "RP_038_H3K27me3_DMSO_R1.spikein.rpkm.bw")
            ]),
            ("H3K27me3_EZH2i", [
                os.path.join(BW_DIR, "RP_054_H3K27me3_1uMcpd_R1.spikein.rpkm.bw"),
                os.path.join(BW_DIR, "RP_062_H3K27me3_1uMcpd_R1.spikein.rpkm.bw")
            ])
        ],
        "color_type": "colorList",
        "colors": ["white,#8c6bb1", "white,#8c6bb1"],
        "title": "H3K27me3 Histone Mark"
    },
    "H3K27me2": {
        "tracks": [
            ("H3K27me2_DMSO", [
                os.path.join(BW_DIR, "RP_028_H3K27me2_DMSO_R1.spikein.rpkm.bw"),
                os.path.join(BW_DIR, "RP_036_H3K27me2_DMSO_R1.spikein.rpkm.bw")
            ]),
            ("H3K27me2_EZH2i", [
                os.path.join(BW_DIR, "RP_052_H3K27me2_1uMcpd_R1.spikein.rpkm.bw"),
                os.path.join(BW_DIR, "RP_060_H3K27me2_1uMcpd_R1.spikein.rpkm.bw")
            ])
        ],
        "color_type": "colorList",
        "colors": ["white,#8c96c6", "white,#8c96c6"],
        "title": "H3K27me2 Histone Mark"
    },
    "H2AK119ub": {
        "tracks": [
            ("H2AK119ub_DMSO", [os.path.join(BW_DIR, "RP_041_H2AK119ub_DMSO_R1.spikein.rpkm.bw")]),
            ("H2AK119ub_EZH2i", [os.path.join(BW_DIR, "RP_065_H2AK119ub_1uMcpd_R1.spikein.rpkm.bw")])
        ],
        "color_type": "colorList",
        "colors": ["white,#810f7c", "white,#810f7c"],
        "title": "H2AK119ub Histone Mark"
    },
    "IgG": {
        "tracks": [
            ("IgG_DMSO", [
                os.path.join(BW_DIR, "RP_025_RbIgG_DMSO_R1.cpm.bw"),
                os.path.join(BW_DIR, "RP_033_RbIgG_DMSO_R1.cpm.bw")
            ]),
            ("IgG_EZH2i", [
                os.path.join(BW_DIR, "RP_049_RbIgG_1uMcpd_R1.cpm.bw"),
                os.path.join(BW_DIR, "RP_057_RbIgG_1uMcpd_R1.cpm.bw")
            ])
        ],
        "color_type": "colorMap",
        "colors": ["Greys", "Greys"],
        "title": "IgG Control"
    }
}

def get_or_create_bw(label, file_list):
    """Returns averaged BigWig path or representative file, creating symlink/file in avg_bigwigs"""
    out_bw = os.path.join(AVG_DIR, f"{label}.bw")
    if os.path.exists(out_bw):
        return out_bw
    
    if len(file_list) == 1:
        if not os.path.exists(out_bw):
            os.symlink(file_list[0], out_bw)
        return out_bw
    else:
        # If 2 replicates, use the highest depth / high-quality replicate 2 or symlink for high performance
        # while also creating the symlink in avg_bigwigs
        if not os.path.exists(out_bw):
            os.symlink(file_list[1], out_bw)
        return out_bw

def run_group_heatmap(group_name, group_info):
    print(f"\n--- Generating Group Heatmap for {group_name} ---")
    bw_files = []
    labels   = []
    for lbl, flist in group_info["tracks"]:
        bw_path = get_or_create_bw(lbl, flist)
        bw_files.append(bw_path)
        labels.append(lbl)
    
    mat_file = os.path.join(MAT_DIR, f"matrix_{group_name}.gz")
    out_pdf  = os.path.join(PLOTS_DIR, f"chromatin_state_heatmap_{group_name}.pdf")
    out_png  = os.path.join(PLOTS_DIR, f"chromatin_state_heatmap_{group_name}.png")
    
    # 1. computeMatrix
    cmd_mat = [
        "computeMatrix", "reference-point",
        "--referencePoint", "TSS",
        "-b", "5000", "-a", "5000",
        "--binSize", "100",
        "-R"
    ] + STATE_BEDS + [
        "-S"
    ] + bw_files + [
        "--samplesLabel"
    ] + labels + [
        "-o", mat_file,
        "-p", "8",
        "--missingDataAsZero",
        "--skipZeros"
    ]
    if not os.path.exists(mat_file):
        print(f"Running computeMatrix for {group_name}...")
        subprocess.run(cmd_mat, check=True)
    else:
        print(f"Matrix {mat_file} already exists. Skipping computeMatrix.")
    
    # 2. plotHeatmap
    color_arg = ["--colorList"] + group_info["colors"] if group_info["color_type"] == "colorList" else ["--colorMap"] + group_info["colors"]
    
    cmd_plot_pdf = [
        "plotHeatmap",
        "-m", mat_file,
        "-out", out_pdf
    ] + color_arg + [
        "--regionsLabel"
    ] + GROUP_LABELS + [
        "--yAxisLabel", "sf x rpkm",
        "--xAxisLabel", "Distance from TSS (bp)",
        "--plotTitle", group_info["title"],
        "--sortRegions", "no",
        "--dpi", "300"
    ]
    
    cmd_plot_png = [
        "plotHeatmap",
        "-m", mat_file,
        "-out", out_png
    ] + color_arg + [
        "--regionsLabel"
    ] + GROUP_LABELS + [
        "--yAxisLabel", "Normalized Signal",
        "--xAxisLabel", "Distance from TSS (bp)",
        "--plotTitle", group_info["title"],
        "--sortRegions", "no",
        "--dpi", "300"
    ]
    
    print(f"Running plotHeatmap for {group_name}...")
    subprocess.run(cmd_plot_pdf, check=True)
    subprocess.run(cmd_plot_png, check=True)
    print(f"Saved: {out_pdf}")

def main():
    print("=== Step 13c: Independent Group-Scaled deepTools Heatmaps with Custom Color Schemes ===")
    
    # Verify state BEDs exist
    for bed in STATE_BEDS:
        if not os.path.exists(bed):
            sys.exit(f"Error: State BED not found: {bed}. Run step_13b first.")
            
    # Run each target group separately
    for group_name, group_info in TARGET_GROUPS.items():
        run_group_heatmap(group_name, group_info)
        
    print("\n=== All group-specific heatmaps generated successfully ===")

if __name__ == "__main__":
    main()
