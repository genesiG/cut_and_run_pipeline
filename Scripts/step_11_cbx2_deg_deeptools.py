#!/usr/bin/env python3

"""
step_11_cbx2_deg_deeptools.py
Approach #2: Submit deepTools computeMatrix + plotProfile for each DEG group
vs. NS background.

Profile plot structure:
  - 3 separate profile plots (one per Group A, B, C)
  - Each plot: 2 region groups — target group BED vs Group D (NS) BED
  - All 4 BigWigs (DMSO Rep1/2, 1uM Rep1/2) on each plot
  - 6 SVGs total: 3 groups × 2 region types (promoters + gene bodies)
"""

import os
import sys

# Allow running from any directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import config

WORKDIR     = config.WORKDIR
BW_DIR      = config.BIGWIGDIR
OUT_DIR     = os.path.join(WORKDIR, "Analysis_Data", "cbx2_deg")
BED_DIR     = os.path.join(OUT_DIR, "bed_files")
DT_DIR      = os.path.join(OUT_DIR, "deeptools")
PLOTS_DIR   = os.path.join(OUT_DIR, "plots")
LOG_DIR     = os.path.join(DT_DIR, "log")

os.makedirs(DT_DIR,    exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(LOG_DIR,   exist_ok=True)

# --- BigWig files ---
BIGWIGS = [
    os.path.join(BW_DIR, "RP_031_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(BW_DIR, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw"),
    os.path.join(BW_DIR, "RP_055_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
    os.path.join(BW_DIR, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw"),
]
BW_STR    = " ".join(BIGWIGS)
SAMPLES   = "DMSO_Rep1 DMSO_Rep2 1uM_Rep1 1uM_Rep2"

# --- Group definitions ---
# Each entry: (group_id, label, color_hex)
GROUPS = [
    ("A", "CBX2-sensitive at DMSO",  "#DAA520"),
    ("B", "EZH2i-sensitive",         "#800080"),
    ("C", "CBX2-sensitive at EZH2i", "#4169E1"),
]
NS_LABEL = "NS background"
NS_COLOR  = "#A0A0A0"

NTHREADS = 6
MEM_MB   = 8192

# Averaged BigWig paths (will be created by bigwigAverage)
BW_AVG_DMSO = os.path.join(DT_DIR, "avg_DMSO.bw")
BW_AVG_1UM  = os.path.join(DT_DIR, "avg_1uM.bw")
BW_DMSO_1 = os.path.join(BW_DIR, "RP_031_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw")
BW_DMSO_2 = os.path.join(BW_DIR, "RP_039_AntiCBX2RP_DMSO_R1.spikein.rpkm.bw")
BW_1UM_1  = os.path.join(BW_DIR, "RP_055_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw")
BW_1UM_2  = os.path.join(BW_DIR, "RP_063_AntiCBX2RP_1uMcpd_R1.spikein.rpkm.bw")

# -----------------------------------------------------------------------
# Build one bsub batch script for all deepTools tasks
# -----------------------------------------------------------------------
batch_lines = [
    "#!/bin/bash",
    f"#BSUB -P {config.PROJECT_NAME}",
    "#BSUB -J cbx2_deg_deeptools",
    f"#BSUB -oo {LOG_DIR}/cbx2_deg_deeptools.log",
    f"#BSUB -eo {LOG_DIR}/cbx2_deg_deeptools.error",
    f"#BSUB -q {config.BSUB_QUEUE}",
    f"#BSUB -n {NTHREADS}",
    "#BSUB -M 49152",
    f"#BSUB -R 'rusage[mem={MEM_MB}] span[hosts=1]'",
    "",
    f"source /home/{config.USERNAME}/miniconda3/etc/profile.d/conda.sh",
    "conda activate chipseq",
    "",
    "# === Step 0: Average replicates ===",
    f"if [ ! -f {BW_AVG_DMSO} ]; then",
    (
        f"    bigwigAverage "
        f"-b {BW_DMSO_1} {BW_DMSO_2} "
        f"-o {BW_AVG_DMSO} "
        f"-p {NTHREADS}"
    ),
    "fi",
    f"if [ ! -f {BW_AVG_1UM} ]; then",
    (
        f"    bigwigAverage "
        f"-b {BW_1UM_1} {BW_1UM_2} "
        f"-o {BW_AVG_1UM} "
        f"-p {NTHREADS}"
    ),
    "fi",
    "",
]

for grp_id, grp_label, grp_color in GROUPS:
    grp_bed_prom = os.path.join(BED_DIR, f"promoters_group{grp_id}.bed")
    ns_bed_prom  = os.path.join(BED_DIR,  "promoters_groupD.bed")
    grp_bed_gene = os.path.join(BED_DIR, f"genebodies_group{grp_id}.bed")
    ns_bed_gene  = os.path.join(BED_DIR,  "genebodies_groupD.bed")

    mat_prom = os.path.join(DT_DIR, f"promoters_group{grp_id}.gz")
    mat_gene = os.path.join(DT_DIR, f"genebodies_group{grp_id}.gz")
    svg_prom = os.path.join(PLOTS_DIR, f"profile_promoters_group{grp_id}.svg")
    svg_gene = os.path.join(PLOTS_DIR, f"profile_genebodies_group{grp_id}.svg")

    # --- Promoter profile ---
    batch_lines += [
        f"# === Group {grp_id}: Promoter profile ===",
        f"if [ ! -f {mat_prom} ]; then",
        (
            f"    computeMatrix reference-point --referencePoint center "
            f"-b 2000 -a 2000 "
            f"-R {grp_bed_prom} {ns_bed_prom} "
            f"-S {BW_AVG_DMSO} {BW_AVG_1UM} "
            f"-o {mat_prom} "
            f"-p {NTHREADS} --missingDataAsZero "
            f"--samplesLabel DMSO_avg 1uM_avg"
        ),
        "fi",
        (
            f"plotProfile "
            f"-m {mat_prom} "
            f"-out {svg_prom} "
            f"--plotTitle 'CBX2 Signal at Promoters: Group {grp_id} vs NS' "
            f"--regionsLabel '{grp_label}' '{NS_LABEL}' "
            f"--colors '{grp_color}' '{NS_COLOR}' "
            f"--plotWidth 5.5 "
            f"--plotType 'se' "
            f"--legendLocation upper-right"
        ),
        "",
    ]

    # --- Gene body profile ---
    batch_lines += [
        f"# === Group {grp_id}: Gene body profile ===",
        f"if [ ! -f {mat_gene} ]; then",
        (
            f"    computeMatrix scale-regions "
            f"--regionBodyLength 8000 -b 4000 -a 4000 "
            f"-R {grp_bed_gene} {ns_bed_gene} "
            f"-S {BW_AVG_DMSO} {BW_AVG_1UM} "
            f"-o {mat_gene} "
            f"-p {NTHREADS} --missingDataAsZero "
            f"--samplesLabel DMSO_avg 1uM_avg"
        ),
        "fi",
        (
            f"plotProfile "
            f"-m {mat_gene} "
            f"-out {svg_gene} "
            f"--plotTitle 'CBX2 Signal over Gene Bodies: Group {grp_id} vs NS' "
            f"--regionsLabel '{grp_label}' '{NS_LABEL}' "
            f"--colors '{grp_color}' '{NS_COLOR}' "
            f"--plotWidth 5.5 "
            f"--plotType 'se' "
            f"--legendLocation upper-right"
        ),
        "",
    ]

batch_script = os.path.join(LOG_DIR, "cbx2_deg_deeptools.batch")
with open(batch_script, "w", encoding="utf-8") as f:
    f.write("\n".join(batch_lines) + "\n")

print(f"Batch script written: {batch_script}")

# Submit
os.system(f"bsub < {batch_script}")
