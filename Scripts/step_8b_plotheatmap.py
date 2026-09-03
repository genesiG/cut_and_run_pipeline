#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_6b_plotheatmap.py - Generate heatmaps (without clustering)
"""

import os
import time
import argparse
import config

### Configuration
avg_summary_plot = "mean"
base_prefix = "retained_h3k27me2"

WINDOW_CONFIGS = [
    (10000, 10000),
    #(5000, 5000),
    #(2000, 2000),
]

plot_type = "se"
color_maps = []
color_list = [
    "'white,#8c96c6'", "'white,#8c96c6'",
    "'white,#8c6bb1'", "'white,#8c6bb1'",
    "'white,#58135e'", "'white,#58135e'",
    "'white,#39107b'", "'white,#39107b'"
]
min_intensity = "0 0 0 0 0 0 0 0"
max_intensity = "25 25 100 100 15 15 25 25"
plot_height = 12
plot_width = 3.5
what_to_show = "plot, heatmap and colorbar"
x_label = "distance (bp)"
y_label = "sf x rpkm"
ref_point_label = "'peak\ncenter'"
sample_labels = ["'H3K27me2\nDMSO'", "'H3K27me2\nEZH2i'",
                "'H3K27me3\nDMSO'", "'H3K27me3\nEZH2i'",
                "'CBX2\nDMSO'", "'CBX2\nEZH2i'",
                "'EZH2\nDMSO'", "'EZH2\nEZH2i'"]
regions_labels = ["'Retained H3K27me2'"]
sort_using_samples = None
plot_title = None
y_min = None
y_max = max_intensity
legend_location = "none"
per_group = False
image_format = "svg"
###

def run_plot_heatmap(upstream_region: int, downstream_region: int):
    half_kb = upstream_region // 1000

    # Simple straight prefix, no clustering added
    output_suffix = f"{half_kb}kb"
    matrix_file = f"{base_prefix}_{half_kb}kb.gz"

    work_dir = config.DEEPTOOLSDIR
    matrix_dir = os.path.join(work_dir, "matrix")
    plot_dir = os.path.join(work_dir, "plots")

    output_prefix = f"{base_prefix}_{output_suffix}"
    matrix_path = os.path.join(matrix_dir, matrix_file)
    regions_dir = os.path.join(work_dir, "regions")

    optional_param = []
    if color_maps:
        optional_param.append(f"--colorMap {' '.join(color_maps)}")
    if color_list:
        optional_param.append(f"--colorList {' '.join(color_list)}")
    if sample_labels:
        optional_param.append(f'--samplesLabel {" ".join(sample_labels)}')
    if regions_labels:
        optional_param.append(f'--regionsLabel {" ".join(regions_labels)}')
    if min_intensity and "auto" not in min_intensity:
        optional_param.append(f"--zMin {min_intensity}")
    if max_intensity and "auto" not in max_intensity:
        optional_param.append(f"--zMax {max_intensity}")
    if x_label:
        optional_param.append(f"-x '{x_label}'")
    if y_label:
        optional_param.append(f"-y '{y_label}'")
    if y_min and "auto" not in y_min:
        optional_param.append(f"--yMin {y_min}")
    if y_max and "auto" not in y_max:
        optional_param.append(f"--yMax {y_max}")
    if per_group:
        optional_param.append("--perGroup")
    if plot_title:
        optional_param.append(f"--plotTitle '{plot_title}'")
    if sort_using_samples is not None:
        optional_param.append(f"--sortUsingSamples {sort_using_samples}")

    optional_params = " \\\n  ".join(optional_param)

    heatmap_file = os.path.join(plot_dir, f"{output_prefix}.svg")
    profile_file = os.path.join(plot_dir, f"{output_prefix}_profile.svg")
    
    job_name = f"{output_prefix}_plotHeatmap"
    log_dir = os.path.join(plot_dir, "log")
    log_out = os.path.join(log_dir, f"{output_prefix}_heatmap.log")
    log_err = os.path.join(log_dir, f"{output_prefix}_heatmap.error")
    batch_file = os.path.join(log_dir, f"{output_prefix}_heatmap.batch")

    os.makedirs(work_dir, exist_ok=True)
    os.makedirs(matrix_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(regions_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    batch_cmd = f"""\
#!/bin/bash
#BSUB -P {config.PROJECT_NAME}
#BSUB -J {job_name}
#BSUB -oo {log_out}
#BSUB -eo {log_err}
#BSUB -q {config.BSUB_QUEUE}
#BSUB -n 8
#BSUB -M 65536
#BSUB -R "rusage [mem=8192] span[hosts=1]"

source {config.HOMEDIR}/miniconda3/etc/profile.d/conda.sh
conda activate chipseq


plotHeatmap -m {matrix_path} \\
   -o {heatmap_file} \\
   --outFileSortedRegions {regions_dir}/{output_prefix}_sorted_regions.bed \\
   --plotType '{plot_type}' \\
   --averageTypeSummaryPlot '{avg_summary_plot}' \\
   --refPointLabel {ref_point_label} \\
   --whatToShow '{what_to_show}' \\
   --heatmapHeight {plot_height} \\
   --heatmapWidth {plot_width} \\
   {optional_params} \\
   --legendLocation {legend_location} \\
   --plotFileFormat {image_format}

"""

    with open(batch_file, 'w', encoding='utf-8') as f:
        f.write(batch_cmd)

    print(f"Submitting {job_name}")
    os.system(f"bsub < {batch_file}")
    time.sleep(1)

def main():
    parser = argparse.ArgumentParser(description="Generate heatmaps (without clustering)")
    parser.add_argument("--base_prefix", type=str, help="Overwrite base_prefix")
    parser.add_argument("--avg_summary", type=str, help="Overwrite avg_summary_plot")
    parser.add_argument("--min_intensity", type=str, help="Overwrite min_intensity")
    parser.add_argument("--max_intensity", type=str, help="Overwrite max_intensity")
    parser.add_argument("--ymin", type=str, help="Overwrite y_min")
    parser.add_argument("--ymax", type=str, help="Overwrite y_max")
    parser.add_argument("--x_label", type=str, help="Overwrite x_label")
    parser.add_argument("--y_label", type=str, help="Overwrite y_label")
    parser.add_argument("--samples_label", type=str, nargs="+", help="Overwrite sample_labels")
    parser.add_argument("--regions_label", type=str, help="Overwrite regions_labels")
    args = parser.parse_args()

    global base_prefix, avg_summary_plot, min_intensity, max_intensity, y_min, y_max, x_label, y_label, sample_labels, regions_labels
    
    if args.samples_label:
        sample_labels = [f"'{label}'" for label in args.samples_label]

    if args.avg_summary:
        avg_summary_plot = args.avg_summary
    
    if args.base_prefix:
        base_prefix = args.base_prefix
    elif args.avg_summary:
        base_prefix = "nuc_spread_enrichment"

    if args.min_intensity: min_intensity = args.min_intensity
    if args.max_intensity: 
        max_intensity = args.max_intensity
        y_max = max_intensity
    if args.ymin: y_min = args.ymin
    if args.ymax: y_max = args.ymax
    if args.x_label: x_label = args.x_label
    if args.y_label: y_label = args.y_label
    if args.regions_label: regions_labels = [f"'{args.regions_label}'"]

    # Log summary of configuration
    print(f"Average summary plot: {avg_summary_plot}")
    print(f"Base prefix: {base_prefix}")
    print(f"Min intensity: {min_intensity}")
    print(f"Max intensity: {max_intensity}")

    for upstream, downstream in WINDOW_CONFIGS:
        run_plot_heatmap(upstream, downstream)

if __name__ == "__main__":
    main()
