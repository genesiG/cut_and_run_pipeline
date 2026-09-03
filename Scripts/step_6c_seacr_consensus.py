#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
step_6c_seacr_consensus.py  —  Consensus peak generation for SEACR outputs.

Strategy
--------
For each group (e.g., K27M_K27me3):
  1. Find all replicate relaxed BED files (from step_6b).
  2. Append the replicate identifier (e.g., 'A') to column 4 of each BED file.
  3. Concatenate and sort them.
  4. Use `bedtools merge -c 4 -o count_distinct` to merge overlapping peaks
     and count the number of distinct replicates supporting each merged region.
  5. Keep only those supported by >= 2 replicates (Majority rules).
  
Output: Analysis_Data/peaks/seacr/<group_label>_consensus.bed
"""

import os
import glob
import subprocess
import re

import config

SEACR_OUTDIR = config.SEACR_OUTDIR

# The groups we care about
GROUPS = [
    "K27M_DMSO_K27me3",
    "K27M_500nM_K27me3",
    "K27MKO_DMSO_K27me3",
    "K27MKO_500nM_K27me3",
]

def main():
    os.makedirs(SEACR_OUTDIR, exist_ok=True)
    
    for group in GROUPS:
        # Find all replicate relaxed beds for this group
        # Pattern: {group}_*_R1_seacr.relaxed.bed
        # Actually in step_6b we named them: clean_sample(bam_name) + ".relaxed.bed"
        # e.g., K27M_K27me3_A_R1.relaxed.bed
        search_pattern = os.path.join(SEACR_OUTDIR, f"{group}_*.relaxed.bed")
        rep_beds = sorted(glob.glob(search_pattern))
        
        if not rep_beds:
            print(f"No BED files found for group {group}")
            continue
            
        print(f"\nProcessing consensus for {group} ({len(rep_beds)} replicates)")
        
        all_bed = os.path.join(SEACR_OUTDIR, f"{group}_all_reps.bed")
        sorted_bed = os.path.join(SEACR_OUTDIR, f"{group}_all_reps.sorted.bed")
        merged_bed = os.path.join(SEACR_OUTDIR, f"{group}_merged.bed")
        consensus_bed = os.path.join(SEACR_OUTDIR, f"{group}_consensus.bed")
        
        # 1. Combine and add replicate ID to column 4
        with open(all_bed, 'w') as out_f:
            for bed in rep_beds:
                # Extract replicate ID, e.g., 'A' from 'K27M_K27me3_A_R1.relaxed.bed'
                basename = os.path.basename(bed)
                # Try to extract the replicate letter (A, B, or C)
                match = re.search(r'_([A-Z])_R1', basename)
                rep_id = match.group(1) if match else basename
                
                with open(bed, 'r') as in_f:
                    for line in in_f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                parts[3] = rep_id
                                out_f.write('\t'.join(parts) + '\n')
                                
        # 2. Sort
        print("  Sorting...")
        subprocess.run(f"sort -k1,1 -k2,2n {all_bed} > {sorted_bed}", shell=True, check=True)
        
        # 3. Merge and count distinct replicates
        print("  Merging and counting...")
        # bedtools merge requires BED3. We merge based on coordinates, and aggregate column 4
        cmd = f"bedtools merge -i {sorted_bed} -c 4 -o count_distinct > {merged_bed}"
        subprocess.run(cmd, shell=True, check=True)
        
        # 4. Filter for >= 2
        print("  Filtering for consensus (>= 2 reps)...")
        passed = 0
        total = 0
        with open(merged_bed, 'r') as in_f, open(consensus_bed, 'w') as out_f:
            for line in in_f:
                total += 1
                parts = line.strip().split('\t')
                count = int(parts[3])
                if count >= 2:
                    out_f.write(line)
                    passed += 1
                    
        print(f"  Consensus peaks: {passed} / {total}")
        
        # Cleanup
        os.remove(all_bed)
        os.remove(sorted_bed)
        os.remove(merged_bed)

if __name__ == "__main__":
    main()
