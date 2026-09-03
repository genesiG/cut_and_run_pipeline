#!/usr/bin/env python3

"""
prepare_samples.py  - Script to automate sample metadata file creation
                    - Run this script to generate the samples.txt and 
                      sample_metadata.txt files that will be used throughout the analysis

INSTRUCTIONS
    - Fastq files need to be in the Original_Data_and_Metadata/Original_Data/fastq_files folder
    - Run the script in the terminal as "python3 step_0c_prepare_samples_file.py"
"""

import os
import re
import config

### Set up variables
data_dir = config.DATADIR
metadata = config.METADATA
seq_extension = config.EXTENSIONS
samples_file = os.path.join(metadata, "samples.txt")
paired_samples = os.path.join(metadata, "paired_samples.txt")
sample_metadata_file = os.path.join(metadata, "sample_metadata.txt")
###

### Antibody and group patterns are defined in config.py (ANTIBODY_PATTERNS, GROUP_PATTERNS)
# Edit config.py to add new antibodies or treatment groups — do not edit here.

def extract_antibody(sample_name):
    """Extract antibody name from sample name using patterns from config.py"""
    for pattern, antibody in config.ANTIBODY_PATTERNS.items():
        if pattern in sample_name:
            return antibody
    return 'Unknown'

def extract_group(sample_name):
    """Extract treatment group from sample name using patterns from config.py"""
    for pattern, group in config.GROUP_PATTERNS.items():
        if pattern in sample_name:
            return group
    return 'Unknown'

### Prepare sample.txt file
# Get a list of all files in the directory
file_list = os.listdir(data_dir)

# Filter the list to include only files (not subdirectories)
file_list = [f for f in file_list if os.path.isfile(os.path.join(data_dir, f))]

# Extract sample names from file names while excluding the custom extensions
sample_names = []
for f in file_list:
    for ext in seq_extension:
        if f.endswith(ext):
            sample_name = f[:-len(ext)]  # Remove the extension from the filename
            sample_names.append(sample_name)
            break  # Stop checking extensions once a match is found

sample_names = sorted(sample_names)

# Create the samples.txt file
with open(samples_file, "w", encoding='utf-8') as file:
    # Write the header
    file.write("sample.name\tsample.control\n")

    # Write sample names and "-"
    for sample_name in sample_names:
        file.write(f"{sample_name}\t-\n")

print(f"samples.txt has been created at {metadata}")
###

### Prepare paired_sample.txt file
if config.IS_PAIRED_END:
    # Create a dictionary to store sample names and their corresponding paired samples
    sample_mapping = {}
    for sample_name in sample_names:
        # Check if this is a reverse read - skip if so
        is_reverse = False
        for rev_pattern in config.REVERSE_PATTERN:
            if rev_pattern in sample_name:
                is_reverse = True
                break

        if is_reverse:
            continue  # skip any reverse read file

        # Find the matching forward pattern and create paired sample name
        for fwd_pattern, rev_pattern in zip(config.FORWARD_PATTERN, config.REVERSE_PATTERN):
            if fwd_pattern in sample_name:
                paired_sample = sample_name.replace(fwd_pattern, rev_pattern)
                sample_mapping[sample_name] = paired_sample
                break

    # Create the paired_samples.txt file
    with open(paired_samples, "w", encoding='utf-8') as file:
        # Write the header
        file.write("sample.name\tsample.control\n")

        # Write paired samples
        for sample_name, paired_sample in sample_mapping.items():
            file.write(f"{sample_name}\t{paired_sample}\n")

    print(f"paired_samples.txt has been created at {metadata}")
###

### Prepare sample_metadata.txt file
# Filter to keep only R1 (forward reads) for paired-end or all samples for single-end
metadata_samples = []
if config.IS_PAIRED_END:
    for sample_name in sample_names:
        # Check if this is a forward read
        is_forward = False
        for fwd_pattern in config.FORWARD_PATTERN:
            if fwd_pattern in sample_name:
                is_forward = True
                break

        if is_forward:
            metadata_samples.append(sample_name)
else:
    metadata_samples = sample_names

# Create the sample_metadata.txt file
with open(sample_metadata_file, "w", encoding='utf-8') as file:
    # Write the header
    file.write("ID\tANTIBODY\tGROUP\n")

    # Write sample metadata
    for sample_name in metadata_samples:
        antibody = extract_antibody(sample_name)
        group = extract_group(sample_name)
        file.write(f"{sample_name}\t{antibody}\t{group}\n")

print(f"sample_metadata.txt has been created at {metadata}")

# Print summary
print("\n=== Summary ===")
print(f"Total samples: {len(sample_names)}")
if config.IS_PAIRED_END:
    print(f"Forward reads (in metadata): {len(metadata_samples)}")
    print(f"Paired samples: {len(sample_mapping)}")
else:
    print(f"Samples in metadata: {len(metadata_samples)}")

# Print antibody distribution
antibody_counts = {}
for sample in metadata_samples:
    ab = extract_antibody(sample)
    antibody_counts[ab] = antibody_counts.get(ab, 0) + 1

print("\nAntibody distribution:")
for ab, count in sorted(antibody_counts.items()):
    print(f"  {ab}: {count}")

# Print group distribution
group_counts = {}
for sample in metadata_samples:
    grp = extract_group(sample)
    group_counts[grp] = group_counts.get(grp, 0) + 1

print("\nGroup distribution:")
for grp, count in sorted(group_counts.items()):
    print(f"  {grp}: {count}")

# Warn about unknowns
unknown_antibodies = [s for s in metadata_samples if extract_antibody(s) == 'Unknown']
if unknown_antibodies:
    print("\n⚠️  WARNING: The following samples have unknown antibodies:")
    for s in unknown_antibodies:
        print(f"    {s}")
    print("  Please update ANTIBODY_PATTERNS in this script.")

unknown_groups = [s for s in metadata_samples if extract_group(s) == 'Unknown']
if unknown_groups:
    print("\n⚠️  WARNING: The following samples have unknown groups:")
    for s in unknown_groups:
        print(f"    {s}")
    print("  Please update GROUP_PATTERNS in this script.")
###
