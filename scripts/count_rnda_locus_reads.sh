#!/bin/bash

P_I_BAM_DIR=$1
P_O_TABLE_FILE=$2

if [ -z "$P_I_BAM_DIR" ] || [ -z "$P_O_TABLE_FILE" ]; then
    echo "ERROR: Missing required arguments."
    echo "USAGE: $0 <input BAM directory> <output table file>"
    exit 1
fi

if [ ! -d "$P_I_BAM_DIR" ]; then
    echo "ERROR: Input BAM directory '$P_I_BAM_DIR' does not exist."
    exit 1
fi

if [ -f "$P_O_TABLE_FILE" ]; then
    echo "WARNING: Output table file '$P_O_TABLE_FILE' already exists and will be overwritten."
    rm -f "$P_O_TABLE_FILE"
fi

file_sep=$'\t'
bam_filter="_dupremoved.sorted.bam"
locus="chrXII:451300-469300"
# Initialize output table with header
echo "Library${file_sep}ReadCount" >> "$P_O_TABLE_FILE"

# Loop through each BAM file in the input directory
for bam_file_path in "$P_I_BAM_DIR"/*$bam_filter; do
    samtools index "$bam_file_path"  # Ensure BAM is indexed
    sample_name=$(basename "$bam_file_path" $bam_filter)
    # Count reads in the BAM file and append to the output table
    read_count=$(samtools view -c "$bam_file_path" "$locus")
    echo -e "$sample_name$file_sep$read_count" >> "$P_O_TABLE_FILE"
done