#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_bam_dir> <output_bam_dir> <regions_bed_file>"
    exit 1
fi

# Assign input arguments to variables
INPUT_BAM_DIR=$1
OUTPUT_BAM_DIR=$2
REGIONS_BED_FILE=$3

# Check if the input BAM directory exists
if [ ! -d "$INPUT_BAM_DIR" ]; then
    echo "Error: Input BAM directory does not exist."
    exit 1
fi

# Check if the BED file exists
if [ ! -f "$REGIONS_BED_FILE" ]; then
    echo "Error: BED file does not exist."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_BAM_DIR"

# Loop through all BAM files in the input directory
for BAM_FILE in "$INPUT_BAM_DIR"/*.bam; do
    if [ -f "$BAM_FILE" ]; then
        # Extract the base name of the BAM file
        BASENAME=$(basename "$BAM_FILE")
        OUTPUT_BAM_FILE="$OUTPUT_BAM_DIR/$BASENAME"

        # Use samtools to filter the BAM file
        samtools view -L "$REGIONS_BED_FILE" -b -U "$OUTPUT_BAM_FILE" "$BAM_FILE" >/dev/null

        if [ $? -eq 0 ]; then
            echo "Filtered $BAM_FILE -> $OUTPUT_BAM_FILE"
        else
            echo "Error filtering $BAM_FILE"
        fi
    fi
done

echo "Filtering complete."