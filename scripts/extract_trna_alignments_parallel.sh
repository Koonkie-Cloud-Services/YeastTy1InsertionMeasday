#!/bin/bash

# Parallelized script to extract alignments from multiple BAM files for tRNA regions
# Requires GNU parallel: sudo apt-get install parallel
# Usage: ./extract_trna_alignments_parallel.sh [output_directory] [num_jobs]

set -euo pipefail

# Set default parameters
OUTPUT_DIR="${1:-data/output/trna_alignments}"
NUM_JOBS="${2:-$(nproc)}"  # Default to number of CPU cores
BED_FILE="data/output/feature_locations/trna_regions.bed"
BAM_PATTERN="data/output/alignments/*_dupremoved_uniq.sorted.bam"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Starting parallelized extraction of tRNA alignments..."
echo "Using $NUM_JOBS parallel jobs"

# Function to extract alignments for one BAM file
extract_alignments() {
    local BAM_FILE="$1"
    local OUTPUT_DIR="$2"
    local BED_FILE="$3"
    
    local BAM_BASE=$(basename "$BAM_FILE" .bam)
    local OUTPUT_FILE="$OUTPUT_DIR/${BAM_BASE}_trna_alignments.sam"
    
    echo "Processing: $BAM_BASE"
    
    # Write SAM header
    samtools view -H "$BAM_FILE" > "$OUTPUT_FILE"
    
    # Extract alignments for all regions and append
    tail -n +2 "$BED_FILE" | while IFS=$'\t' read -r chr start end name score strand; do
        samtools view "$BAM_FILE" "${chr}:${start}-${end}" >> "$OUTPUT_FILE" 2>/dev/null || true
    done
    
    # Count alignments
    local ALIGNMENT_COUNT=$(grep -v '^@' "$OUTPUT_FILE" | wc -l)
    echo "Completed: $BAM_BASE with $ALIGNMENT_COUNT alignments"
}

# Export function for parallel
export -f extract_alignments

# Run extraction in parallel
ls $BAM_PATTERN | parallel -j $NUM_JOBS extract_alignments {} $OUTPUT_DIR $BED_FILE

echo "All extractions completed!"
echo "Output files are in: $OUTPUT_DIR"