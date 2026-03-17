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
    if [ -f $OUTPUT_FILE ]; then
        echo "Skipping $BAM_BASE"
        return
    fi
    
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

# Function to process SAM files and extract adjacent sequences
process_adjacent_sequences() {
    local SAM_FILE="$1"
    local BED_FILE="$2"
    local OUTPUT_DIR="$3"
    
    local SAM_BASE=$(basename "$SAM_FILE" .sam)
    local TSV_OUTPUT="$OUTPUT_DIR/${SAM_BASE}_adjacent_sequences.tsv"
    
    echo "Processing adjacent sequences for: $SAM_BASE"
    
    # Create TSV header
    echo -e "read_name\talignment_start\ttrna_gene\tadjacent_sequence\tadjacent_sequence_length" > "$TSV_OUTPUT"
    
    # Load tRNA regions into associative arrays for quick lookup
    declare -A trna_starts trna_ends trna_names trna_chrs
    while IFS=$'\t' read -r chr start end name score strand; do
        if [[ "$chr" != "chr" ]]; then  # Skip header
            local key="${chr}:${start}-${end}"
            trna_starts["$key"]=$start
            trna_ends["$key"]=$end
            trna_names["$key"]=$name
            trna_chrs["$key"]=$chr
        fi
    done < "$BED_FILE"
    
    # Process each alignment in SAM file
    grep -v '^@' "$SAM_FILE" | while IFS=$'\t' read -r qname flag rname pos mapq cigar rnext pnext tlen seq qual rest; do
        # Parse alignment start position (convert from 1-based to 0-based for calculations)
        local align_start=$((pos - 1))
        local align_end=$((align_start + ${#seq}))
        
        # Find overlapping tRNA gene
        local overlapping_gene=""
        local gene_start=0
        local gene_end=0
        
        for key in "${!trna_chrs[@]}"; do
            if [[ "${trna_chrs[$key]}" == "$rname" ]]; then
                local g_start=${trna_starts[$key]}
                local g_end=${trna_ends[$key]}
                
                # Check for overlap (even partial)
                if [[ $align_start -lt $g_end && $align_end -gt $g_start ]]; then
                    overlapping_gene=${trna_names[$key]}
                    gene_start=$g_start
                    gene_end=$g_end
                    break
                fi
            fi
        done
        
        # Calculate adjacent sequence
        local adjacent_seq=""
        if [[ -n "$overlapping_gene" ]]; then
            # Calculate overlap boundaries relative to read
            local overlap_start_in_read=$((gene_start > align_start ? gene_start - align_start : 0))
            local overlap_end_in_read=$((gene_end < align_end ? gene_end - align_start : ${#seq}))
            
            # Extract sequences before and after the overlap
            local upstream_seq=""
            local downstream_seq=""
            
            if [[ $overlap_start_in_read -gt 0 ]]; then
                upstream_seq=${seq:0:$overlap_start_in_read}
            fi
            
            if [[ $overlap_end_in_read -lt ${#seq} ]]; then
                downstream_seq=${seq:$overlap_end_in_read}
            fi
            
            adjacent_seq="${upstream_seq};${downstream_seq}"
            len_adj=$(( ${#adjacent_seq} - 1 ))
        else
            # No overlap found, keep entire sequence
            adjacent_seq="$seq"
            overlapping_gene="NO_OVERLAP"
            len_adj=${#adjacent_seq}
        fi
        
        # Output to TSV
        echo -e "$qname\t$pos\t$overlapping_gene\t$adjacent_seq\t$len_adj" >> "$TSV_OUTPUT"
    done
    
    echo "Adjacent sequences processed: $TSV_OUTPUT"
}

# Export function for parallel
export -f process_adjacent_sequences

# Run extraction in parallel
ls $BAM_PATTERN | parallel -j $NUM_JOBS extract_alignments {} $OUTPUT_DIR $BED_FILE

echo "All extractions completed!"

# Process adjacent sequences for all SAM files
echo "Processing adjacent sequences..."
ls "$OUTPUT_DIR"/*_trna_alignments.sam | parallel -j $NUM_JOBS process_adjacent_sequences {} $BED_FILE $OUTPUT_DIR

echo "Adjacent sequence processing completed!"
echo "Output files are in: $OUTPUT_DIR"



# Calculate and print average percentage of specific sequences in the alignments
echo "Calculating average percentage of specific sequences in alignments..."
MOTIF_SUMMARY_FILE="$OUTPUT_DIR/motif_summary.tsv"
MOTIF="AACATGG"
# Create reverse complement of MOTIF
MOTIF_RC=$(echo "$MOTIF" | tr 'ATCG' 'TAGC' | rev)
ls $OUTPUT_DIR/*_trna_alignments.sam > /dev/null 2>&1 || { echo "No alignment files found for motif calculation."; exit 1; }
echo -e "Filename\tTotal_Alignments\tPercentage_$MOTIF" > "$MOTIF_SUMMARY_FILE"
for f in $OUTPUT_DIR/*_trna_alignments.sam; do
    fname=$(basename "$f")
    echo "Processing motif calculation for $fname... "

    lines=$( grep -cv "^@" $f )
    # grep "$MOTIF\|$MOTIF_RC" -c $f || true
    hits=$(grep "$MOTIF\|$MOTIF_RC" -c $f || true)
    if [ $hits -eq 0 ] || [ $lines -eq 0 ]; then
        perc_hits=0
    else
        perc_hits=$(echo "scale=4;100 * $hits / $lines" | bc)
    fi
    echo -e "$fname\t$lines\t$perc_hits" >> "$MOTIF_SUMMARY_FILE"
    echo "Completed motif calculation for $fname: $hits / $lines"
done

echo -e "\nAverage percentage of $MOTIF and its reverse complement across all alignments:"
awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' "$MOTIF_SUMMARY_FILE"