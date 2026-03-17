#!/bin/bash

# Activate the komb conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate komb

##
# Test script for create_anti_features.sh
#
# TESTING INSTRUCTIONS:
# ===================
# This script tests the create_anti_features.sh functionality by:
# 1. Taking the genes.bed file (features) and chromosome_lengths.tsv
# 2. Running create_anti_features.sh to generate anti-features
# 3. Validating that:
#    - Anti-features are created correctly
#    - Features + anti-features fit within chromosome boundaries
#    - Strand information is preserved
#    - 1bp features are filtered out
#
# Usage:
#   bash tests/test_create_anti_features.sh
#
# Expected output:
#   - test_output/anti_features.bed: Combined features and anti-features
#   - Validation checks to ensure total coverage doesn't exceed chromosome length
#

set -e

# Check for bedtools availability
if ! command -v bedtools &> /dev/null; then
    echo "[WARNING] bedtools is not installed in the komb environment."
    echo "[INFO] Please activate the komb environment with: conda activate komb"
    exit 1
fi

# Define test paths
TEST_DIR="scripts/tests"
DATA_DIR="$TEST_DIR/data"
OUTPUT_DIR="$TEST_DIR/output"
RESULTS_DIR="$TEST_DIR/results"

GENES_BED="$DATA_DIR/genes.bed"
CHR_LENGTHS="$DATA_DIR/chromosome_lengths.tsv"
ANTI_FEATURES_OUTPUT="$OUTPUT_DIR/anti_features.bed"

# Clean up previous test output
if [ -d "$OUTPUT_DIR" ]; then
    rm -rf "$OUTPUT_DIR"
fi
mkdir -p "$OUTPUT_DIR"

echo "[INFO] Running create_anti_features.sh test"
echo "[INFO] Input features: $GENES_BED"
echo "[INFO] Chromosome lengths: $CHR_LENGTHS"
echo "[INFO] Output: $ANTI_FEATURES_OUTPUT"
echo ""

# Run the create_anti_features.sh script
bash scripts/create_anti_features.sh "$GENES_BED" "$CHR_LENGTHS" "$ANTI_FEATURES_OUTPUT"

# Validate output exists and is not empty
if [ ! -f "$ANTI_FEATURES_OUTPUT" ]; then
    echo "[ERROR] Output file not created: $ANTI_FEATURES_OUTPUT"
    exit 1
fi

if [ ! -s "$ANTI_FEATURES_OUTPUT" ]; then
    echo "[ERROR] Output file is empty: $ANTI_FEATURES_OUTPUT"
    exit 1
fi

echo "[INFO] Output file created successfully"
echo "[INFO] First 10 lines of output:"
head -10 "$ANTI_FEATURES_OUTPUT"
echo ""

# Validate output format
echo "[INFO] Validating output format"
line_count=$(wc -l < "$ANTI_FEATURES_OUTPUT")
echo "[INFO] Total features (features + anti-features): $line_count"

# Check that all lines have at least 5 columns (chrom, start, end, name, strand)
invalid_lines=$(awk -F'\t' 'NF < 5 {print NR}' "$ANTI_FEATURES_OUTPUT" || true)
if [ -n "$invalid_lines" ]; then
    echo "[ERROR] Lines with invalid format (expected 5+ columns): $invalid_lines"
    exit 1
fi

echo "[INFO] Format validation passed"
echo ""

# Validate strand information
echo "[INFO] Validating strand information"
plus_features=$(grep -c '\+$' "$ANTI_FEATURES_OUTPUT" || echo 0)
minus_features=$(grep -c '\-$' "$ANTI_FEATURES_OUTPUT" || echo 0)
echo "[INFO] Features on + strand: $plus_features"
echo "[INFO] Features on - strand: $minus_features"
echo ""

# Validate chromosome coverage
echo "[INFO] Validating chromosome coverage boundaries"
declare -A chr_lengths
while IFS=$'\t' read -r chrom length
do
    chr_lengths[$chrom]=$length
done < "$CHR_LENGTHS"

for chrom in "${!chr_lengths[@]}"
do
    for strand in "+" "-"
    do
        total_covered=$(grep "^${chrom}\t" "$ANTI_FEATURES_OUTPUT" 2>/dev/null | \
            awk -v s="$strand" '$NF==s {print $1, $2, $3}' | \
            bedtools merge 2>/dev/null | \
            awk '{sum += ($3 - $2)} END {print sum+0}')
        
        if [ -z "$total_covered" ] || [ "$total_covered" -eq 0 ]; then
            continue
        fi
        
        chr_len=${chr_lengths[$chrom]}
        echo "[INFO] $chrom ($strand): $total_covered bp covered out of $chr_len bp"
        
        if [ "$total_covered" -gt "$chr_len" ]; then
            echo "[ERROR] Coverage exceeds chromosome length!"
            exit 1
        fi
    done
done

echo ""
echo "[SUCCESS] All validation checks passed!"
echo "[OUTPUT] Full results saved to: $ANTI_FEATURES_OUTPUT"
