#!/bin/bash

##              ##
# LOGGING FUNCTION #
##              ##
log() {
    local level=$1
    shift
    local message="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case "${level^^}" in
        INFO)
            echo "[$timestamp] [INFO] $message"
            ;;
        WARNING|WARN)
            echo "[$timestamp] [WARNING] $message" >&2
            ;;
        ERROR)
            echo "[$timestamp] [ERROR] $message" >&2
            ;;
        *)
            echo "[$timestamp] [UNKNOWN] $message"
            ;;
    esac
}
    
##                                    ##
# VALIDATE FEATURE COVERAGE BOUNDARIES #
##                                    ##

usage() {
    echo "Usage: $0 <chrom_lengths.tsv> <features.bed>" >&2
}

if [ "$#" -ne 2 ]; then
    log ERROR "incorrect number of arguments provided"
    usage
    exit 1
fi

chrom_lengths="$1"
bed_path="$2"

log INFO "validating that features and anti-features do not exceed chromosome lengths"

# Check for bedtools availability
if ! command -v bedtools &> /dev/null; then
    log WARNING "bedtools was not found in the current environment. Activating 'komb'."
    # Activate the komb conda environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate komb
fi

# Read chromosome lengths into an associative array
declare -A chr_lengths
while IFS=$'\t' read -r chrom length
do
	chr_lengths[$chrom]=$length
done < $chrom_lengths

# Use bedtools merge to properly calculate total covered length per chromosome-strand
for chrom in "${!chr_lengths[@]}"
do
	for strand in "+" "-"
	do
		# Extract features for this chromosome and strand, merge overlaps, and sum
		total_covered=$(grep "^${chrom}\t" $bed_path | \
			awk -v s="$strand" '$NF==s {print $1, $2, $3}' | \
			bedtools merge 2>/dev/null | \
			awk '{sum += ($3 - $2)} END {print sum+0}')
		
		if [ -z "$total_covered" ] || [ "$total_covered" -eq 0 ]; then
			continue
		fi
		
		chr_len=${chr_lengths[$chrom]}
		if [ -z "$chr_len" ]; then
			log ERROR "chromosome '$chrom' not found in chromosome lengths file"
			exit 1
		fi
		
		if [ "$total_covered" -gt "$chr_len" ]; then
			log ERROR "Features and anti-features exceed chromosome length: ${chrom} (${strand}) - covered: $total_covered, chromosome length: $chr_len"
			exit 1
		fi
	done
done

log INFO "Feature boundary validation passed"