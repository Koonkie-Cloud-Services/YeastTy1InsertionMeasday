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

if [ "$#" -ne 3 ]; then
    log ERROR "unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <feature BED> <chromosome lengths> <output BED>"
    exit
fi

feature_bed=$1
chrom_lengths=$2
output_bed=$3

if [ -f $output_bed ]; then
	rm $output_bed
fi

cat $feature_bed | \
	grep '\+$' | \
	bedtools sort -g $chrom_lengths -i - | \
	bedtools complement -i - -g $chrom_lengths | \
	awk -v id="anti-feature" -v strand="+" '{ OFS="\t"; print $0,id,strand }' \
	>> "$output_bed.tmp"

cat $feature_bed | \
	grep '\-$' | \
	bedtools sort -g $chrom_lengths -i - | \
	bedtools complement -i - -g $chrom_lengths | \
	awk -v id="anti-feature" -v strand="-" '{ OFS="\t"; print $0,id,strand }' \
	>> "$output_bed.tmp"

# Sort the features by their position, remove anti-features with a length of 1 bp
bedtools sort -g $chrom_lengths -i "$output_bed.tmp" | \
	awk -F"\t" '{ OFS="\t"; if ($3-$2 > 1) print $0 }' \
	> $output_bed
rm "$output_bed.tmp"

log INFO "Finished creating anti-features file"

bash scripts/validate_anti_feature_coords.sh "$chrom_lengths" "$output_bed"