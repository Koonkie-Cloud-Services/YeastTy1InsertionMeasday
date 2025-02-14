#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "[ERROR] unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <feature BED> <chromosome lengths> <output BED>"
    exit
fi

output_bed=$3

if [ -f $output_bed ]; then
	rm $output_bed
fi

cat $1 | \
	grep '\+$' | \
	bedtools complement -i - -g $2 | \
	awk -v id="anti-feature" -v strand="+" '{ OFS="\t"; print $0,id,strand }' \
	>> "$output_bed.tmp"

cat $1 | \
	grep '\-$' | \
	bedtools complement -i - -g $2 | \
	awk -v id="anti-feature" -v strand="-" '{ OFS="\t"; print $0,id,strand }' \
	>> "$output_bed.tmp"

# Sort the features by their position, remove anti-features with a length of 1 bp
bedtools sort -g $2 -i "$output_bed.tmp" | \
	awk -F"\t" '{ OFS="\t"; if ($3-$2 > 1) print $0 }' \
	> $output_bed
rm "$output_bed.tmp"

echo "[INFO] Finished creating anti-features file"
