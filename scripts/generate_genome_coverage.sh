#!/bin/bash

output_file_path=$1

if [ "$#" -ne 1 ]; then
    echo "ERROR: unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <output_file_name>"
    exit
fi

# Define directory paths
working_dir="data"
bam_suffix="dupremoved.sorted"
bam_dir="$working_dir/BAM_alignments_to_yeast_reference"
bed_dir="$working_dir/BEDs/"
sample_name_map="$working_dir/hiseq_sample_id_map.tsv"
features_file="$working_dir/feature_locations/all_features.bed"
chr_lengths_file="$working_dir/feature_locations/chromosome_lengths.tsv"


if [ ! -d $working_dir ]; then
	echo "ERROR: working directory '$working_dir' doesn't exist"
	exit
fi

if [ -f $output_file_path ]; then
	echo "ERROR: remove output file '$output_file_path' to prevent overwriting"
	exit
fi

mkdir -p $bed_dir
touch $output_file_path

while IFS=$'\t' read -r sampleId strainId sugar uniqueId prefix fastqName || [ -n "$fastqName" ]
do
	# Skip the header line
	if [[ $sampleId == "sampleId" ]]; then
		continue
	fi

	# Define file paths for sample
	old_bam="${bam_dir}/${prefix}_${bam_suffix}.bam"
	new_bam_path="${bam_dir}/${uniqueId}_${bam_suffix}.bam"
	bed_path="${bed_dir}/${uniqueId}_${bam_suffix}_feature_coverage.bed"

	if [ ! -f $old_bam ]; then
		echo "ERROR: Unable to find input BAM file '$old_bam'"
		exit
	fi

	echo "[INFO] Renaming $uniqueId BAM: '$old_bam'"
	cp $old_bam $new_bam_path

	echo "[INFO] calculating feature coverage for sample: $uniqueId"
	bedtools coverage \
		-a $features_file \
		-b $new_bam_path \
		-g $chr_lengths_file \
		-bed -sorted -counts \
		> $bed_path

	# Write to the output file path
	cat $bed_path | \
		awk -v id="$uniqueId" '{ OFS="\t"; print id,$0 }' >> \
		$output_file_path
done<$sample_name_map
