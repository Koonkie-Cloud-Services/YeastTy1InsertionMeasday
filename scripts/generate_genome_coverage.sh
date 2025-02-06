#!/bin/bash

sample_name_map=$1
output_file_name=$2

if [ "$#" -ne 2 ]; then
    echo "ERROR: unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <sample_name_map> <output_file_name>"
    exit
fi

# Define directory paths
working_dir="data"
bam_suffix="dupremoved.sorted"
bam_dir="$working_dir/BAM_alignments_to_yeast_reference"
bed_dir="$working_dir/BEDs/"
features_file="$working_dir/feature_locations/all_features.bed"


if [ ! -d $working_dir ]; then
	echo "ERROR: working directory '$working_dir' doesn't exist"
	exit
fi

echo "Rename BAM files"
while IFS=$'\t' read -r sampleId strainId sugar uniqueId prefix fastqName
do
	# Skip the header line
	if [[ $sampleId == "sampleId" ]]; then
		continue
	fi

	old_bam="${bam_dir}/${prefix}_${bam_suffix}.bam"
	new_bam="${bam_dir}/${uniqueId}_${bam_suffix}.bam"
	if [ ! -f $old_bam ]; then
		echo "ERROR: Unable to find input BAM file '$old_bam'"
		exit
	fi

	cp $old_bam $new_bam
done<$sample_name_map

echo "Create BED files"
mkdir -p $bed_dir
while IFS=$'\t' read -r sampleId strainId sugar uniqueId prefix fastqName
do
        # Skip the header line
        if [[ $sampleId == "sampleId" ]]; then
                continue
        fi

	echo $uniqueId
	bam_file="${bam_dir}/${uniqueId}_${bam_suffix}.bam"
	bedtools coverage \
		-a $features_file \
		-b $bam_file \
		-bed -sorted -counts > \
		"${bed_dir}/${uniqueId}_${bam_suffix}_feature_coverage.bed"
done<$sample_name_map
