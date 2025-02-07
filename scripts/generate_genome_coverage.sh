#!/bin/bash

if [ "$#" -ne 0 ]; then
    echo "ERROR: unexpected number of arguments provided"
    echo -e "USAGE:\n$0"
    exit
fi

##                          ##
# DEFINE ANALYSIS PARAMETERS #
##                          ##
# Define directory paths
p_working_dir="data"
p_input_dir="$p_working_dir/input"
p_output_dir="$p_working_dir/output"

# Analysis parameters
p_bam_suffix="dupremoved.sorted"

# Input paths
p_i_bam_dir="$p_input_dir/BAMs"
p_i_features_dir="$p_input_dir/feature_locations"
p_i_sample_name_map="$p_input_dir/hiseq_sample_id_map.tsv"
p_i_features_file="$p_i_features_dir/all_features.bed"
p_i_trna_locations_file="$p_i_features_dir/trna_regions.edited.bed"
p_i_chr_lengths_file="$p_i_features_dir/chromosome_lengths.tsv"

# Output paths
p_o_bam_dir="$p_output_dir/BAMs"
p_o_bed_dir="$p_output_dir/BEDs"
p_o_features_dir="$p_output_dir/feature_locations"
p_o_features_file="$p_o_features_dir/all_features.bed"
p_o_pos_features_file="$p_o_features_dir/all_positive_features.bed"
p_o_neg_features_file="$p_o_features_dir/all_negative_features.bed"
p_o_anti_features_file="$p_o_features_dir/all_anti-features.bed"
p_o_pos_anti_features_file="$p_o_features_dir/all_positive_anti-features.bed"
p_o_neg_anti_features_file="$p_o_features_dir/all_negative_anti-features.bed"
p_o_trna_locations_file="$p_o_features_dir/trna_regions.bed"
p_o_coverage_file="$p_output_dir/genome_coverage_counts_${p_bam_suffix}.bed"


##                    ##
# PARAMETER VALIDATION #
##                    ##

if [ ! -d $working_dir ]; then
	echo "[ERROR] working directory '$working_dir' doesn't exist"
	exit
fi

if [ -f $p_o_coverage_file ]; then
	echo "[WARNING] removing output file '$p_o_coverage_file'"
	rm $p_o_coverage_file
        touch $p_o_coverage_file
fi

for dir_path in $p_output_dir $p_o_bed_dir $p_o_bam_dir $p_o_features_dir
do
	if [ ! -d $dir_path ]; then
		mkdir $dir_path
	fi
done

##              ##
# BEGIN ANALYSIS #
##              ##

echo "[INFO] converting arabic to roman chromosome names"
python scripts/arabic_to_roman.py $p_i_features_file $p_o_features_file
python scripts/arabic_to_roman.py $p_i_trna_locations_file $p_o_trna_locations_file


echo "[INFO] creating strand-wise feature files"
cat $p_o_features_file |
	grep '\+$' > $p_o_pos_features_file
cat $p_o_features_file |
	grep '\-$' > $p_o_neg_features_file

echo "[INFO] creating anti-features file '$p_o_anti_features_file'"
scripts/create_anti_features.sh $p_o_features_file $p_i_chr_lengths_file $p_o_anti_features_file

echo "[INFO] creating strand-wise anti-feature files"
cat $p_o_anti_features_file |
	grep '\+$' > $p_o_pos_anti_features_file
cat $p_o_anti_features_file |
	grep '\-$' > $p_o_neg_anti_features_file


while IFS=$'\t' read -r sampleId strainId sugar uniqueId prefix fastqName || [ -n "$fastqName" ]
do
	# Skip the header line
	if [[ $sampleId == "sampleId" ]]; then
		continue
	fi

	# Define file paths for sample
	old_bam="${p_i_bam_dir}/${prefix}_${p_bam_suffix}.bam"
	new_bam_path="${p_o_bam_dir}/${uniqueId}_${p_bam_suffix}.bam"
	bed_path="${p_o_bed_dir}/${uniqueId}_${p_bam_suffix}_feature_coverage.bed"

	if [ ! -f $old_bam ]; then
		echo "[ERROR] Unable to find input BAM file '$old_bam'"
		exit
	fi

	echo "[INFO] Renaming $uniqueId BAM: '$old_bam'"
	cp $old_bam $new_bam_path

	if [ -f $bed_path ]; then
		rm $bed_path
	fi
	touch $bed_path

	echo "[INFO] calculating feature coverage for sample: $uniqueId"
	# positive-features
	bedtools bamtobed -i $new_bam_path | \
		grep '\+$'  | \
		bedtools coverage \
		-a $p_o_pos_features_file \
		-b - \
		-g $p_i_chr_lengths_file \
		-sorted -counts \
		>> $bed_path

	# negative-features
	bedtools bamtobed -i $new_bam_path | \
		grep '\-$'  | \
		bedtools coverage \
		-a $p_o_neg_features_file \
		-b - \
		-g $p_i_chr_lengths_file \
		-sorted -counts \
		>> $bed_path

	# positive-antifeatures
	bedtools bamtobed -i $new_bam_path | \
		grep '\+$'  | \
		bedtools coverage \
		-a $p_o_pos_anti_features_file \
		-b - \
		-g $p_i_chr_lengths_file \
		-sorted -counts \
		>> $bed_path

	# negative-antifeatures
	bedtools bamtobed -i $new_bam_path | \
		grep '\-$'  | \
		bedtools coverage \
		-a $p_o_neg_anti_features_file \
		-b - \
		-g $p_i_chr_lengths_file \
		-sorted -counts \
		>> $bed_path

	# Write to the output file path
	cat $bed_path | \
		awk -v id="$uniqueId" '{ OFS="\t"; print id,$0 }' >> \
		$p_o_coverage_file
done<$p_i_sample_name_map

echo "[INFO] Genome coverage information written to '$p_o_coverage_file'"
echo "[INFO] Workflow completed"
