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
p_tmp_dir="$p_working_dir/.tmp"

# Analysis parameters
p_bam_suffix="dupremoved_uniq.sorted"

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
p_o_trna_locations_file="$p_o_features_dir/trna_regions.bed"
p_o_coverage_file="$p_output_dir/genome_coverage_counts_${p_bam_suffix}.bed"

# Intermediate file paths
t_features_file="$p_tmp_dir/features.bed"
t_pos_features_file="$p_tmp_dir/positive_features.bed"
t_neg_features_file="$p_tmp_dir/negative_features.bed"
t_anti_features_file="$p_tmp_dir/anti-features.bed"
t_pos_anti_features_file="$p_tmp_dir/positive_anti-features.bed"
t_neg_anti_features_file="$p_tmp_dir/negative_anti-features.bed"


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

if [ -d $p_tmp_dir ]; then
	rm -r $p_tmp_dir
fi

for dir_path in $p_output_dir $p_tmp_dir $p_o_bed_dir $p_o_bam_dir $p_o_features_dir
do
	if [ ! -d $dir_path ]; then
		mkdir $dir_path
	fi
done

##              ##
# BEGIN ANALYSIS #
##              ##

echo "[INFO] converting arabic to roman chromosome names"
python scripts/arabic_to_roman.py $p_i_features_file $t_features_file
python scripts/arabic_to_roman.py $p_i_trna_locations_file $p_o_trna_locations_file


echo "[INFO] creating strand-wise feature files"
cat $t_features_file |
	grep '\+$' > $t_pos_features_file
cat $t_features_file |
	grep '\-$' > $t_neg_features_file

echo "[INFO] creating anti-features file '$t_anti_features_file'"
scripts/create_anti_features.sh $t_features_file $p_i_chr_lengths_file $t_anti_features_file

echo "[INFO] creating strand-wise anti-feature files"
cat $t_anti_features_file |
	grep '\+$' >> $t_pos_features_file
cat $t_anti_features_file |
	grep '\-$' >> $t_neg_features_file

echo "[INFO] sorting strand-wise feature files"
tmp_file_path="$p_tmp_dir/tmp"
for ff in $t_pos_features_file $t_neg_features_file
do
	mv $ff $tmp_file_path
	bedtools sort \
		-i $tmp_file_path \
		-g $p_i_chr_lengths_file >$ff
done
rm $tmp_file_path

echo "[INFO] creating the combined BED file of features and anti-features"
cat $t_pos_features_file $t_neg_features_file | \
	bedtools sort -g $p_i_chr_lengths_file > $p_o_features_file

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
		-a $t_pos_features_file \
		-b - \
		-g $p_i_chr_lengths_file \
		-sorted -counts \
		>> $bed_path

	# negative-features
	bedtools bamtobed -i $new_bam_path | \
		grep '\-$'  | \
		bedtools coverage \
		-a $t_neg_features_file \
		-b - \
		-g $p_i_chr_lengths_file \
		-sorted -counts \
		>> $bed_path

	# Write to the output file path
	echo "[INFO] sorting and concatenating BED file '$bed_path'"
	bedtools sort -g $p_i_chr_lengths_file -i $bed_path | \
		awk -v id="$uniqueId" '{ OFS="\t"; print id,$0 }' >> \
		$p_o_coverage_file
done<$p_i_sample_name_map


# TODO: summarise number of reads in BAMs and BEDs

echo "[INFO] Genome coverage information written to '$p_o_coverage_file'"
echo "[INFO] Workflow completed"
