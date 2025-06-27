#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "ERROR: unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <working directory> <features BED> [BAM suffix] [input directory]"
    exit
fi

##                          ##
# DEFINE ANALYSIS PARAMETERS #
##                          ##
# Define directory paths
p_working_dir=$1
p_output_dir="$p_working_dir/output"
p_tmp_dir="$p_working_dir/.tmp"

# Analysis parameters
p_i_features_file=$2
if [ ! -f "$p_i_features_file" ]; then
	echo "ERROR: input features file '${p_i_features_file}' doesn't exist"
	exit
fi

p_bam_suffix=$3
if [ -z "$p_bam_suffix" ]; then
	p_bam_suffix="dupremoved_uniq.sorted"
fi

if [ -z $4 ]; then
	p_input_dir="$p_working_dir/input"
else
	p_input_dir=$4
fi
FLOWCELL="2"

# Input paths
p_i_bam_dir="$p_input_dir/alignments"
p_i_features_dir="$p_input_dir/feature_locations"
p_i_sample_name_map="$p_input_dir/hiseq_sample_id_map.tsv"
p_i_chr_lengths_file="$p_input_dir/chromosome_lengths.tsv"
p_features_prefix=$( basename $p_i_features_file | sed 's/.bed$//g' )

# Output paths
p_o_bam_dir="$p_output_dir/BAMs"
p_o_bed_dir="$p_output_dir/BEDs"
p_o_features_dir="$p_output_dir/feature_locations"
p_o_features_file="$p_o_features_dir/${p_features_prefix}_features.bed"
p_o_coverage_file="$p_output_dir/${p_features_prefix}_coverage_${p_bam_suffix}.bed"

# Intermediate file paths
t_features_file="$p_tmp_dir/features.bed"
t_pos_features_file="$p_tmp_dir/positive_features.bed"
t_neg_features_file="$p_tmp_dir/negative_features.bed"


##                    ##
# PARAMETER VALIDATION #
##                    ##

if [ ! -d $p_working_dir ]; then
	echo "[ERROR] working directory '$p_working_dir' doesn't exist"
	exit
fi

if [ ! -d $p_input_dir ]; then
	echo "[ERROR] input directory '$p_input_dir' doesn't exist"
	exit
fi

if [ -f $p_o_coverage_file ]; then
	echo "[WARNING] output file '$p_o_coverage_file' already exists. Returning existing file."
	exit 0
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


echo "[INFO] creating strand-wise feature files from '$t_features_file'"
cat $t_features_file |
	grep '\+$' > $t_pos_features_file
cat $t_features_file |
	grep '\-$' > $t_neg_features_file

echo "[INFO] sorting strand-wise feature files"
tmp_file_path="$p_tmp_dir/tmp"
for ff in $t_pos_features_file $t_neg_features_file
do
	mv $ff $tmp_file_path
	bedtools sort \
		-i $tmp_file_path \
		-g $p_i_chr_lengths_file >$ff
	if [ $? -ne 0 ]; then
		echo "[ERROR] unable to sort file '$ff'"
		exit
	fi
	# Remove the temporary file
	rm $tmp_file_path
done

echo "[INFO] creating the combined BED file of features"
cat $t_pos_features_file $t_neg_features_file | \
	bedtools sort -g $p_i_chr_lengths_file > $p_o_features_file

while IFS=$'\t' read -r sampleId strainId sugar uniqueId runId barcode libraryName fastqName || [ -n "$fastqName" ]
do
	# Skip the header line
	if [[ $sampleId == "sampleId" ]]; then
		continue
	fi

	# Define file paths for sample
	old_bam="${p_i_bam_dir}/${runId}_${FLOWCELL}_1_${barcode}_${p_bam_suffix}.bam"
	new_bam_path="${p_o_bam_dir}/${uniqueId}_${p_bam_suffix}.bam"
	tmp_bed_path="${p_tmp_dir}/${uniqueId}_${p_bam_suffix}_feature_coverage.bed"
	bed_path="${p_o_bed_dir}/${uniqueId}_${p_bam_suffix}_feature_coverage.bed"

	if [ ! -f $old_bam ]; then
		echo "[ERROR] Unable to find input BAM file '$old_bam'"
		exit
	fi

	echo -e "[INFO] Renaming $uniqueId BAM:\nsrc:\t'$old_bam'\ndst:\t'$new_bam_path'"
	cp $old_bam $new_bam_path

	for file_path in $tmp_bed_path $bed_path
	do
		if [ -f $file_path ]; then
			rm $file_path
		fi
		touch $file_path
	done

	echo "[INFO] calculating feature coverage for sample: $uniqueId"
	bedtools genomecov \
		-ibam $new_bam_path \
		-d \
		-strand + \
		-fs 1 | \
		awk -v strand="+" -F"\t" '{OFS="\t"; print $1,$2,$3,strand }' \
		>> $tmp_bed_path

	bedtools genomecov \
		-ibam $new_bam_path \
		-d \
		-strand - \
		-fs 1 | \
		awk -v strand="-" -F"\t" '{OFS="\t"; print $1,$2,$3,strand }' \
		>> $tmp_bed_path

	# Intersect the coverage with the features
	cat $tmp_bed_path | \
	awk -F"\t" '{ OFS="\t"; print $1,$2,$2,$3,$4 }' | \
	bedtools intersect \
	-a $p_o_features_file -b - \
	-wa -wb | \
	awk -F"\t" '{ OFS="\t"; if ( $5 == $10 ) print $1,$7,$5,$4,$9 }' > $bed_path

	# Write to the output file path
	echo "[INFO] concatenating BED file '$bed_path'"
	cat $bed_path | \
		awk -F"\t" -v id="$uniqueId" '{ OFS="\t"; if ($NF > 0) print id,$0 }' >> \
		$p_o_coverage_file

done<$p_i_sample_name_map

if [ ! -s "$p_o_coverage_file" ]; then
	echo "[ERROR] Output file '$p_o_coverage_file' is empty or does not exist."
	exit 1
fi

# TODO: summarise number of reads in BAMs and BEDs

echo "[INFO] Genome coverage information written to '$p_o_coverage_file'"
echo "[INFO] Workflow completed"
