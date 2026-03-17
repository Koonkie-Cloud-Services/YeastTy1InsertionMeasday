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

if [ "$#" -lt 1 ]; then
    log ERROR "unexpected number of arguments provided"
    echo -e "USAGE:\n$0 <working directory> [BAM suffix]"
    exit
fi

##                          ##
# DEFINE ANALYSIS PARAMETERS #
##                          ##
# Define directory paths
p_working_dir=$1
p_input_dir="$p_working_dir/input"
p_output_dir="$p_working_dir/output"
p_tmp_dir="$p_working_dir/.tmp"

# Analysis parameters
p_bam_suffix=$2
if [ -z "$p_bam_suffix" ]; then
	p_bam_suffix="dupremoved_uniq.sorted"
fi

p_i_features_file=$3
if [ ! -f "$p_i_features_file" ]; then
	log ERROR "input features file '${p_i_features_file}' doesn't exist"
	exit
fi

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
t_anti_features_file="$p_tmp_dir/anti-features.bed"
t_pos_anti_features_file="$p_tmp_dir/positive_anti-features.bed"
t_neg_anti_features_file="$p_tmp_dir/negative_anti-features.bed"


##                    ##
# PARAMETER VALIDATION #
##                    ##

if [ ! -d $p_working_dir ]; then
	log ERROR "working directory '$p_working_dir' doesn't exist"
	exit
fi

if [ ! -d $p_input_dir ]; then
	log ERROR "input directory '$p_input_dir' doesn't exist"
	exit
fi

if [ -f $p_o_coverage_file ]; then
	log WARNING "output file '$p_o_coverage_file' already exists. Returning existing file."
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

log INFO "converting arabic to roman chromosome names"
python scripts/arabic_to_roman.py $p_i_features_file $t_features_file


log INFO "creating strand-wise feature files"
cat $t_features_file |
	grep '\+$' > $t_pos_features_file
cat $t_features_file |
	grep '\-$' > $t_neg_features_file

log INFO "creating anti-features file '$t_anti_features_file'"
scripts/create_anti_features.sh $t_features_file $p_i_chr_lengths_file $t_anti_features_file

log INFO "creating strand-wise anti-feature files"
cat $t_anti_features_file |
	grep '\+$' >> $t_pos_features_file
cat $t_anti_features_file |
	grep '\-$' >> $t_neg_features_file

log INFO "sorting strand-wise feature files"
tmp_file_path="$p_tmp_dir/tmp"
for ff in $t_pos_features_file $t_neg_features_file
do
	mv $ff $tmp_file_path
	bedtools sort \
		-i $tmp_file_path \
		-g $p_i_chr_lengths_file >$ff
	if [ $? -ne 0 ]; then
		log ERROR "Command failed while processing file '$tmp_file_path'"
		exit 1
	fi
done
rm $tmp_file_path

log INFO "creating the combined BED file of features and anti-features"
cat $t_pos_features_file $t_neg_features_file | \
	bedtools sort -g $p_i_chr_lengths_file > $p_o_features_file

while IFS=$'\t' read -r sampleId strainId sugar uniqueId runId barcode libraryName fastqName || [ -n "$fastqName" ]
do
	# Skip the header line
	if [[ $sampleId == "sampleId" ]]; then
		continue
	fi

	# Define file paths for sample
	old_bam="${p_i_bam_dir}/${runId}_${barcode}_${p_bam_suffix}.bam"
	new_bam_path="${p_o_bam_dir}/${uniqueId}_${p_bam_suffix}.bam"
	bed_path="${p_o_bed_dir}/${uniqueId}_${p_bam_suffix}_feature_coverage.bed"

	if [ ! -f $old_bam ]; then
		log ERROR "Unable to find input BAM file '$old_bam'"
		exit
	fi

	log INFO "Renaming $uniqueId BAM: '$old_bam'"
	cp $old_bam $new_bam_path

	if [ -f $bed_path ]; then
		rm $bed_path
	fi
	touch $bed_path

	log INFO "calculating feature coverage for sample: $uniqueId"
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
	log INFO "sorting and concatenating BED file '$bed_path'"
	bedtools sort -g $p_i_chr_lengths_file -i $bed_path | \
		awk -v id="$uniqueId" '{ OFS="\t"; print id,$0 }' >> \
		$p_o_coverage_file
done<$p_i_sample_name_map

if [ ! -s "$p_o_coverage_file" ]; then
	log ERROR "Output file '$p_o_coverage_file' is empty or does not exist."
	exit 1
fi

# TODO: summarise number of reads in BAMs and BEDs

log INFO "Genome coverage information written to '$p_o_coverage_file'"
log INFO "Workflow completed"
