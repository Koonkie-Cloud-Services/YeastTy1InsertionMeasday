#!/bin/bash

# Wrapper function to call a command and exit if it fails
run_command() {
    local cmd="$*"
    echo "Running command: $cmd"
    eval "$cmd"
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo "Error: Command '$cmd' failed with exit code $exit_code" >&2
        exit $exit_code
    fi
}

working_dir="data"
N_CPUS=4

if [ -z $( which bedtools) ]; then
	echo "[ERROR] bedtools executable not in path"
	exit
fi

if [ ! -d $input_dir ]; then
	echo "[ERROR] input directory '$input_dir' doesn't exist"
	exit 1
fi

input_dir="$working_dir/input"
output_dir="$working_dir/output"
p_o_bam_dir="${output_dir}/alignments"
p_o_aln_dir="${output_dir}/merged_alignments"
p_o_aln_flt_dir="${output_dir}/filtered_alignments"
p_o_coverage_dir="${output_dir}/genome_coverage"
p_o_features_dir="${output_dir}/feature_locations"
p_i_trna_locations_file="$input_dir/trna_regions.edited.bed"
p_i_filter_regions_file="$input_dir/filter_regions.bed"
p_o_trna_locations_file="$p_o_features_dir/trna_regions.bed"

if [ ! -d $output_dir ]; then
	mkdir $output_dir
fi

for dir_path in $p_o_features_dir $p_o_coverage_dir
do
	if [ ! -d $dir_path ]; then
		mkdir $dir_path
	fi
done

# The following command writes BAM files to $output_dir/alignments
./scripts/align_reads_to_reference_genome.sh $working_dir $N_CPUS
./scripts/merge_bam_files.sh $p_o_bam_dir $input_dir/BAMs
./scripts/filter_bam_files.sh $p_o_aln_dir $p_o_aln_flt_dir $p_i_filter_regions_file
ln -s "../../$p_o_aln_flt_dir" "$input_dir/alignments"
if [ ! -d "$input_dir/alignments" ] || [ -z "$(ls -A $input_dir/alignments/*.bam 2>/dev/null)" ]; then
	echo "[ERROR] No BAM files found in '$input_dir/alignments'."
	exit 1
fi

# Download genome annotation files from NCBI
./scripts/download_reference_genome_data.sh $p_o_features_dir

python scripts/arabic_to_roman.py $p_i_trna_locations_file $p_o_trna_locations_file

# Convert the GFF features to BED, add upstream+downstream as needed
run_command python scripts/reformat_genome_features_to_bed.py \
	--gff_in "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.gff" \
	--bed_out "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.bed" \
	--buffer_features "$input_dir/buffer_regions.tsv"

run_command python scripts/subset_features_for_polymerase_analysis.py \
	--gff_in "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.gff" \
	--xlsx_in "${input_dir}/custom_feature_annotations.xlsx" \
	--bed_out_prefix "${p_o_features_dir}/polymerase_analysis" \
	--buffer_length_bp 1000 \
	--adjacent_buffer_length 1500 \
	--sheet_name "polymerase_analysis" \
	--drop_unmatched \
	--nucleosome_data "${input_dir}/13059_2018_1398_MOESM2_ESM.xlsx"

# Calculate genome coverage for all features
ln -s $PWD/$input_dir $p_o_coverage_dir
ln -s $PWD/$p_o_features_dir $p_o_coverage_dir/$( basename $input_dir)
./scripts/calculate_feature_coverage.sh \
	$p_o_coverage_dir "dupremoved.sorted" "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.bed"
./scripts/calculate_feature_coverage.sh \
	$p_o_coverage_dir "dupremoved_uniq.sorted" "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.bed"
# Calculate coverage of features used in polymerase analysis
./scripts/calculate_feature_coverage.sh \
	$p_o_coverage_dir "dupremoved.sorted" "${p_o_features_dir}/polymerase_analysis.RNA_Pol_3.bed"
./scripts/calculate_position_coverage.sh \
	$p_o_coverage_dir "${p_o_features_dir}/polymerase_analysis.RNA_Pol_2.bed" "dupremoved.sorted"
./scripts/count_rnda_locus_reads.sh \
	$p_o_aln_dir "${p_o_coverage_dir}/rDNA_locus_counts.tsv"
unlink $p_o_coverage_dir/$( basename $input_dir)

# Render the RMarkdown report to create figures
Rscript -e "rmarkdown::render('./Ty1_insertion_report.Rmd')"
