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

if [ -z $( which bedtools) ]; then
	echo "[ERROR] bedtools executable not in path"
	exit
fi

input_dir="$working_dir/input"
output_dir="$working_dir/output"
p_o_coverage_dir="${output_dir}/genome_coverage"
p_o_features_dir="${output_dir}/feature_locations"
p_i_trna_locations_file="$input_dir/trna_regions.edited.bed"
p_o_trna_locations_file="$p_o_features_dir/trna_regions.bed"

if [ -d $output_dir ]; then
	rm -r $output_dir
fi
mkdir $output_dir

for dir_path in $p_o_features_dir $p_o_coverage_dir
do
	mkdir $dir_path
done

# Download input data from NCBI
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
./scripts/calculate_position_coverage.sh \
	$p_o_coverage_dir "${p_o_features_dir}/polymerase_analysis.RNA_Pol_2.bed" "dupremoved.sorted"
./scripts/calculate_feature_coverage.sh \
	$p_o_coverage_dir "dupremoved.sorted" "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.bed"
./scripts/calculate_feature_coverage.sh \
	$p_o_coverage_dir "dupremoved_uniq.sorted" "$p_o_features_dir/Saccharomyces_cerevisiae_S288C.bed"
# Calculate coverage of features used in polymerase analysis
./scripts/calculate_feature_coverage.sh \
	$p_o_coverage_dir "dupremoved.sorted" "${p_o_features_dir}/polymerase_analysis.RNA_Pol_3.bed"
unlink $p_o_coverage_dir/$( basename $input_dir)

# Render the RMarkdown report to create figures
Rscript -e "rmarkdown::render('./Ty1_insertion_report.Rmd')"
