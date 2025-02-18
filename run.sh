#!/bin/bash

working_dir="data"

if [ -z $( which bedtools) ]; then
	echo "[ERROR] bedtools executable not in path"
	exit
fi

input_dir="$working_dir/input"
output_dir="$working_dir/output"
genome_cov_dir="${working_dir}/output/genome_coverage"
feature_loc_dir="${working_dir}/output/feature_locations"

if [ -d $output_dir ]; then
	rm -r $output_dir
fi
mkdir $output_dir

for dir_path in $feature_loc_dir $genome_cov_dir
do
	mkdir $dir_path
done

# Download input data from NCBI
./scripts/download_reference_genome_data.sh $feature_loc_dir

# Convert the GFF features to BED, add upstream+downstream as needed
python scripts/reformat_genome_features_to_bed.py \
	--gff_in "$feature_loc_dir/Saccharomyces_cerevisiae_S288C.gff" \
	--bed_out "$feature_loc_dir/Saccharomyces_cerevisiae_S288C.bed" \
	--buffer_features "$input_dir/buffer_regions.tsv"

# Calculate genome coverage for features
ln -s $PWD/$input_dir $genome_cov_dir
ln -s $PWD/$feature_loc_dir $genome_cov_dir/$( basename $input_dir)
./scripts/generate_genome_coverage.sh $genome_cov_dir "dupremoved.sorted"
./scripts/generate_genome_coverage.sh $genome_cov_dir "dupremoved_uniq.sorted"
unlink $genome_cov_dir/$( basename $input_dir)

# Render the RMarkdown report to create figures
Rscript -e "rmarkdown::render('./Ty1_insertion_report.Rmd')"
