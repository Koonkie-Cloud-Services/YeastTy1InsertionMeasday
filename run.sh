#!/bin/bash

./scripts/download_reference_genome_data.sh
python ./scripts/reformat_genbank_features_to_bed.py
./scripts/generate_genome_coverage.sh

Rscript -e "rmarkdown::render('./Ty1_insertion_report.Rmd')"
