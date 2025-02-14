#!/bin/bash

output_dir="data/output/feature_locations/"
file_prefix="Saccharomyces_cerevisiae_S288C"

feature_file_path="$output_dir/$file_prefix.gff"
wget \
	--output-document $feature_file_path.gz \
	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
gunzip $feature_file_path.gz
