#!/bin/bash

output_dir=$1

if [ ! -d $output_dir ]; then
	echo "[ERROR] output directory '$output_dir' doesn't exist"
	exit
fi

file_prefix="Saccharomyces_cerevisiae_S288C"
ncbi_ftp_path="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic"
sgd_ftp_path="http://sgd-archive.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source"

feature_file_path="$output_dir/$file_prefix.gff"
genbank_file_path="$output_dir/$file_prefix.gbff"

if [ ! -f $feature_file_path ]; then
	wget \
	--output-document $feature_file_path.gz \
	${ncbi_ftp_path}.gff.gz

	# Decompress
	gunzip $feature_file_path.gz
fi

# if [ ! -f $genbank_file_path ]; then
# 	wget \
# 	--output-document $genbank_file_path.gz \
# 	${ncbi_ftp_path}.gbff.gz
# 
# 	# Decompress
# 	gunzip $genbank_file_path.gz
# fi

# for i in $(seq -f "%02g" 1 16); do
#   echo $i
#   curl \
#     ${sgd_ftp_path}/chr$i.gbf > \
#     "$output_dir/${file_prefix}_chr${i}.gbff"
# done
