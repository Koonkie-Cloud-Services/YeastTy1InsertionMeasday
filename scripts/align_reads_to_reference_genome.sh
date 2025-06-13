#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
	echo "Usage: $0 <working_directory>"
	echo "This script aligns reads to a reference genome and performs various processing steps."
	exit 0
fi

##
#  FUNCTION DEFINITIONS
##
run() {
	local cmd="$@"
	echo "Running: $cmd"
	eval "$cmd"
	if [[ $? -ne 0 ]]; then
		echo "Error: Command failed - $cmd"
		exit 1
	fi
}

find_sample_names() {
	# Function to find all fastq.gz files and extract sample names
	local raw_fastq_dir=$1
	local sample_names=()

	for file in "${raw_fastq_dir}"/*.fastq.gz; do
		if [[ -f "$file" ]]; then
			local sample_name=$(basename "$file" .fastq.gz | sed 's/_100bp.concat_chastity_passed//g')
			sample_names+=("$sample_name")
		fi
	done

	# Check if any sample names were found
	if [[ ${#sample_names[@]} -eq 0 ]]; then
		echo "Error: No .fastq.gz files found in ${raw_fastq_dir}."
		exit 1
	fi

	echo "${sample_names[@]}"
}

process_sample() {
	local sample_id=$1
	local working_dir=$2
	local raw_fastq_dir=$3
	local bwa_index=$4
	local nslots=$5

	fastq_r1=$( find "${raw_fastq_dir}" -name "${sample_id}*.fastq.gz")

	if [[ ! -f "${fastq_r1}" ]]; then
		echo "Error: Fastq file ${fastq_r1} not found in ${raw_fastq_dir}."
		exit
	fi

	echo $'\tRemoving AACA sequence'
	run cutadapt -g "^AACA" \
		-o "${working_dir}/${sample_id}.5primetrimmed.fastq" \
		--discard-untrimmed \
		"${fastq_r1}" \
		> "${working_dir}/${sample_id}.ssbstats.txt"

	echo $'\tRemoving the adapter'
	run cutadapt -a "TAGTCCCTTAAGCGGAGCCCTATAG" \
		-o "${working_dir}/${sample_id}.5primetrimmed.3primetrimmed.fastq" \
		-m 20 \
		"${working_dir}/${sample_id}.5primetrimmed.fastq" \
		> "${working_dir}/${sample_id}.adapterstats.txt"

	echo $'\tBWA mapping'
	run bwa mem "$bwa_index" -M -t "$nslots" -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:illumina" \
		"${working_dir}/${sample_id}.5primetrimmed.3primetrimmed.fastq" > "${working_dir}/${sample_id}.sam"

	echo $'\tSamtools sorting BAM'
	run samtools view -F 0x900 -@ "$nslots" -ubS "${working_dir}/${sample_id}.sam" | \
		samtools sort - -@ "$nslots" -T "${working_dir}/${sample_id}.sorted" -O BAM -o "${working_dir}/${sample_id}.sorted.bam"

	echo $'\tIndexing BAM'
	run samtools index "${working_dir}/${sample_id}.sorted.bam"
	run samtools flagstat "${working_dir}/${sample_id}.sorted.bam" > "${working_dir}/${sample_id}.sorted.bam.flagstat"

	echo $'\tRemoving duplicates to make a mapped.dedup.bam'
	local tmpdir="${working_dir}/picardtmp/"
	rm -rf "$tmpdir"
	mkdir "$tmpdir"
	run java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
		I="${working_dir}/${sample_id}.sorted.bam" \
		REMOVE_DUPLICATES=true TMP_DIR="$tmpdir" \
		M="${working_dir}/${sample_id}_uniq_DuplicateResults.txt" \
		O="${working_dir}/${sample_id}_dupremoved.sorted.bam"

	echo $'\tIndexing deduplicated BAM'
	run samtools index "${working_dir}/${sample_id}_dupremoved.sorted.bam"
	run samtools flagstat "${working_dir}/${sample_id}_dupremoved.sorted.bam" > "${working_dir}/${sample_id}_dupremoved.sorted.bam.flagstat"

	echo $'\tRemoving non-unique mapping reads'
	run samtools view -q 30 "${working_dir}/${sample_id}.sorted.bam" -o "${working_dir}/${sample_id}_uniq.sorted.bam"
	run samtools index "${working_dir}/${sample_id}_uniq.sorted.bam"
	run samtools flagstat "${working_dir}/${sample_id}_uniq.sorted.bam" > "${working_dir}/${sample_id}_uniq.sorted.bam.flagstat"

	echo $'\tPicard MarkDuplicates for unique reads'
	rm -rf "$tmpdir"
	mkdir "$tmpdir"
	run java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
		I="${working_dir}/${sample_id}_uniq.sorted.bam" \
		REMOVE_DUPLICATES=true TMP_DIR="$tmpdir" \
		M="${working_dir}/${sample_id}_uniq_DuplicateResults.txt" \
		O="${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam"

	run samtools index "${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam"
	run samtools flagstat "${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam" > "${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam.flagstat"

	echo "Sample $sample_id processing completed at $(date)"
	echo "----------------------------------------"
	echo -e "\n"
}
##
# END FUNCTION DEFINITIONS
##

WORKING_DIR=$1
N_CORE=$2

if [ -z "$WORKING_DIR" ]; then
	echo "Usage: $0 <working_directory> [<number_of_cores>]"
	exit 1
fi

if [ -z "$N_CORE" ]; then
	N_CORE=$(nproc)  # Default number of cores if not specified
fi

# Inputs
INPUT_DIR="${WORKING_DIR}/data/input"
RAW_FASTQ_DIR="${INPUT_DIR}/FASTQ"
GENOME_DIR="${INPUT_DIR}/genome"
# Outputs
OUTPUT_DIR="${WORKING_DIR}/data/output"
ALIGN_DIR="${OUTPUT_DIR}/alignments"
# Hard-coded parameters
REF_GENOME_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
BWA_INDEX="${GENOME_DIR}/S288C_reference_sequence_R64-2-1_20150113.fsa"

# Check if required executables are available
for tool in cutadapt bwa samtools; do
	command -v $tool >/dev/null 2>&1 || { echo "Error: $tool is not installed or not in PATH."; exit 1; }
done

for dir_path in "$OUTPUT_DIR" "$ALIGN_DIR"; do
	if [[ ! -d "$dir_path" ]]; then
		echo "Creating directory: $dir_path"
		mkdir -p "$dir_path"
	fi
done

echo "Primary Analysis Started"
date

SAMPLE_NAMES=($(find_sample_names "${RAW_FASTQ_DIR}"))
echo "Number of samples found: ${#SAMPLE_NAMES[@]}"
echo "Found sample names:"
for sample in "${SAMPLE_NAMES[@]}"; do
	echo "$sample"
done

# Check if BWA index file exists
if [[ ! -f "$BWA_INDEX" ]]; then
	GENOME_DIR="$OUTPUT_DIR/genome"
	# Update BWA_INDEX variable to point to the newly created index
	BWA_INDEX="${GENOME_DIR}/reference_genome"
	mkdir -p "$GENOME_DIR"

	echo "Downloading reference genome..."
	curl -o "${GENOME_DIR}/reference_genome.fna.gz" "$REF_GENOME_URL"

	echo "Building BWA index..."
	run bwa index -p "${BWA_INDEX}" "${GENOME_DIR}/reference_genome.fna.gz"
fi

for sample_id in "${SAMPLE_NAMES[@]}"; do
	echo "Processing sample: $sample_id"
	process_sample "$sample_id" "$ALIGN_DIR" "$RAW_FASTQ_DIR" "$BWA_INDEX" $N_CORE
done

echo "Primary Analysis Completed at $(date)"
echo "----------------------------------------"
echo $'\n'