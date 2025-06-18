#!/bin/bash

##
#  FUNCTION DEFINITIONS
##
log_message() {
  local level="$1"
  shift
  local message="$@"
  case "$level" in
    info)
      echo "INFO: $message" >&2
      ;;
    warning)
      echo "WARNING: $message" >&2
      ;;
    error)
      echo "ERROR: $message" >&2
      ;;
    debug)
      echo "DEBUG: $message" >>$LOG_FILE_PATH
      ;;
    *)
      echo "INVALID LEVEL: $message"
      ;;
  esac
}

run_command() {
    local cmd="$*"
    log_message debug "Running command: $cmd"
    eval "$cmd"
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        log_message error "Command '$cmd' failed with exit code $exit_code" >&2
        exit $exit_code
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
		 log_message error "No .fastq.gz files found in ${raw_fastq_dir}."
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

	log_message info $'\tRemoving AACA sequence'
	run_command cutadapt -g "^AACA" \
		-o "${working_dir}/${sample_id}.5primetrimmed.fastq" \
		--discard-untrimmed \
		"${fastq_r1}" \
		1> "${working_dir}/${sample_id}.ssbstats.txt"

	log_message info $'\tRemoving the adapter'
	run_command cutadapt -a "TAGTCCCTTAAGCGGAGCCCTATAG" \
		-o "${working_dir}/${sample_id}.5primetrimmed.3primetrimmed.fastq" \
		-m 20 \
		"${working_dir}/${sample_id}.5primetrimmed.fastq" \
		> "${working_dir}/${sample_id}.adapterstats.txt"

	log_message info $'\tBWA mapping'
	run_command bwa mem "$bwa_index" -M -t "$nslots" \
		-R "'@RG\tID:${sample_id}\tSM:${sample_id}\tPL:illumina'" \
		"${working_dir}/${sample_id}.5primetrimmed.3primetrimmed.fastq" \
		-o "${working_dir}/${sample_id}.sam"

	log_message info $'\tSamtools sorting BAM'
	run_command samtools view -F 0x900 -@ "$nslots" -ubS "${working_dir}/${sample_id}.sam" | \
		samtools sort - -@ "$nslots" -T "${working_dir}/${sample_id}.sorted" \
		-O BAM -o "${working_dir}/${sample_id}.sorted.bam"

	log_message info $'\tIndexing BAM'
	run_command samtools index "${working_dir}/${sample_id}.sorted.bam"
	run_command samtools flagstat "${working_dir}/${sample_id}.sorted.bam" > "${working_dir}/${sample_id}.sorted.bam.flagstat"

	log_message info $'\tRemoving duplicates to make a mapped.dedup.bam'
	local tmpdir="${working_dir}/picardtmp/"
	rm -rf "$tmpdir"
	mkdir "$tmpdir"
	run_command picard MarkDuplicates \
		I="${working_dir}/${sample_id}.sorted.bam" \
		REMOVE_DUPLICATES=true TMP_DIR="$tmpdir" \
		M="${working_dir}/${sample_id}_uniq_DuplicateResults.txt" \
		O="${working_dir}/${sample_id}_dupremoved.sorted.bam"

	log_message info $'\tIndexing deduplicated BAM'
	run_command samtools index "${working_dir}/${sample_id}_dupremoved.sorted.bam"
	run_command samtools flagstat "${working_dir}/${sample_id}_dupremoved.sorted.bam" > "${working_dir}/${sample_id}_dupremoved.sorted.bam.flagstat"

	log_message info $'\tRemoving non-unique mapping reads'
	run_command samtools view -q 30 "${working_dir}/${sample_id}.sorted.bam" -o "${working_dir}/${sample_id}_uniq.sorted.bam"
	run_command samtools index "${working_dir}/${sample_id}_uniq.sorted.bam"
	run_command samtools flagstat "${working_dir}/${sample_id}_uniq.sorted.bam" > "${working_dir}/${sample_id}_uniq.sorted.bam.flagstat"

	log_message info $'\tPicard MarkDuplicates for unique reads'
	rm -rf "$tmpdir"
	mkdir "$tmpdir"
	run_command picard MarkDuplicates \
		I="${working_dir}/${sample_id}_uniq.sorted.bam" \
		REMOVE_DUPLICATES=true TMP_DIR="$tmpdir" \
		M="${working_dir}/${sample_id}_uniq_DuplicateResults.txt" \
		O="${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam"

	run_command samtools index "${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam"
	run_command samtools flagstat "${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam" > "${working_dir}/${sample_id}_dupremoved_uniq.sorted.bam.flagstat"

	log_message info $'\tCleaning up intermediate files'
	rm -f "${working_dir}/${sample_id}.sam"
	rm -f "${working_dir}/${sample_id}.5primetrimmed.fastq"
	rm -f "${working_dir}/${sample_id}.5primetrimmed.3primetrimmed.fastq"
	rm -rf "$tmpdir"

	echo "Sample $sample_id processing completed at $(date)"
	echo "----------------------------------------"
	echo -e "\n"
}
##
# END FUNCTION DEFINITIONS
##

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
	echo -e "Usage: $0 <working_directory>"
	
	echo """
Description:
	This script aligns reads to a reference genome and performs various processing steps.

	The working directory should contain the following structure:
		<working_directory>/data/input/FASTQ - Contains .fastq.gz files
		<working_directory>/data/input/genome - Contains the reference genome
		<working_directory>/data/output - Output directory for results

	The script will create the following directories:
		<working_directory>/data/output/alignments - For alignment results
		<working_directory>/data/output/genome - For the reference genome (if none provided)
	"""
	exit 0
fi


WORKING_DIR=$1
N_CORE=$2
LOG_FILE_PATH="${WORKING_DIR}/align_reads_to_reference_genome.log"

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
BWA_INDEX="${GENOME_DIR}/reference_genome"

# Check if required executables are available
for tool in cutadapt bwa samtools; do
	command -v $tool >/dev/null 2>&1 || { log_message error "$tool is not installed or not in PATH."; exit 1; }
done

for dir_path in "$OUTPUT_DIR" "$ALIGN_DIR"; do
	if [[ ! -d "$dir_path" ]]; then
		log_message info "Creating directory: $dir_path"
		mkdir -p "$dir_path"
	fi
done

log_message info "Primary Analysis Started"
date

SAMPLE_NAMES=($(find_sample_names "${RAW_FASTQ_DIR}"))
if [[ ${#SAMPLE_NAMES[@]} -eq 0 ]]; then
	log_message error "No sample names found. Exiting."
	exit 1
fi
log_message info "Number of samples found: ${#SAMPLE_NAMES[@]}"
log_message info "Found sample names:"
for sample in "${SAMPLE_NAMES[@]}"; do
	echo "$sample"
done

# Check if BWA index file exists
if [[ ! -f "$BWA_INDEX.sa" ]]; then
	GENOME_DIR="$OUTPUT_DIR/genome"
	# Update BWA_INDEX variable to point to the newly created index
	BWA_INDEX="${GENOME_DIR}/reference_genome"
	mkdir -p "$GENOME_DIR"

	log_message info "Downloading reference genome..."
	run_command curl -o "${GENOME_DIR}/reference_genome.fna.gz" "$REF_GENOME_URL"

	log_message info "Building BWA index..."
	run_command bwa index -p "${BWA_INDEX}" "${GENOME_DIR}/reference_genome.fna.gz"
fi

for sample_id in "${SAMPLE_NAMES[@]}"; do
	log_message info "Processing sample: $sample_id"
	# Skip processing if the output file already exists
	output_file="${ALIGN_DIR}/${sample_id}_dupremoved_uniq.sorted.bam"
	if [[ -f "$output_file" ]]; then
		log_message info "Output file for sample $sample_id already exists. Skipping processing."
		continue
	fi
	process_sample "$sample_id" "$ALIGN_DIR" "$RAW_FASTQ_DIR" "$BWA_INDEX" $N_CORE
done

log_message info """
Primary Analysis Completed at $(date)
----------------------------------------
"""