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
      echo "DEBUG: $message" >>$P_I_LOG_FILE_PATH
      ;;
    *)
      echo "INVALID LEVEL: $message"
      ;;
  esac
}

# Inputs
P_I_BAM_DIRS="$@"
P_I_INPUT_DIR="data/input"
P_I_SAMPLE_ID_MAP="$P_I_INPUT_DIR/hiseq_sample_id_map.tsv"
P_I_LOG_FILE_PATH="./merge_bam_files.log"
# Outputs
P_O_WORKING_DIR="data"
P_O_OUTPUT_DIR="${P_O_WORKING_DIR}/output"
P_O_ALIGN_DIR="${P_O_OUTPUT_DIR}/merged_alignments"

if [ ! -d "$P_O_ALIGN_DIR" ]; then
    mkdir -p "$P_O_ALIGN_DIR"
fi

# Combine BAM files from the same sample
while IFS=$'\t' read -r sampleId strainId sugar uniqueId runId barcode libraryName fastqName || [ -n "$fastqName" ]
do
	# Skip the header line
	if [[ $sampleId == "sampleId" ]]; then
		continue
	fi

    # Get all BAM files for the sample across multiple directories
    sample_bam_files=()
    for bam_dir in $P_I_BAM_DIRS; do
        for suffix in "_dupremoved.sorted" "_dupremoved_uniq.sorted"; do
            sample_bam_files+=("${bam_dir}/${runId}_*${barcode}*${suffix}.bam")
        done
    done

	if [[ ${#sample_bam_files[@]} -gt 0 ]]; then
		log_message info "Combining BAM files for sample: $sampleId"

		for suffix in "_dupremoved.sorted" "_dupremoved_uniq.sorted"; do
            # Filter BAM files by suffix
            filtered_bam_files=()
            for bam_file in "${sample_bam_files[@]}"; do
                if [[ $bam_file == *"${suffix}.bam" ]]; then
                    filtered_bam_files+=("$bam_file")
                fi
            done

			new_bam_name="${runId}_${barcode}${suffix}.bam"
			new_bam_filepath="${P_O_ALIGN_DIR}/${new_bam_name}"
			if [ -f "${new_bam_filepath}" ]; then
				log_message info "Merged BAM file already exists: ${new_bam_filepath}"
				continue
			fi

			samtools merge \
                -o - \
                ${filtered_bam_files[@]} | \
            samtools sort \
                -o "${new_bam_filepath}" \
                -T "${P_O_ALIGN_DIR}/.tmp" \
                --threads "$(nproc)" \
                -
			log_message info "Combined BAM file created: ${new_bam_filepath}"
			done
	else
		log_message warning "No BAM files found for sample: $sampleId"
	fi
done < "$P_I_SAMPLE_ID_MAP"