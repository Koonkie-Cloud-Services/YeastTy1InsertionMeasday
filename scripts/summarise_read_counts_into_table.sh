#!/bin/bash

# This bash script processes sequencing data and summarizes read counts at different processing stages. It:

# 1. **Reads a sample map file** containing library metadata (barcode, sample info)
# 2. **Constructs library names** by combining sample fields with flowcell number
# 3. **Finds alignment files** for each barcode in the output directory
# 4. **Counts reads** at multiple stages:
#     - Mapped reads (BAM file)
#     - Deduplicated reads
#     - Unique/high-quality reads
# 5. **Outputs a TSV table** with one row per sample showing read counts through the pipeline

# It serves as a quality control summary tracking how many reads survive each processing step (trimming, alignment, deduplication, filtering).

SAMPLE_MAP_FILE="data/input/hiseq_sample_id_map.tsv"
OUTPUT_FILE="data/output/summarised_read_counts.tsv"

FASTQ_DIR="data/input/FASTQ"
OUTPUT_DIR="data/output"
FC_ONE_ALN_DIR="data/input/BAMs"
FC_TWO_ALN_DIR="$OUTPUT_DIR/alignments"

if [ -f "$OUTPUT_FILE" ]; then
    rm "$OUTPUT_FILE"
fi
touch "$OUTPUT_FILE"
echo -e "Unique_ID\tFlowcell\tMapped_Reads\tDeduplicated_Reads\tUnique_Reads" >> "$OUTPUT_FILE"

declare -a FLOWCELL_NUMS=("1" "2")
while read line; do
    for flowcell in "${FLOWCELL_NUMS[@]}"; do
        # Extract the sample name and read count from each line
        library_name=$(echo "$line" | awk -F"\t" '{OFS="_"; print $5,"'$flowcell'","1",$6}')
        barcode=$(echo "$line" | awk -F"\t" '{print $6}')
        unique_id=$(echo "$line" | awk -F"\t" '{print $4}')

        # Check if the library name is empty
        if [ -z "$library_name" ]; then
            echo "Skipping empty library name in line: $line"
            continue
        elif [[ "$barcode" == "barcode" ]]; then
            continue
        fi

        >&2 echo "Processing library: $library_name with barcode: $barcode"

        if [ "$flowcell" == "1" ]; then
            alignment_dir="$FC_ONE_ALN_DIR"
        else
            alignment_dir="$FC_TWO_ALN_DIR"
        fi
        alignment_files="$(find "$alignment_dir"/ -name "*${barcode}*")"

        if [ -z "$alignment_files" ]; then
            >&2 echo "No alignment files found for barcode: $barcode in directory: $alignment_dir"
            exit 1
        fi

        # SKIPPED: count reads from FASTQ files - missing files from flowcell 1
        # in_read_count=$(find "$FASTQ_DIR" -name "${library_name}_*.fastq.gz" -exec zcat {} + | wc -l | awk '{ print $1/4 }')

        # Count reads from trimming stats files
        primer_read_count="$(find "$alignment_dir" -name "*${barcode}.ssbstats.txt" -exec grep '^Reads written' {} + | awk -F: '{ print $2 }' | awk '{ print $1 }' | sed 's/,//g')"
        adapter_read_count="$(find "$alignment_dir" -name "*${barcode}.adapterstats.txt" -exec grep "^Reads written" {} + | awk -F: '{ print $2 }' | awk '{ print $1 }' | sed 's/,//g')"

        if [ -z "$primer_read_count" ] || [ -z "$adapter_read_count" ]; then
            >&2 echo "Warning: Missing read counts for barcode: $barcode"
            primer_read_count="N/A"
            adapter_read_count="N/A"
        fi

        sorted_bam="$(find "$alignment_dir"/ -name "*${barcode}.sorted.bam")"
        dedup_bam="$(find "$alignment_dir"/ -name "*${barcode}*_dupremoved.sorted.bam")"
        uniq_bam="$(find "$alignment_dir"/ -name "*${barcode}*_dupremoved_uniq.sorted.bam")"

        # Count alignments
        if [ -z "$sorted_bam" ]; then
            mapped_read_count="0"
            echo "No sorted BAM file found for barcode: $barcode, setting mapped read count to 0"
        else
            mapped_read_count="$(samtools view -c "$sorted_bam")"
        fi

        if [ -z "$dedup_bam" ]; then
            dedup_read_count="0"
            echo "No deduplicated BAM file found for barcode: $barcode, setting dedup read count to 0"
        else
            dedup_read_count="$(samtools flagstat "$dedup_bam" | grep "mapped" | head -n 1 | awk '{ print $1 }')"
        fi

        if [ -z "$uniq_bam" ]; then
            uniq_read_count="0"
            echo "No unique BAM file found for barcode: $barcode, setting unique read count to 0"
        else
            uniq_read_count="$(samtools flagstat "$uniq_bam" | grep "mapped" | head -n 1 | awk '{ print $1 }')"
        fi

        # Output the read counts for this sample
        echo -e "${unique_id}\t${flowcell}\t${mapped_read_count}\t${dedup_read_count}\t${uniq_read_count}" >> "$OUTPUT_FILE"
    done
done < "$SAMPLE_MAP_FILE"