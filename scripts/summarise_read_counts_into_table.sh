#!/bin/bash

SAMPLE_MAP_FILE="data/input/hiseq_sample_id_map.tsv"
OUTPUT_FILE="data/output/summarised_read_counts.tsv"
FLOWCELL_NUM="2"

FASTQ_DIR="data/input/FASTQ"
OUTPUT_DIR="data/output"

while read line; do
    # Extract the sample name and read count from each line
    library_name=$(echo "$line" | awk -F"\t" '{OFS="_"; print $5,"'$FLOWCELL_NUM'","1",$6}')
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

    alignment_files="$(find "$OUTPUT_DIR/alignments" -name "*${barcode}*")"

    if [ -z "$alignment_files" ]; then
        >&2 echo "No alignment files found for barcode: $barcode"
        exit 1
    fi

    # Count reads from FASTQ files and alignment stats
    in_read_count=$(find "$FASTQ_DIR" -name "${library_name}_*.fastq.gz" -exec zcat {} + | wc -l | awk '{ print $1/4 }')
    primer_read_count="$(find "$OUTPUT_DIR/alignments" -name "*${barcode}.ssbstats.txt" -exec grep '^Reads written' {} + | awk -F: '{ print $2 }' | awk '{ print $1 }' | sed 's/,//g')"
    adapter_read_count="$(find "$OUTPUT_DIR/alignments" -name "*${barcode}.adapterstats.txt" -exec grep "^Reads written" {} + | awk -F: '{ print $2 }' | awk '{ print $1 }' | sed 's/,//g')"

    if [ -z "$primer_read_count" ] || [ -z "$adapter_read_count" ]; then
        >&2 echo "Error: Missing read counts for barcode: $barcode"
        exit 1
    fi

    # Count alignments
    mapped_read_count="$(find "$OUTPUT_DIR/alignments" -name "*${barcode}.sorted.bam.flagstat" -exec head -n 7 {} + | tail -n 1 | awk '{ print $1 }')"
    dedup_read_count="$(find "$OUTPUT_DIR/alignments" -name "*${barcode}_dupremoved.sorted.bam.flagstat" -exec head -n 7 {} + | tail -n 1 | awk '{ print $1 }')"
    uniq_read_count="$(find "$OUTPUT_DIR/alignments" -name "*${barcode}_dupremoved_uniq.sorted.bam.flagstat" -exec head -n 7 {} + | tail -n 1 | awk '{ print $1 }')"

    echo -e "${unique_id}\t${in_read_count}\t${primer_read_count}\t${adapter_read_count}\t${mapped_read_count}\t${dedup_read_count}\t${uniq_read_count}"
done < "$SAMPLE_MAP_FILE" > "$OUTPUT_FILE"