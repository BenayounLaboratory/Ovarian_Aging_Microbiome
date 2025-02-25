#!/bin/bash

INPUT_DIR="/Volumes/MKim_3/MM_FMT_Metagenomics/2_Unmapped_fq"
OUTPUT_DIR="/Volumes/MKim_3/MM_FMT_Metagenomics/7_humann3/humann3_outputs"
THREADS=4

mkdir -p "$OUTPUT_DIR"

LOG_FILE="$OUTPUT_DIR/humann3_run.log"

# Clear the log file
> "$LOG_FILE"

echo "Starting HUMAnN3 analysis on paired-end reads at $(date)" | tee -a "$LOG_FILE"

for FILE1 in "$INPUT_DIR"/*.1.fq.gz; do

    BASENAME=$(basename "$FILE1" .1.fq.gz)
    FILE2="$INPUT_DIR/${BASENAME}.2.fq.gz" 
    
    if [[ -f "$FILE2" ]]; then
        SAMPLE_OUTPUT_DIR="$OUTPUT_DIR/$BASENAME"

        mkdir -p "$SAMPLE_OUTPUT_DIR"

        echo "Processing $BASENAME with paired-end files $FILE1 and $FILE2..." | tee -a "$LOG_FILE"

        # Run HUMAnN3
        humann --input "$FILE1" \
               --input "$FILE2" \
               --output "$SAMPLE_OUTPUT_DIR" \
               --threads "$THREADS" 2>> "$LOG_FILE"

        if [ $? -eq 0 ]; then
            echo "Successfully processed $BASENAME" | tee -a "$LOG_FILE"
        else
            echo "Error processing $BASENAME. Check the log for details." | tee -a "$LOG_FILE"
        fi
    else
        echo "Warning: Missing reverse read for $FILE1. Skipping this sample." | tee -a "$LOG_FILE"
    fi
done

echo "HUMAnN3 paired-end processing completed at $(date)" | tee -a "$LOG_FILE"

# Merge pathway abundance tables
echo "Merging pathway abundance tables..." | tee -a "$LOG_FILE"
humann_join_tables --input "$OUTPUT_DIR" --output "$OUTPUT_DIR/merged_pathabundance.tsv" --file_name pathabundance 2>> "$LOG_FILE"

# Merge gene family tables
echo "Merging gene family tables..." | tee -a "$LOG_FILE"
humann_join_tables --input "$OUTPUT_DIR" --output "$OUTPUT_DIR/merged_genefamilies.tsv" --file_name genefamilies 2>> "$LOG_FILE"

# Normalize the merged tables
echo "Normalizing merged tables..." | tee -a "$LOG_FILE"
humann_renorm_table --input "$OUTPUT_DIR/merged_pathabundance.tsv" --output "$OUTPUT_DIR/normalized_pathabundance.tsv" --units cpm 2>> "$LOG_FILE"
humann_renorm_table --input "$OUTPUT_DIR/merged_genefamilies.tsv" --output "$OUTPUT_DIR/normalized_genefamilies.tsv" --units cpm 2>> "$LOG_FILE"

echo "HUMAnN3 analysis completed. Results are in $OUTPUT_DIR" | tee -a "$LOG_FILE"
