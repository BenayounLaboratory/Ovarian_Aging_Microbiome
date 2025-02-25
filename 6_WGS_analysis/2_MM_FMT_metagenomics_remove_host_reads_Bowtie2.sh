#!/bin/bash

#############################
# Menopause-microbiome project
# WGS metagenomics - FMT cohort
# Remove host DNA
# Bowtie2 v. 2.4.4 & mm10
#############################

BOWTIE2_INDEX="/Users/minhookim/Data/01_Reference/Bowtie2/mm10/mm10"
INPUT_DIR="/Users/minhookim/Data/Metagenomics/1_Trimmed_fq"
OUTPUT_DIR="/Users/minhookim/Data/Metagenomics/2_Unmapped_fq"
THREADS=4

mkdir -p "$OUTPUT_DIR"

for R1 in "$INPUT_DIR"/*_R1_paired.fq.gz; do

    R2="${R1/_R1_paired/_R2_paired}"

    SAMPLE=$(basename "$R1" | sed 's/_R1_paired.*//')

    if [[ ! -f "$R2" ]]; then
        echo "Error: Missing R2 file for $SAMPLE. Skipping..."
        continue
    fi

    UNMAPPED_R1="$OUTPUT_DIR/${SAMPLE}_R1_unmapped.fq.gz"
    UNMAPPED_R2="$OUTPUT_DIR/${SAMPLE}_R2_unmapped.fq.gz"
    SAM_OUT="$OUTPUT_DIR/${SAMPLE}.sam"

    echo "Processing $SAMPLE with paired-end alignment..."
    bowtie2 -x "$BOWTIE2_INDEX" -1 "$R1" -2 "$R2" \
        --threads "$THREADS" --very-sensitive --un-conc-gz "$OUTPUT_DIR/${SAMPLE}_unmapped.fq.gz" \
        -S "$SAM_OUT"
    echo "Finished processing $SAMPLE."
done
