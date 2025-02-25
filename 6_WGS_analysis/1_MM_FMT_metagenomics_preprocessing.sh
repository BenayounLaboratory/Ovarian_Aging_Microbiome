#!/bin/bash

#############################
# Menopause-microbiome project
# Pre-process WGS metagenomics - FMT cohort
# Trimmomatic v. 0.39
#############################

# Define paths
TRIMMOMATIC_JAR="/Users/minhookim/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/Users/minhookim/Programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
INPUT_DIR="/Users/minhookim/Data/00_NGS_Rawdata_files/Metagenomics/usftp21.novogene.com/01.RawData"
OUTPUT_DIR="/Users/minhookim/Data/Metagenomics/1_Trimmed_fq"

mkdir -p "$OUTPUT_DIR"

# Trimmomatic parameters
THREADS=4
TRIMMOMATIC_PARAMS="ILLUMINACLIP:$ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50"

# Loop through each folder
for SAMPLE_DIR in MHK271 MHK258 MHK259 MHK260 MHK261 MHK262 MHK263 MHK264 MHK266 MHK267 MHK268 MHK269 MHK270; do

    R1=$(find "$INPUT_DIR/$SAMPLE_DIR" -type f -name "*_1.fq.gz")
    R2=$(find "$INPUT_DIR/$SAMPLE_DIR" -type f -name "*_2.fq.gz")

    # Define output files
    R1_PAIRED="$OUTPUT_DIR/${SAMPLE_DIR}_R1_paired.fq.gz"
    R1_UNPAIRED="$OUTPUT_DIR/${SAMPLE_DIR}_R1_unpaired.fq.gz"
    R2_PAIRED="$OUTPUT_DIR/${SAMPLE_DIR}_R2_paired.fq.gz"
    R2_UNPAIRED="$OUTPUT_DIR/${SAMPLE_DIR}_R2_unpaired.fq.gz"

    # Run Trimmomatic
    echo "Processing $SAMPLE_DIR..."
    java -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" \
        "$R1" "$R2" \
        "$R1_PAIRED" "$R1_UNPAIRED" \
        "$R2_PAIRED" "$R2_UNPAIRED" \
        $TRIMMOMATIC_PARAMS
    echo "Finished processing $SAMPLE_DIR."
done

echo "All samples processed!"