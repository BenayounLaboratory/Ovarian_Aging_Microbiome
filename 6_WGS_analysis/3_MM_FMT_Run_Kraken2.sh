#!/bin/bash

#SBATCH --account=bbenayou_34
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=48:00:00

# Set the paths and parameters
KRAKEN_DB="/project/bbenayou_34/kim/MM_FMT_Metagenomics/Kraken2/Kraken2DB"
THREADS=24
INPUT_DIR="/project/bbenayou_34/kim/MM_FMT_Metagenomics/1_Clean_fq"
OUTPUT_DIR="/project/bbenayou_34/kim/MM_FMT_Metagenomics/Kraken2/Kraken_outputs"

# Build DB
/home1/minhooki/Programs/kraken2-master/kraken2-build --standard --threads 24 --db /project/bbenayou_34/kim/MM_FMT_Metagenomics/Kraken2/Kraken2DB

mkdir -p ${OUTPUT_DIR}

for FORWARD in ${INPUT_DIR}/*_unmapped.1.fq.gz; do

  SAMPLE=$(basename ${FORWARD} _unmapped.1.fq.gz)
  REVERSE="${INPUT_DIR}/${SAMPLE}_unmapped.2.fq.gz"
  
  kraken2 --db ${KRAKEN_DB} \
          --threads ${THREADS} \
          --report ${OUTPUT_DIR}/${SAMPLE}.kreport \
          --paired ${FORWARD} ${REVERSE} > ${OUTPUT_DIR}/${SAMPLE}.kraken
done
