#!/bin/bash

#SBATCH --account=bbenayou_34
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=48:00:00

input_dir="/project/bbenayou_34/kim/MM_FMT_Metagenomics/Kraken2/Kraken_outputs"
database_dir="/project/bbenayou_34/kim/MM_FMT_Metagenomics/Kraken2/Kraken2DB"

for kreport_file in "$input_dir"/*.kreport; do

    sample=$(basename "$kreport_file" .kreport)
    
    /home1/minhooki/Programs/Bracken-master/bracken \
        -d "$database_dir" \
        -i "$kreport_file" \
        -o "${input_dir}/${sample}.bracken" \
        -w "${input_dir}/${sample}.bracken.report" \
        -r 150
done