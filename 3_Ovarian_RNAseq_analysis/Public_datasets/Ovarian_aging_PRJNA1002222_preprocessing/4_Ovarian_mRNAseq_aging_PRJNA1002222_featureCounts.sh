#!/usr/bin/bash

####################################
# Menopause-microbiome project
# Ovarian aging mRNA-seq - PRJNA1002222
# https://www.nature.com/articles/s43587-023-00532-9
# featureCounts
# featureCounts v2.0.3
####################################

#featureCounts

cd ./mRNA-seq/Ovarian_aging_PRJNA1002222/03_STAR

bam_names=$(ls -p | grep "Aligned.sortedByCoord.out.bam" | tr '\n' ' ')

featureCounts -O -p -M --fraction -a ./GCF_000001635.27_GRCm39_genomic.gtf -o ./mRNA-seq/MM_FMT_ovary_RNAseq/04_featureCounts_output/MM_PRJNA1002222_ovary_aging_mRNAseq_featureCounts_output.txt $bam_names