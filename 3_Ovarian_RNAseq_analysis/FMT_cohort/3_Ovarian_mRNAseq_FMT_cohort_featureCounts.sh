#!/usr/bin/bash

####################################
# Menopause-microbiome project
# Ovarian mRNA-seq - FMT cohort (FMT-YF vs. FMT-EF)
# featureCounts
# featureCounts v2.0.3
####################################

#featureCounts

cd ./mRNA-seq/MM_FMT_ovary_RNAseq/03_STAR

bam_names=$(ls -p | grep "Aligned.sortedByCoord.out.bam" | tr '\n' ' ')

featureCounts -O -p -M --fraction -a ./GCF_000001635.27_GRCm39_genomic.gtf -o ./mRNA-seq/MM_FMT_ovary_RNAseq/04_featureCounts_output/MM_FMT_ovary_mRNAseq_featureCounts_output.txt $bam_names