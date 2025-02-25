#!/usr/bin/bash

####################################
# Menopause-microbiome project
# Ovarian aging mRNA-seq - CRA003645
# https://academic.oup.com/hmg/article/30/21/1941/6302459
# featureCounts
# featureCounts v2.0.3
####################################

#featureCounts

cd ./mRNA-seq/Ovarian_aging_CRA003645/03_STAR

bam_names=$(ls -p | grep "Aligned.sortedByCoord.out.bam" | tr '\n' ' ')

featureCounts -O -p -M --fraction -a ./GCF_000001635.27_GRCm39_genomic.gtf -o ./mRNA-seq/MM_FMT_ovary_RNAseq/04_featureCounts_output/MM_CRA003645_ovary_aging_mRNAseq_featureCounts_output.txt $bam_names