#!/bin/bash

####################################
# Menopause-microbiome project
# Ovarian aging mRNA-seq - CRA003645
# https://academic.oup.com/hmg/article/30/21/1941/6302459
# Re-name file names to match metadata
####################################

while IFS=$'\t' read -r original new; do
  mv "$original" "$new"
done < ./MM_CRA003645_ovary_aging_mRNAseq_filenames_mapping.txt

echo "Renaming complete."