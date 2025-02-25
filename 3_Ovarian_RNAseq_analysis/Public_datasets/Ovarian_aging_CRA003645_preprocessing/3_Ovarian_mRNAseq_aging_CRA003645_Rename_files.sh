#!/bin/bash

####################################
# Menopause-microbiome project
# Ovarian aging mRNA-seq - CRA003645
# https://academic.oup.com/hmg/article/30/21/1941/6302459
# Re-name file names for compatibility with STAR
####################################

cd ./mRNA-seq/Ovarian_aging_CRA003645/02_Trimmed_reads

for file in *_f1_hardtrim.fq.gz; do
  # Use parameter expansion to generate the new name
  newname="${file/_f1_/_1_}"
  echo "Renaming $file to $newname"
  mv "$file" "$newname"
done

for file in *_r2_hardtrim.fq.gz; do
  newname="${file/_r2_/_2_}"
  echo "Renaming $file to $newname"
  mv "$file" "$newname"
done