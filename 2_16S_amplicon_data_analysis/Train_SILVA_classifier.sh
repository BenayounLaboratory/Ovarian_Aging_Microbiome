######################################
# Menopause-microbiome project
# Train SILVA classifier for taxonomy analysis
# QIIME2 v 2023.7
######################################

# FWD_341F: CCTAYGGGRBGCASCAG
# REV_806R: GGACTACNVGGGTWTCTAAT

# conda activate qiime2-2023.7

######################################
# Import data
######################################

# SILVA

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ./SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna  \
  --output-path silva_99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ./SILVA_132_QIIME_release/taxonomy/16S_only/99/taxonomy_all_levels.txt  \
  --output-path silva_99_ref-taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences silva_99_otus.qza \
  --p-f-primer CCTAYGGGRBGCASCAG \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-min-length 100 \
  --p-max-length 500 \
  --o-reads silva_99_ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva_99_ref-seqs.qza  \
  --i-reference-taxonomy silva_99_ref-taxonomy.qza \
  --o-classifier silva_classifier.qza
