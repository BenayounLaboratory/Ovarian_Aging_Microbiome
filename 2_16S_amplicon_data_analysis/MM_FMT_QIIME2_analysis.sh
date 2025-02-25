######################################
# Amplicon analysis of Aging cohort males
# By cohort analysis
# QIIME2 v 2023.7
######################################

# FWD_341F: CCTAYGGGRBGCASCAG
# REV_806R: GGACTACNVGGGTWTCTAAT

# conda activate qiime2-2023.7

######################################
# Import data
######################################

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./1_Manifest_and_metadata/MM_FMT_pe-33-manifest.txt \
  --output-path MM_FMT.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Demultiplex - summary

qiime demux summarize \
  --i-data MM_FMT.qza  \
  --o-visualization ./2_Demux_summary/MM_FMT_demux_summarize.qzv

# Denoise - DADA2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs MM_FMT.qza \
  --p-trunc-len-f 226 \
  --p-trunc-len-r 224 \
  --o-table ./3_Feature_table_and_rep_seq_table/MM_FMT_table.qza  \
  --o-representative-sequences ./3_Feature_table_and_rep_seq_table/MM_FMT_rep-seqs.qza \
  --o-denoising-stats ./3_Feature_table_and_rep_seq_table/MM_FMT_denoising-stats.qza  \
  --p-n-threads 4
  
# Summarize feature table and rep seqs

qiime feature-table summarize \
  --i-table ./3_Feature_table_and_rep_seq_table/MM_FMT_table.qza  \
  --o-visualization ./3_Feature_table_and_rep_seq_table/MM_FMT_table.qzv  \
  --m-sample-metadata-file ./1_Manifest_and_metadata/MM_FMT_sample-metadata.txt
  
qiime feature-table tabulate-seqs \
  --i-data ./3_Feature_table_and_rep_seq_table/MM_FMT_rep-seqs.qza \
  --o-visualization ./3_Feature_table_and_rep_seq_table/MM_FMT_rep-seqs.qzv
  

######################################
# Alpha and beta diversity analysis
######################################

# Generating a phylogenetic tree for diversity analysis

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ./3_Feature_table_and_rep_seq_table/MM_FMT_rep-seqs.qza  \
  --o-alignment ./4_Alpha_and_beta_diversity_analysis/MM_FMT_rep-seqs_aligned.qza \
  --o-masked-alignment ./4_Alpha_and_beta_diversity_analysis/MM_FMT_masked-aligned-rep-seqs.qza \
  --o-tree ./4_Alpha_and_beta_diversity_analysis/MM_FMT_unrooted-tree.qza \
  --o-rooted-tree ./4_Alpha_and_beta_diversity_analysis/MM_FMT_rooted-tree.qza


qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ./4_Alpha_and_beta_diversity_analysis/MM_FMT_rooted-tree.qza  \
  --i-table ./3_Feature_table_and_rep_seq_table/MM_FMT_table.qza  \
  --p-sampling-depth 2000 \
  --m-metadata-file ./1_Manifest_and_metadata/MM_FMT_sample-metadata.txt \
  --output-dir ./4_Alpha_and_beta_diversity_analysis/MM_FMT_diversity-core-metrics-phylogenetic

  
######################################
# Taxonomy 
######################################

# SILVA

qiime feature-classifier classify-sklearn \
  --i-classifier ../Common_files/silva_classifier.qza  \
  --i-reads MM_FMT_rep-seqs.qza  \
  --o-classification silva_MM_FMT_taxonomy.qza

  
######################################
# Differential abundance analysis - aldex2
######################################

qiime aldex2 aldex2 \
    --i-table ./3_Feature_table_and_rep_seq_table/MM_FMT_table.qza  \
    --m-metadata-file ./1_Manifest_and_metadata/MM_FMT_sample-metadata.txt  \
    --m-metadata-column age_cat \
    --output-dir ./6_Differential_abundance_and_functional_prediction_analysis/aldex2_AC_males_cohort16


######################################
# Functional abundance prediction analysis - picrust2
######################################

qiime picrust2 full-pipeline \
   --i-table ./3_Feature_table_and_rep_seq_table/MM_FMT_table.qza  \
   --i-seq ./3_Feature_table_and_rep_seq_table/MM_FMT_rep-seqs.qza  \
   --output-dir ./6_Differential_abundance_and_functional_prediction_analysis/q2-picrust2_output_FMT \
   --p-placement-tool sepp \
   --p-threads 3 \
   --p-hsp-method pic \
   --p-max-nsti 2 \
   --verbose   
