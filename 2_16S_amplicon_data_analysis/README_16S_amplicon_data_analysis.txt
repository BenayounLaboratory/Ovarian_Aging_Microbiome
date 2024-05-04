README -- 16S amplicon data analysis
############################

	- MM_AC_females_QIIME2_analysis.sh  : Bash script QIIME2 (v 2023.7) analysis for aging cohort females
	- MM_AC_males_QIIME2_analysis.sh    : Bash script QIIME2 (v 2023.7) analysis for aging cohort males
	- MM_VCD_QIIME2_analysis.sh         : Bash script QIIME2 (v 2023.7) analysis for VCD cohort
	- MM_FMT_QIIME2_analysis.sh         : Bash script QIIME2 (v 2023.7) analysis for FMT cohort
        
        - 1_Manifest_and_metadata                                      : Metadata and manifest files for QIIME2 analysis
        - 2_Demux_summary                                              : De-multiplex summary files from QIIME2 analysis
        - 3_Feature_table_and_rep_seq_table                            : Feature table, representative sequence table and denoising statistics from QIIME2 analysis
        - 4_Alpha_and_beta_diversity_analysis                          :
                                                                           - 0_Generate_plots_R: R scripts for alpha and beta diversity plot generation
        - 5_Taxonomy                                                   : Taxonomy files from QIIME2 analysis
        - 6_Differential_abundance_and_functional_prediction_analysis  :   - Differential_abundance_aldex
                                                                                - aldex2 outputs
                                                                                - R scripts to analyze aldex2 results and plot generation
                                                                           - Functional_prediction_picrust2
                                                                                - picrust2 outputs
                                                                                - R scripts to analyze picrust2 results and plot generation
        - 7_PCA_Post_CLR_transformation                                : R scripts for CLR-transformation of features and PCA plot generation
        - 8_Mediation_analysis                                         : R script for mediation analysis of microbial genera and ovarian transcriptome dataset
        - Common_files                                                 : SILVA reference files used for QIIME2 analysis