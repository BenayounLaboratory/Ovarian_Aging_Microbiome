README -- Ovarian RNAseq analysis
############################

	- FMT_cohort      : 1_Ovarian_mRNAseq_FMT_cohort_QC_and_preprocessing.sh  : Bash script for QC and pre-processing FASTQ files using FASTX Toolkit and Trimgalore
                            2_Ovarian_mRNAseq_FMT_cohort_STAR.sh                  : Bash script for mapping using STAR (GRCm39)
                            3_Ovarian_mRNAseq_FMT_cohort_featureCounts.sh         : Bash script for featureCounts
                            4_Ovarian_mRNAseq_FMT_cohort_DESeq2.R                 : R script for DESeq2 analysis
                            5_Ovarian_mRNAseq_FMT_cohort_GO.R                     : R script for GO analysis and plot generation
                            6_Ovarian_mRNAseq_FMT_cohort_Deconvolution.R          : R script for deconvolution using ovarian scRNAseq dataset (Broad Institute Single Cell Portal: SCP1914) and plot generation
                            7_Ovarian_mRNAseq_FMT_cohort_TF_candidate_GSEA.R      : R script for GSEA analysis of candidate transcription factors
                           
                            MM_FMT_ovary_mRNAseq_featureCounts_output.txt         : Output text file from featureCounts

                            DESeq2                                                : Output files from DESeq2 analysis
                            GO                                                    : Output files from GO analysis
                            Deconvolution                                         : Output files from deconvolution analysis
                            TF_peaks_files                                        : Transcription factor peak files for GSEA analysis

        - Public_datasets : 
                            - Ovarian_aging_CRA003645_preprocessing                            
                                 1_Ovarian_mRNAseq_aging_CRA003645_download_fastq_files.sh  : Bash script for downloading FASTQ files
                                 2_Ovarian_mRNAseq_aging_CRA003645_QC_and_preprocessing.sh  : Bash script for QC and pre-processing FASTQ files using FASTX Toolkit and Trimgalore
                                 3_Ovarian_mRNAseq_aging_CRA003645_Rename_files.sh          : Bash script for re-naming files to be compatible with STAR
                                 4_Ovarian_mRNAseq_aging_CRA003645_STAR.sh                  : Bash script for mapping using STAR (GRCm39)
                                 5_Ovarian_mRNAseq_aging_CRA003645_featureCounts.sh         : Bash script for featureCounts
                                 6_Ovarian_mRNAseq_aging_CRA003645_rename_files.sh          : Bash script for re-naming files to match sample name
                                 7_Ovarian_mRNAseq_aging_CRA003645_DESeq2.R                 : R script for DESeq2 analysis
                           
                                 MM_CRA003645_ovary_aging_mRNAseq_featureCounts_output.txt  : Output text file from featureCounts
                                 MM_CRA003645_ovary_aging_mRNAseq_filenames_mapping.txt     : Mapping file used to re-name file names

                                 DESeq2                                                     : Output files from DESeq2 analysis

                            - Ovarian_aging_PRJNA1002222_preprocessing
                                 1_Ovarian_mRNAseq_aging_PRJNA1002222_download_fastq_files.sh  : Bash script for downloading FASTQ files
                                 2_Ovarian_mRNAseq_aging_PRJNA1002222_QC_and_preprocessing.sh  : Bash script for QC and pre-processing FASTQ files using FASTX Toolkit and Trimgalore
                                 3_Ovarian_mRNAseq_aging_PRJNA1002222_STAR.sh                  : Bash script for mapping using STAR (GRCm39)
                                 4_Ovarian_mRNAseq_aging_PRJNA1002222_featureCounts.sh         : Bash script for featureCounts
                                 5_Ovarian_mRNAseq_aging_PRJNA1002222_rename_files.sh          : Bash script for re-naming files to match sample name

                                 MM_PRJNA1002222_ovary_aging_mRNAseq_featureCounts_output.txt  : Output text file from featureCounts

                                 DESeq2                                                     : Output files from DESeq2 analysis

                            - Combine_CRA003645_and_PRJNA1002222
                                2024-04-08_MM_Filter_ovarian_aging_sig_genes_from_public_data.R : R script used to filter significantly differentially expressed genes
                                2024-04-08_MM_Ovarian_aging_public_data_GSEA.R                  : R script used to perform GSEA

