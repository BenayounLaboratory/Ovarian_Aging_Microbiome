README -- WGS analysis (revision)
############################

	- 1_MM_FMT_metagenomics_preprocessing.sh                   : Script for pre-processing WGS dataset
        - 2_MM_FMT_metagenomics_remove_host_reads_Bowtie2.sh       : Script for removing host reads
        - 3_MM_FMT_Run_Kraken2.sh                                  : Script for running Kraken2
        - 4_MM_FMT_Run_Bracken.sh                                  : Script for running Bracken
        - 5_Run_humann3.sh                                         : Script for running HUMAnN3
        - 6_MM_FMT_metagenomics_analysis.R                         : R script for differential abundance analysis
        - 7_MM_FMT_metagenomics_mediation_analysis.R               : R script for mediation analysis

        - WGS_hormone_spearman                                     : Scripts for WGS top 20 species vs. serum hormone (AMH, FSH & INHBA) spearman correlation analysis