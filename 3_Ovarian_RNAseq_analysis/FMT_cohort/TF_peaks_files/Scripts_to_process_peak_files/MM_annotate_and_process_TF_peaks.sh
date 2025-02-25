#!/usr/bin/bash

####################################
# Menopause-microbiome project
# Ovarian mRNA-seq - FMT cohort (FMT-YF vs. FMT-EF)
# Process and annotate downloaded peak files using HOMER
# HOMER v4.11.1
####################################

# HOMER

# Nucks1

annotatePeaks.pl ../Public_data_TF/Nucks1/9_Downloaded_peak_info/48959_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Nucks1/Nucks1.peaks.txt -annStats ../Public_data_TF/Nucks1/annotation_stats_Nucks1_HOMER_annotated_HOMER_peaks.txt

# Foxo1

annotatePeaks.pl ../Public_data_TF/Foxo1/2_Downloaded_peak_info/92461_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Foxo1/Foxo1.peaks.txt -annStats ../Public_data_TF/Foxo1/annotation_stats_Foxo1_HOMER_annotated_HOMER_peaks.txt

# Ikzf1

annotatePeaks.pl ../Public_data_TF/Ikzf1/1_Downloaded_peak_info/57856_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Ikzf1/Ikzf1.peaks.txt -annStats ../Public_data_TF/Ikzf1/annotation_stats_Ikzf1_HOMER_annotated_HOMER_peaks.txt

# Lyl1

annotatePeaks.pl ../Public_data_TF/Lyl1/1_Downloaded_peak_info/85097_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Lyl1/Lyl1.peaks.txt -annStats ../Public_data_TF/Lyl1/annotation_stats_Lyl1_HOMER_annotated_HOMER_peaks.txt

# Mllt3

annotatePeaks.pl ../Public_data_TF/Mllt3/1_Downloaded_peak_info/52855_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Mllt3/Mllt3.peaks.txt -annStats ../Public_data_TF/Mllt3/annotation_stats_Mllt3_HOMER_annotated_HOMER_peaks.txt

# Ncoa1

annotatePeaks.pl ../Public_data_TF/Ncoa1/1_Downloaded_peak_info/44213_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Ncoa1/Ncoa1.peaks.txt -annStats ../Public_data_TF/Ncoa1/annotation_stats_Ncoa1_HOMER_annotated_HOMER_peaks.txt

# Usp7

annotatePeaks.pl ../Public_data_TF/Usp7/1_Downloaded_peak_info/72322_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Usp7/Usp7.peaks.txt -annStats ../Public_data_TF/Usp7/annotation_stats_Usp7_HOMER_annotated_HOMER_peaks.txt

#Ing1

annotatePeaks.pl ../Public_data_TF/Ing1/1_Downloaded_peak_info/44360_peaks.bed ./GCF_000001635.27_GRCm39_genomic.gtf > ../Public_data_TF/Ing1/Ing1.peaks.txt -annStats ../Public_data_TF/Ing1/annotation_stats_Ing1_HOMER_annotated_HOMER_peaks.txt
