#!/usr/bin/bash

####################################
# Menopause-microbiome project
# Ovarian aging mRNA-seq - CRA003645
# https://academic.oup.com/hmg/article/30/21/1941/6302459
# STAR alignment - GRCm39
# STAR 2.7.9a
####################################

Indexed_Genome_Dir=./STAR_index/GRCm39

FastqDir=./mRNA-seq/Ovarian_aging_CRA003645/01_Rawdata
FastxDir=./mRNA-seq/Ovarian_aging_CRA003645/02_Trimmed_reads

trimgalout=$FastxDir/"quality"

Output_directory=./mRNA-seq/Ovarian_aging_CRA003645/03_STAR

cd $trimgalout
	
for f in $(find "." -name '*_1_hardtrim_val_1.fq.gz')
	do
    	f2=$(basename "${f}" | sed 's/_1_hardtrim_val_1\.fq\.gz/_2_hardtrim_val_2\.fq\.gz/g');
    	inf2="${trimgalout}/${f2}"
    	of=$(basename "${f}" | sed 's/_1_hardtrim_val_1\.fq\.gz//g');
    	oFname="${Output_directory}/${of}"
    	STAR --genomeDir $Indexed_Genome_Dir --readFilesIn $f $inf2 --readFilesCommand gzcat --runThreadN 2 --outFilterMultimapNmax 200 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --alignEndsProtrude 10 ConcordantPair --limitGenomeGenerateRAM 60000000000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $oFname
	done
			
	echo "Done Mapping Reads to Genome with STAR"
	
fi