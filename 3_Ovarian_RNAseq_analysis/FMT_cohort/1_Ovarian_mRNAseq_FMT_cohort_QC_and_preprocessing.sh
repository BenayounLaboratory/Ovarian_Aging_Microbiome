#!/usr/bin/bash

####################################
# Menopause-microbiome project
# Ovarian mRNA-seq - FMT cohort (FMT-YF vs. FMT-EF)
# Hard trim & QC fastq reads
# FASTX Toolkit 0.0.14 & Trim Galore 0.6.6
####################################

# trimmer

FastqDir=./mRNA-seq/MM_FMT_ovary_RNAseq/01_Rawdata
FastxDir=./mRNA-seq/MM_FMT_ovary_RNAseq/02_Trimmed_reads

cd $FastqDir

for f in $(find "." -name '*.fq.gz')
do
        f2=$(basename "${f}");
        of=$(basename "${f}" | sed 's/\.fq\.gz/_hardtrim\.fq\.gz/g');

        gzcat $f2 | fastx_trimmer -f 7 -l 100 -z -i - -Q33 -o $fastxdir/$of
        
        echo "Done hard-trimming $f2"
        
done

echo "Done with hard trimming!"


# trim_galore

cd $FastxDir

trimgalout=$FastxDir/"quality"

mkdir $trimgalout

for f in $(find "." -name '*1_hardtrim.fq.gz')
	do
        f1=$(basename "${f}");
        f2=$(basename "${f}" | sed 's/1_hardtrim\.fq\.gz/2_hardtrim\.fq\.gz/g');
		trim_galore --paired $fastxdir/$f1 $fastxdir/$f2 -o $trimgalout  
	done
fi

echo "Done trimming adapters!"