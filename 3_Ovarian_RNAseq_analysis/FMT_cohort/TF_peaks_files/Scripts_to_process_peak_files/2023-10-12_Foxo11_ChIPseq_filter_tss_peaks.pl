# 2023-10-12
# Filter TSS-promoter peaks from Foxo1 ChIP-seq data

use strict;
use warnings;

my $HOMERinfile = "/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/Public_data_TF/Foxo1/Foxo1.peaks.txt";

my $HOMERoutfile = "/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/Public_data_TF/Processed_annotated_peak_files/2023-10-12_Downloaded_Foxo1_ChIPseq_HOMER_tss_peaks.bed";

my $HOMERgeneout = "/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/Public_data_TF/Processed_annotated_peak_files/2023-10-12_Downloaded_Foxo1_ChIPseq_HOMER_tss_peaks_gene_names.txt";

my $count = 0;

# Print header
# Print line if line contains "TSS" and "protein-coding"

# Process HOMER peak file
open(IN, "<$HOMERinfile") or die "Can't open $HOMERinfile\n";
open(OUT, ">$HOMERoutfile") or die "Can't open $HOMERoutfile\n";
open(OUTgene, ">$HOMERgeneout") or die "Can't open $HOMERgeneout\n";

while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /TSS/) {
        if ($line =~ /protein-coding/) {
            print OUT "$line\n";
            # Extract gene name
            my @line = split(/\t/, $line);
            my $gene = $line[15];
            $gene =~ s/.*\((.*)\)/$1/;
            print OUTgene "$gene\n";
            $count++;
        }
    }
}
close IN;
close OUT;
close OUTgene;

print "Number of TSS peaks: $count\n";
# Number of TSS peaks: 1773
