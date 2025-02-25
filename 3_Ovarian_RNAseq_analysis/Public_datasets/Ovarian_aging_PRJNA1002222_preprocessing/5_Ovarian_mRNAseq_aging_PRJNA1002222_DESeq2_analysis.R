options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(fields)
library(VennDiagram)

# rm(list = ls())

###############################################
# Menopause-microbiome project
# Public ovarian aging bulk mRNA-seq analysis - PRJNA1002222
# Pre-processing and DESeq2
# Compare gene expression patterns between PRJNA1002222 and FMT ovarian bulk mRNA-seq datasets
###############################################

##### read in count matrix #####

my.ovary.counts1 <- read.table("../MM_PRJNA1002222_ovary_aging_mRNAseq_featureCounts_output.txt",
                               sep="\t", header=T)

my.ovary.counts2 <- my.ovary.counts1[,c(1, 7:14)]

colnames(my.ovary.counts2) <- c("GeneName", 
                                "12M_3", "12M_2", "12M_1",
                                "2M_5", "2M_4", "2M_3", 
                                "2M_2", "2M_1")

my.ovary.counts <- my.ovary.counts2[,c("GeneName",
                                       "2M_1",
                                       "2M_2",
                                       "2M_3",
                                       "2M_4",
                                       "2M_5",
                                       "12M_1",
                                       "12M_2",
                                       "12M_3")]

my.Age <- c(rep("2M", 5), rep("12M", 3))

rownames(my.ovary.counts) <- my.ovary.counts$GeneName

# get the genes with no reads out
my.good <- which(apply(my.ovary.counts[,-1]>0, 1, sum) >= 5) # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.ovary.counts[my.good,-1] # 21709  genes

################################################################################
# 1. Run SVA to remove unwanted variation
# build design matrix

dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                         age = my.Age)

# Set null and alternative models (ignore batch)
mod1 = model.matrix(~ age, data = dataDesign)
n.sv.be = num.sv(my.filtered.matrix, mod1, method="be")

# apply SVAseq algortihm
my.svseq = svaseq(as.matrix(my.filtered.matrix), mod1, n.sv=n.sv.be, constant = 0.1)

# remove RIN and SV, preserve age and sex
my.clean <- removeBatchEffect(log2(my.filtered.matrix + 0.1), 
                              batch=NULL, 
                              design=mod1[,1:2])

# delog and round data for DEseq2 processing
my.filtered.sva <- round(2^my.clean-0.1)
write.table(my.filtered.sva, file = paste(Sys.Date(),"MM_PRJNA1002222_STAR_postSVA_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

################################################################################
# 2. DESeq2 on cleaned data

my.outprefix <- paste(Sys.Date(),"MM_PRJNA1002222_DESeq2_Analysis_STAR_post_SVA",sep="_")

# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.sva ), 
                         age = my.Age)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                              colData = dataDesign,
                              design = ~ age)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix,"dispersion_plot.pdf",sep="_")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

save(tissue.cts, file = paste(my.outprefix, "normalized_counts.RData", sep = "_"))

# color-code 
my.colors <- rep("deeppink1",10)
my.colors[grep("12M",colnames(tissue.cts))] <- "deeppink4"

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.1,0.1),
     ylim = c(-0.1,0.1),
     cex.lab = 1.5,
     cex.axis = 1.5)
legend("bottomleft", legend = c("YF (3m)", "EF (12m)"), col = c("deeppink1", "deeppink4"), pch = 16)
dev.off()

# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

################################################################################
# 3. Differential gene expression from DESeq2 output

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"MM_PRJNA1002222_DESeq2_aging",sep="_")

res.Age <- results(dds.deseq, contrast = c("age", "12M", "2M"))       # OF / YF

### get the heatmap of aging changes at FDR5; exclude NA
res.Age <- res.Age[!is.na(res.Age$padj),]

genes.Age <- rownames(res.Age)[res.Age$padj < 0.05]
my.num.Age <- length(genes.Age)

my.heatmap.out <- paste(my.outprefix, "Age_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Age significant (FDR<5%), ", my.num.Age, " genes", sep="")
pheatmap(tissue.cts[genes.Age,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.Age, file = paste(Sys.Date(),"MM_PRJNA1002222_Age_Ovaries_Age.RData", sep ="_"))

# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2_SVA.txt",sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.Age <- paste(my.outprefix,"Age_all_genes_statistics.txt",sep = "_")
write.table(res.Age, file = my.out.stats.Age , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.Age <- paste(my.outprefix,"Age_FDR5_genes_statistics.txt",sep = "_")
write.table(res.Age[genes.Age,], file = my.out.fdr5.Age, sep = "\t" , row.names = T, quote=F)

################################################################################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()
