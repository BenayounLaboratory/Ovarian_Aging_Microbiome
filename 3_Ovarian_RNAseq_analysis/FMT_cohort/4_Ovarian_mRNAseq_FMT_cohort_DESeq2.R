setwd('/Users/minhookim/Dropbox/Benayoun_lab/Projects/Menopause_microbiome_project/Data/For_Github/Ovarian_mRNAseq/FMT_cohort/DESeq2')
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
library(DescTools)

###############################################
# Menopause-Microbiome project
# Ovarian mRNA-seq - FMT cohort (FMT-YF vs. FMT-EF)
# Pre-processing and DESeq2
###############################################

##### read in count matrix #####

my.ovary.counts1 <- read.table("../MM_FMT_ovary_mRNAseq_featureCounts_output.txt", sep="\t", header=T)
my.ovary.counts2 <- my.ovary.counts1[,-c(2:6)]

colnames(my.ovary.counts2) <- c("GeneName", 
                                "FMT_YF_1",
                                "FMT_YF_2",
                                "FMT_EF_1",
                                "FMT_EF_2",
                                "FMT_EF_3",
                                "FMT_YF_3",
                                "FMT_YF_4",
                                "FMT_YF_5",
                                "FMT_YF_6",
                                "FMT_YF_7",
                                "FMT_EF_4",
                                "FMT_EF_5",
                                "FMT_EF_6",
                                "FMT_EF_7")

my.ovary.counts <- my.ovary.counts2[,c("GeneName",
                                       "FMT_YF_1",
                                       "FMT_YF_2",
                                       "FMT_YF_3",
                                       "FMT_YF_4",
                                       "FMT_YF_5",
                                       "FMT_YF_6",
                                       "FMT_YF_7",
                                       "FMT_EF_1",
                                       "FMT_EF_2",
                                       "FMT_EF_3",
                                       "FMT_EF_4",
                                       "FMT_EF_5",
                                       "FMT_EF_6",
                                       "FMT_EF_7")]

my.RINs <- c(7.7, 7, 6.9, 7.4, 5.9, 7.3, 7.5, 7.3, 7.6, 7.1, 7.7, 7.4, 7.6, 7.6)
my.FMT <- c(rep("YF", 7), rep("EF", 7))

rownames(my.ovary.counts) <- my.ovary.counts$GeneName

# Get the genes with no reads out
my.good <- which(apply(my.ovary.counts[,-1]>0, 1, sum) >= 6) # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.ovary.counts[my.good,-1] # 29426  genes

################################################################################
# 1. Run SVA to remove unwanted variation
# build design matrix

dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                         fmt = my.FMT,
                         rin = my.RINs)

# Set null and alternative models (ignore batch)
mod1 = model.matrix(~ fmt + rin, data = dataDesign)
n.sv.be = num.sv(my.filtered.matrix, mod1, method="be")    # 3

# Apply SVAseq algortihm
my.svseq = svaseq(as.matrix(my.filtered.matrix), mod1, n.sv=n.sv.be, constant = 0.1)

# remove RIN and SV, preserve age and sex
my.clean <- removeBatchEffect(log2(my.filtered.matrix + 0.1), 
                              batch=NULL, 
                              covariates=cbind(my.svseq$sv,my.RINs),
                              design=mod1[,1:2])

# De-log and round data for DEseq2 processing
my.filtered.sva <- round(2^my.clean-0.1)
write.table(my.filtered.sva, file = paste(Sys.Date(),"MM_FMT_Ovaries_STAR_postSVA_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Check count value range
range(my.filtered.sva)
# [1] 0.000000e+00 3.021079e+13

# Winsorize data - remove extreme values
winsorized_data <- apply(my.filtered.sva, 2, Winsorize, probs = c(0, 0.999))
range(winsorized_data)
# [1]          0 1258429556

# Convert all values in the matrix to rounded integers
integer_winsorized_data <- matrix(as.integer(round(winsorized_data)), 
                                  nrow=nrow(winsorized_data))

colnames(integer_winsorized_data) <- colnames(winsorized_data)
rownames(integer_winsorized_data) <- rownames(winsorized_data)

# Check the result
print(integer_winsorized_data)

# Save data
write.table(integer_winsorized_data, file = paste(Sys.Date(),"MM_FMT_Ovaries_STAR_postSVA_Counts_winsorized.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

################################################################################
# 2. DESeq2 on cleaned data

my.outprefix <- paste(Sys.Date(),"MM_FMT_Ovaries_DESeq2_Analysis_STAR_post_SVA",sep="_")

# design matrix
dataDesign = data.frame( row.names = colnames( integer_winsorized_data ), 
                         fmt = my.FMT,
                         rin = my.RINs)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = integer_winsorized_data,
                              colData = dataDesign,
                              design = ~ fmt)

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
my.colors <- rep("#9C509F", 14)
my.colors[grep("FMT_EF",colnames(tissue.cts))] <- "#F99F27"

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
     xlim = c(-0.02,0.02),
     ylim = c(-0.025,0.025),
     cex.lab = 1.5,
     cex.axis = 1.5)
legend("bottomleft", legend = c("FMT_YF", "FMT_EF"), col = c("#9C509F", "#F99F27"), pch = 16)
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts,1,var) > 0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot_PC1_vs_PC2.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     xlim = c(-175,175),
     ylim = c(-175,175),
     cex.lab = 1.5,
     cex.axis = 1.5) 
legend("bottomleft", legend = c("FMT_YF", "FMT_EF"), col = c("#9C509F", "#F99F27"), pch = 16)
dev.off()

x <- my.pca$x[,1]
y <- my.pca$x[,3]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot_PC1_vs_PC3.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC3 (', round(100*my.summary$importance[,3][2],1),"%)", sep=""),
     xlim = c(-175,175),
     ylim = c(-175,175),
     cex.lab = 1.5,
     cex.axis = 1.5) 
legend("bottomleft", legend = c("FMT_YF", "FMT_EF"), col = c("#9C509F", "#F99F27"), pch = 16)
dev.off()

# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

################################################################################
# 3. Differential gene expression from DESeq2 output

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"MM_FMT_Ovaries_DESeq2_model_with_FMT",sep="_")

res.FMT <- results(dds.deseq, contrast = c("fmt", "EF", "YF"))       # EF / YF

### get the heatmap of aging changes at FDR5; exclude NA
res.FMT <- res.FMT[!is.na(res.FMT$padj),]

genes.FMT <- rownames(res.FMT)[res.FMT$padj < 0.05]
my.num.FMT <- length(genes.FMT)

my.heatmap.out <- paste(my.outprefix, "FMT_Heatmap_FDR5.pdf", sep = "_")
pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("FMT significant (FDR<5%), ", my.num.FMT, " genes", sep="")
pheatmap(tissue.cts[genes.FMT,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, cellwidth = 15)
dev.off()

save(res.FMT, file = paste(Sys.Date(),"MM_FMT_Ovaries_FMT.RData", sep ="_"))

# output result tables of combined analysis to text files
my.out.ct.mat <- paste(my.outprefix,"_log2_counts_matrix_DEseq2_SVA.txt",sep = "_")
write.table(tissue.cts, file = my.out.ct.mat , sep = "\t" , row.names = T, quote=F)

my.out.stats.FMT <- paste(my.outprefix,"FMT_all_genes_statistics.txt",sep = "_")
write.table(res.FMT, file = my.out.stats.FMT , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.FMT <- paste(my.outprefix,"FMT_FDR5_genes_statistics.txt",sep = "_")
write.table(res.FMT[genes.FMT,], file = my.out.fdr5.FMT, sep = "\t" , row.names = T, quote=F)

################################################################################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()