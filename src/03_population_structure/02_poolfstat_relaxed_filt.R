# Calculate Fst statistics and generate PCA of all SNPs

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('poolfstat', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'maps', 'mapdata', 'ggrepel', 'pheatmap'))
library(poolfstat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(maps) 
library(mapdata)
library(ggrepel)
library(pheatmap)

# ================================================================================== #

# Read in population names
pops <- read.table("guide_files/N.canaliculata_relaxed_filt_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh)
pooldata_relaxed_filt <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean_relaxed_filt/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 2, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 8,277,206 SNPs for 19 Pools

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

# ================================================================================== #

# Generate Folders and files

# Make data directory
data_dir="data/processed/pop_structure"
if (!dir.exists(data_dir)) {dir.create(data_dir)}

# Make data directory for fst
data_dir_fst="data/processed/pop_structure/Fst"
if (!dir.exists(data_dir_fst)) {dir.create(data_dir_fst)}

# Make output directory
output_dir="output/figures/pop_structure"
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# Make output directory for poolfstat
output_dir_poolfstat="output/figures/pop_structure/poolfstat"
if (!dir.exists(output_dir_poolfstat)) {dir.create(output_dir_poolfstat)}

# ================================================================================== #

# Save pooldata
save(pooldata_relaxed_filt, file="data/processed/pop_structure/pooldata_relaxed_filt.RData")
# Reload pooldata
load("data/processed/pop_structure/pooldata_relaxed_filt.RData")

# ================================================================================== #

# Estimate genome wide Fst 

# Use computeFST function (default to using the Anova method)
pooldata_relaxed_filt.fst <- computeFST(pooldata_relaxed_filt,verbose=FALSE)
pooldata_relaxed_filt.fst$Fst 
#     ####0.5801261

# Block-Jackknife estimation of Fst standard error and confidence intervals
pooldata_relaxed_filt.fst.bjack <- computeFST(pooldata_relaxed_filt, nsnp.per.bjack.block = 1000, verbose=FALSE)
pooldata_relaxed_filt.fst.bjack$Fst
#   Estimate  bjack mean  bjack s.e.     CI95inf     CI95sup 
#      ####0.580126095 0.582574371 0.001481405 0.579670817 0.585477924 

# Compute multi-locus Fst over sliding window of SNPs
pooldata_relaxed_filt.fst.sliding.window <- computeFST(pooldata_relaxed_filt, sliding.window.size=1000)

# Plot sliding window
pdf("output/figures/pop_structure/poolfstat/fst.sliding.window_relaxed_filt.pdf", width = 14, height = 8)
par(mar=c(4.5, 5, 5, 1))
plot(pooldata_relaxed_filt.fst.sliding.window$sliding.windows.fvalues$CumMidPos/1e6, 
pooldata_relaxed_filt.fst.sliding.window$sliding.windows.fvalues$MultiLocusFst,
xlab="Cumulated Position (in Mb)", ylab="Multi-locus Fst", pch=16, cex.lab=2)
abline(h=pooldata_relaxed_filt.fst.sliding.window$Fst,lty=2, lwd=4, col="red") # Dashed red line indicates the estimated overall genome-wide Fst
dev.off()

# ================================================================================== #

# Estimate pairwise-population Fst

# Use compute.pairwiseFST function
pooldata_relaxed_filt.pairwisefst <- compute.pairwiseFST(pooldata_relaxed_filt, verbose=FALSE)

# Graph heatmap
pdf("output/figures/pop_structure/poolfstat/heatmap.pdf", width = 8, height = 8)
heatmap(pooldata_relaxed_filt.pairwisefst, cexRow=1.7, cexCol=1.7)
dev.off()

# Block-Jackknife estimation of pairwise Fst standard error and confidence intervals
pooldata_relaxed_filt.pairwisefst.bjack <- compute.pairwiseFST(pooldata_relaxed_filt, nsnp.per.bjack.block = 1000, verbose=FALSE)

# Estimated pairwise Fst are stored in the slot values: 5 first estimated pairwise
head(pooldata_relaxed_filt.pairwisefst.bjack@values)

# Graph estimated pairwise-population FST with their 95% confidence intervals 
pdf("output/figures/pop_structure/poolfstat/pairwise_Fst_relaxed_filt.pdf", width = 8, height = 8)
plot(pooldata_relaxed_filt.pairwisefst.bjack, cex=0.5)
dev.off()

# Save pairwise matrix as xls for follow-up analyses
pooldata_relaxed_filt.pairwisefst.matrix <- pooldata_relaxed_filt.pairwisefst@PairwiseFSTmatrix
pooldata_relaxed_filt.pairwisefst.matrix <- as.data.frame(pooldata_relaxed_filt.pairwisefst.matrix)
write.csv(pooldata_relaxed_filt.pairwisefst.matrix, "data/processed/pop_structure/Fst/pooldata_relaxed_filt.pairwisefst.matrix.csv")

# Graph Fst heatmap with pheatmap 
pdf("output/figures/pop_structure/poolfstat/pheatmap_Fst_relaxed_filt.pdf", width = 8, height = 8)
pheatmap(pooldata_relaxed_filt.pairwisefst.matrix, border_color = "black", fontsize_col = 16, fontsize_row = 16)
dev.off()
pdf("output/figures/pop_structure/poolfstat/pheatmap_Fst_alt_relaxed_filt.pdf", width = 8, height = 8)
pheatmap(pooldata_relaxed_filt.pairwisefst.matrix, border_color = "black", fontsize_col = 16, fontsize_row = 16, color=rev(viridis(n=361, alpha=1,begin=0, end=1, option="viridis")))
dev.off()

# ================================================================================== #

# Principle Components Analysis with randomallele.pca

#PCA on the read count data (the object)
pooldata_relaxed_filt.pca = randomallele.pca(pooldata_relaxed_filt, main="Read Count data")

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
colors.reorder <- mycolors[c(19,2,3,4,11,5,1,10,17,14,6,8,7,13,9,18,16,12,15)]

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC1_PC2_relaxed_filt.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata_relaxed_filt.pca$pop.loadings[,1],pooldata_relaxed_filt.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata_relaxed_filt.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata_relaxed_filt.pca$perc.var[2],2),"%)"),
col="black", bg=colors.reorder, pch=21, cex = 3, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC3_PC4_relaxed_filt.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata_relaxed_filt.pca$pop.loadings[,3],pooldata_relaxed_filt.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata_relaxed_filt.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata_relaxed_filt.pca$perc.var[4],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC5 and PC6
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC5_PC6.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata_relaxed_filt.pca$pop.loadings[,5],pooldata_relaxed_filt.pca$pop.loadings[,6],
xlab=paste0("PC",5," (",round(pooldata_relaxed_filt.pca$perc.var[5],2),"%)"),
ylab=paste0("PC",6," (",round(pooldata_relaxed_filt.pca$perc.var[6],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plot names on PCAs

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC1_PC2_names.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata_relaxed_filt.pca$pop.loadings[,1],pooldata_relaxed_filt.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata_relaxed_filt.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata_relaxed_filt.pca$perc.var[2],2),"%)"))
text(pooldata_relaxed_filt.pca$pop.loadings[,1], pooldata_relaxed_filt.pca$pop.loadings[,2], pooldata_relaxed_filt@poolnames, cex=0.5, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC3_PC4_names.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata_relaxed_filt.pca$pop.loadings[,3],pooldata_relaxed_filt.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata_relaxed_filt.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata_relaxed_filt.pca$perc.var[4],2),"%)"))
text(pooldata_relaxed_filt.pca$pop.loadings[,3], pooldata_relaxed_filt.pca$pop.loadings[,4], pooldata_relaxed_filt@poolnames, cex=0.5, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

