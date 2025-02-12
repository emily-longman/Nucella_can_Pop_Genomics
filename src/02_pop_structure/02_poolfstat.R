# Calculate Fst statistics and gnerate PCA of all SNPs

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
install.packages(c('poolfstat', 'ggplot2', 'RColorBrewer'))
library(poolfstat)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/pop_structure/guide_files/Nucella_pops.list", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf.gz" ,poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 10, min.rc = 1, max.cov.per.pool = 200, min.maf = 0.01, nlines.per.readblock = 1e+06)

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

# ================================================================================== #

# Estimate genome wide Fst 

# Use computeFST function
pooldata.fst <- computeFST(pooldata,verbose=FALSE)
pooldata.fst$Fst 

# Block-Jackknife estimation of Fst standard error and confidence intervals
pooldata.fst.bjack <- computeFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)
pooldata.fst.bjack$Fst

# Compute multi-locus Fst over sliding window of SNPs
pooldata.fst.sliding.window <- computeFST(pooldata, sliding.window.size=100)

# Plot sliding window 
pdf("output/figures/fst.sliding.window.pdf", width = 10, height = 10)
plot(pooldata.fst.sliding.window$sliding.windows.fvalues$CumMidPos/1e6, 
pooldata.fst.sliding.window$sliding.windows.fvalues$MultiLocusFst,
xlab="Cumulated Position (in Mb)", ylab="Muli-locus Fst",pch=16)
abline(h=pooldata.fst.sliding.window$Fst,lty=2) # Dashed line indicates the estimated overall genome-wide Fst
dev.off()

# ================================================================================== #

# Estimate pairwise-population Fst

# Use compute.pairwiseFST function
pooldata.pairwisefst <- compute.pairwiseFST(pooldata, verbose=FALSE)

# Graph heatmap
pdf("output/figures/heatmap.pdf", width = 10, height = 10)
heatmap(pooldata.pairwisefst)
dev.off()

# Block-Jackknife estimation of pairwise Fst standard error and confidence intervals
pooldata.pairwisefst.bjack <- compute.pairwiseFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)

# Estimated pairwise Fst are stored in the slot values: 5 first estimated pairwise
head(pooldata.pairwisefst.bjack@values)

# Graph estimated pairwise-population FST with their 95% confidence intervals 
pdf("output/figures/pairwise_Fst.pdf", width = 10, height = 10)
plot(pooldata.pairwisefst.bjack)
dev.off()

# ================================================================================== #

# Principle Components Analysis with randomallele.pca

#PCA on the read count data (the object)
pooldata.pca = randomallele.pca(pooldata, col=1:19, pch=21, main="Read Count data")

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
colors.reorder <- mycolors[c(4, 11, 5, 1, 10, 17, 8, 18, 16, 12, 13, 6, 15, 14, 3, 2, 7, 19, 9)]

# Plotting PC1 and PC2
pdf("output/figures/PCA_all_SNPs.pdf", width = 10, height = 10)
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"),
col=colors.reorder, pch=16, cex = 3, main="Read Count data")
text(pooldata.pca$pop.loadings[,1], pooldata.pca$pop.loadings[,2], pooldata@poolnames)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()
