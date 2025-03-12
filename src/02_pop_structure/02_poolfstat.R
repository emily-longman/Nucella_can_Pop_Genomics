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
install.packages(c('poolfstat', 'ggplot2', 'RColorBrewer', 'WriteXLS'))
library(poolfstat)
library(ggplot2)
library(RColorBrewer)
library(WriteXLS)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh - # variants: 26,206,958)
#pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf", 
#poolsizes=rep(40,19), poolnames=pops$V1, 
#min.cov.per.pool = 15, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 11,656,080 SNPs for 19 Pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh - # variants: 14,897,468)
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 8,277,206 SNPs for 19 Pools

#pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ50_maxmissing0.9.recode.vcf", 
#poolsizes=rep(40,19), poolnames=pops$V1, 
#min.cov.per.pool = 15, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 9,804,957 SNPs for 19 Pools

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

# ================================================================================== #

# Estimate genome wide Fst 

# Use computeFST function (default to using the Anova method)
pooldata.fst <- computeFST(pooldata,verbose=FALSE)
pooldata.fst$Fst 
# 0.5801261

# Block-Jackknife estimation of Fst standard error and confidence intervals
pooldata.fst.bjack <- computeFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)
pooldata.fst.bjack$Fst
#   Estimate  bjack mean  bjack s.e.     CI95inf     CI95sup 
# 0.580126095 0.582574371 0.001481405 0.579670817 0.585477924 

# Compute multi-locus Fst over sliding window of SNPs
pooldata.fst.sliding.window <- computeFST(pooldata, sliding.window.size=1000)

# Plot sliding window
pdf("output/figures/pop_structure/poolfstat/fst.sliding.window.pdf", width = 8, height = 8)
plot(pooldata.fst.sliding.window$sliding.windows.fvalues$CumMidPos/1e6, 
pooldata.fst.sliding.window$sliding.windows.fvalues$MultiLocusFst,
xlab="Cumulated Position (in Mb)", ylab="Multi-locus Fst", pch=16)
abline(h=pooldata.fst.sliding.window$Fst,lty=2, col="red") # Dashed red line indicates the estimated overall genome-wide Fst
dev.off()

# ================================================================================== #

# Estimate pairwise-population Fst

# Use compute.pairwiseFST function
pooldata.pairwisefst <- compute.pairwiseFST(pooldata, verbose=FALSE)

# Graph heatmap
pdf("output/figures/pop_structure/poolfstat/heatmap.pdf", width = 8, height = 8)
heatmap(pooldata.pairwisefst)
dev.off()

# Block-Jackknife estimation of pairwise Fst standard error and confidence intervals
pooldata.pairwisefst.bjack <- compute.pairwiseFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)

# Estimated pairwise Fst are stored in the slot values: 5 first estimated pairwise
head(pooldata.pairwisefst.bjack@values)

# Graph estimated pairwise-population FST with their 95% confidence intervals 
pdf("output/figures/pop_structure/poolfstat/pairwise_Fst.pdf", width = 8, height = 8)
plot(pooldata.pairwisefst.bjack, cex=0.5)
dev.off()

# Save pairwise matrix 
pooldata.pairwisefst.matrix <- pooldata.pairwisefst@PairwiseFSTmatrix
pooldata.pairwisefst.matrix <- as.data.frame(pooldata.pairwisefst.matrix)
WriteXLS(pooldata.pairwisefst.matrix, "data/processed/pop_structure/Fst/pooldata.pairwisefst.xls")

# ================================================================================== #

# Principle Components Analysis with randomallele.pca

#PCA on the read count data (the object)
pooldata.pca = randomallele.pca(pooldata, main="Read Count data")

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
colors.reorder <- mycolors[c(19,2,3,4,11,5,1,10,17,14,6,8,7,13,9,18,16,12,15)]

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC1_PC2.pdf", width = 8, height = 8)
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"),
col="black", bg=colors.reorder, pch=21, cex = 3, main="Read Count data")
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC3_PC4.pdf", width = 8, height = 8)
pca <- plot(pooldata.pca$pop.loadings[,3],pooldata.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata.pca$perc.var[4],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, main="Read Count data")
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC5 and PC6
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC5_PC6.pdf", width = 8, height = 8)
pca <- plot(pooldata.pca$pop.loadings[,5],pooldata.pca$pop.loadings[,6],
xlab=paste0("PC",5," (",round(pooldata.pca$perc.var[5],2),"%)"),
ylab=paste0("PC",6," (",round(pooldata.pca$perc.var[6],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, main="Read Count data")
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plot names on PCAs

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC1_PC2_names.pdf", width = 8, height = 8)
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"))
text(pooldata.pca$pop.loadings[,1], pooldata.pca$pop.loadings[,2], pooldata@poolnames, cex=0.5)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC3_PC4_names.pdf", width = 8, height = 8)
pca <- plot(pooldata.pca$pop.loadings[,3],pooldata.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata.pca$perc.var[4],2),"%)"))
text(pooldata.pca$pop.loadings[,3], pooldata.pca$pop.loadings[,4], pooldata@poolnames, cex=0.5)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# ================================================================================== #

library(maps) 
library(mapdata)
library(ggrepel)

# Read in metadata 
metadata <- read.csv("data/processed/pop_structure/guide_files/Populations_metadata.csv", header=T)

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Combine PC loadings with metadata
metadata$PC1 <- round(pooldata.pca$pop.loadings[,1],3)
metadata$PC2 <- round(pooldata.pca$pop.loadings[,2],3)
metadata$PC3 <- round(pooldata.pca$pop.loadings[,3],3)
metadata$PC4 <- round(pooldata.pca$pop.loadings[,4],3)
metadata$PC5 <- round(pooldata.pca$pop.loadings[,5],3)
metadata$PC6 <- round(pooldata.pca$pop.loadings[,6],3)


pdf("output/figures/pop_structure/poolfstat/Map_PC1.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC1), shape = 21, size = 5) + 
  scale_fill_gradient(low = "firebrick", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 1 projections onto Map")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC2.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC2), shape = 21, size = 5) + 
  scale_fill_gradient(low = "cyan1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 2 projections onto Map")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC3.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC3), shape = 21, size = 5) + 
  scale_fill_gradient(low = "orchid", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 3 projections onto Map")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC4.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC4), shape = 21, size = 5) + 
  scale_fill_gradient(low = "gold", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 4 projections onto Map")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC5.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC5), shape = 21, size = 5) + 
  scale_fill_gradient(low = "darkolivegreen2", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 5 projections onto Map")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC6.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC6), shape = 21, size = 5) + 
  scale_fill_gradient(low = "darkorange1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 6 projections onto Map")
dev.off()


save.image("data/processed/pop_structure/N.can_poolfstat.RData")

# ================================================================================== #

# Compute parameters Fs, F3, F4, D, heterozygosities

# Estimation of f-statistics on Pool-Seq data (with computation of Dstat)
pooldata.fstats <- compute.fstats(pooldata, nsnp.per.bjack.block = 1000, computeDstat = TRUE)

# Look at stats:
head(pooldata.fstats@f2.values, 3) # 3 first f2
head(pooldata.fstats@fst.values, 3) # 3 first Fst
head(pooldata.fstats@divergence, 3) # 3 first pairwise genetic divergence
head(pooldata.fstats@f3star.values, 3) # 3 first f3*
head(pooldata.fstats@f4values, 3) # 3 first f4 
head(pooldata.fstats@Dstat.values, 3) # 3 first D 

# Plot heterozygosities
pdf("output/figures/pop_structure/poolfstat/Heterozygosities.pdf", width = 6, height = 6)
plot(pooldata.fstats, stat.name="heterozygosities", main="Heterozygosities")
dev.off()

# ================================================================================== #

# Save heterozygosities data 
pooldata.fstats.f3.matrix <- pooldata.fstats@f3star.values
pooldata.fstats.f3.matrix <- as.data.frame(pooldata.fstats.f3.matrix)
WriteXLS(pooldata.fstats.f3.matrix, "data/processed/pop_structure/Fst/pooldata.fstats.f3_filt.xls")


# Save heterozygosities data 
pooldata.fstats.het.matrix <- pooldata.fstats@heterozygosities
pooldata.fstats.het.matrix <- as.data.frame(pooldata.fstats.het.matrix)
WriteXLS(pooldata.fstats.het.matrix, "data/processed/pop_structure/Fst/pooldata.fstats.het_filt.xls")

