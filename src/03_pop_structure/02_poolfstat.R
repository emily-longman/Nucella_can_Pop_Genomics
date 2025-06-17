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

# Save pooldata
save(pooldata, file="data/processed/pop_structure/pooldata.RData")
# Reload pooldata
load("data/processed/pop_structure/pooldata.RData")

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
pdf("output/figures/pop_structure/poolfstat/fst.sliding.window.pdf", width = 14, height = 8)
par(mar=c(4.5, 5, 5, 1))
plot(pooldata.fst.sliding.window$sliding.windows.fvalues$CumMidPos/1e6, 
pooldata.fst.sliding.window$sliding.windows.fvalues$MultiLocusFst,
xlab="Cumulated Position (in Mb)", ylab="Multi-locus Fst", pch=16, cex.lab=2)
abline(h=pooldata.fst.sliding.window$Fst,lty=2, lwd=4, col="red") # Dashed red line indicates the estimated overall genome-wide Fst
dev.off()

# ================================================================================== #

# Estimate pairwise-population Fst

# Use compute.pairwiseFST function
pooldata.pairwisefst <- compute.pairwiseFST(pooldata, verbose=FALSE)

# Graph heatmap
pdf("output/figures/pop_structure/poolfstat/heatmap.pdf", width = 8, height = 8)
heatmap(pooldata.pairwisefst, cexRow=1.7, cexCol=1.7)
dev.off()

# Block-Jackknife estimation of pairwise Fst standard error and confidence intervals
pooldata.pairwisefst.bjack <- compute.pairwiseFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)

# Estimated pairwise Fst are stored in the slot values: 5 first estimated pairwise
head(pooldata.pairwisefst.bjack@values)

# Graph estimated pairwise-population FST with their 95% confidence intervals 
pdf("output/figures/pop_structure/poolfstat/pairwise_Fst.pdf", width = 8, height = 8)
plot(pooldata.pairwisefst.bjack, cex=0.5)
dev.off()

# Save pairwise matrix as xls for follow-up analyses
pooldata.pairwisefst.matrix <- pooldata.pairwisefst@PairwiseFSTmatrix
pooldata.pairwisefst.matrix <- as.data.frame(pooldata.pairwisefst.matrix)
write.csv(pooldata.pairwisefst.matrix, "data/processed/pop_structure/Fst/pooldata.pairwisefst.csv")

# Graph Fst heatmap with pheatmap 
pdf("output/figures/pop_structure/poolfstat/pheatmap_Fst.pdf", width = 8, height = 8)
pheatmap(pooldata.pairwisefst.matrix, border_color = "black", fontsize_col = 16, fontsize_row = 16)
dev.off()
pdf("output/figures/pop_structure/poolfstat/pheatmap_Fst_alt.pdf", width = 8, height = 8)
pheatmap(pooldata.pairwisefst.matrix, border_color = "black", fontsize_col = 16, fontsize_row = 16, color=rev(viridis(n=361, alpha=1,begin=0, end=1, option="viridis")))
dev.off()

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
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"),
col="black", bg=colors.reorder, pch=21, cex = 3, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC3_PC4.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata.pca$pop.loadings[,3],pooldata.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata.pca$perc.var[4],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC5 and PC6
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC5_PC6.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata.pca$pop.loadings[,5],pooldata.pca$pop.loadings[,6],
xlab=paste0("PC",5," (",round(pooldata.pca$perc.var[5],2),"%)"),
ylab=paste0("PC",6," (",round(pooldata.pca$perc.var[6],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plot names on PCAs

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC1_PC2_names.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"))
text(pooldata.pca$pop.loadings[,1], pooldata.pca$pop.loadings[,2], pooldata@poolnames, cex=0.5, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/poolfstat/PCA_all_SNPs_PC3_PC4_names.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata.pca$pop.loadings[,3],pooldata.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata.pca$perc.var[4],2),"%)"))
text(pooldata.pca$pop.loadings[,3], pooldata.pca$pop.loadings[,4], pooldata@poolnames, cex=0.5, cex.lab = 1.75)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# ================================================================================== #


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


# Graph Projections of PCs

# Projection of PC 1
pdf("output/figures/pop_structure/poolfstat/Map_PC1.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC1), shape = 21, size = 5) + 
  scale_fill_gradient(low = "firebrick", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 12) + ggtitle("PC 1 Projections")
dev.off()

# Projection of PC 2
pdf("output/figures/pop_structure/poolfstat/Map_PC2.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC2), shape = 21, size = 5) + 
  scale_fill_gradient(low = "cyan1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 2 Projections")
dev.off()

# Projection of PC 3
pdf("output/figures/pop_structure/poolfstat/Map_PC3.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC3), shape = 21, size = 5) + 
  scale_fill_gradient(low = "orchid", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 3 Projections")
dev.off()

# Projection of PC 4
pdf("output/figures/pop_structure/poolfstat/Map_PC4.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC4), shape = 21, size = 5) + 
  scale_fill_gradient(low = "gold", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 4 Projections")
dev.off()

# Projection of PC 5
pdf("output/figures/pop_structure/poolfstat/Map_PC5.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC5), shape = 21, size = 5) + 
  scale_fill_gradient(low = "darkolivegreen2", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 5 Projections")
dev.off()

# Projection of PC 6
pdf("output/figures/pop_structure/poolfstat/Map_PC6.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC6), shape = 21, size = 5) + 
  scale_fill_gradient(low = "darkorange1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 6 Projections")
dev.off()


#######

# Alternative graphing of projections of PCs

# Projection of PC1
pdf("output/figures/pop_structure/poolfstat/Map_PC1_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC1), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(paste0("PC 1 Projections (",round(pooldata.pca$perc.var[1],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()
# Projection of PC2
pdf("output/figures/pop_structure/poolfstat/Map_PC2_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC2), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(paste0("PC 2 Projections (",round(pooldata.pca$perc.var[2],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()
# Projection of PC3
pdf("output/figures/pop_structure/poolfstat/Map_PC3_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC3), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(paste0("PC 3 Projections (",round(pooldata.pca$perc.var[3],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()
# Projection of PC4
pdf("output/figures/pop_structure/poolfstat/Map_PC4_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC4), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(paste0("PC 4 Projections (",round(pooldata.pca$perc.var[4],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()
# Projection of PC5
pdf("output/figures/pop_structure/poolfstat/Map_PC5_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC5), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(paste0("PC 5 Projections (",round(pooldata.pca$perc.var[5],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()
# Projection of PC6
pdf("output/figures/pop_structure/poolfstat/Map_PC6_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC6), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  ggtitle(paste0("PC 6 Projections (",round(pooldata.pca$perc.var[6],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()

# Coordinates of site labels
sites <- data.frame(
  longitude = c(-124.0593, -124.0848, -124.1148, -124.4015, -124.5647, -124.2529, -124.0809, -123.7895, -123.8036, 
                -123.2551, -123.0740, -122.3976, -121.9537, -121.9290, -121.3187, -121.2868, -120.8838, -120.6399, -120.6157),
  latitude = c(44.83777, 44.50540, 44.24999, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461, 39.28090, 38.51198, 38.31900, 
               37.18506, 36.51939, 36.44750, 35.72893, 35.66549, 35.28994, 34.88117, 34.73024),
  site.abrev = c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

lat.site.labels <- c(44.83777+0.04, 44.50540, 44.24999-0.04, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461, 39.28090, 38.51198+0.04, 38.31900-0.07, 
               37.18506, 36.51939+0.2, 36.44750-0.1, 35.72893+0.14, 35.66549-0.07, 35.28994-0.02, 34.88117, 34.73024-0.14)

long.site.labels.abrev <- c(-124.0593-0.7, -124.0848-0.8, -124.1148-0.75, -124.4015-0.8, -124.5647-0.75, -124.2529-0.8, 
                      -124.0809-0.75, -123.7895-0.7, -123.8036-0.7, -123.2551-0.75, -123.0740-0.85, -122.3976-0.8, 
                      -121.9537-0.65, -121.9290-0.8, -121.3187-0.8, -121.2868-0.75, -120.8838-0.8, -120.6399-0.8, -120.6157-0.8)

# Projection of PC1 with site codes
pdf("output/figures/pop_structure/poolfstat/Map_PC1_alt_site_labels.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC1), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  geom_text(data=sites, aes(long.site.labels.abrev, lat.site.labels, label=site.abrev)) +
  ggtitle(paste0("PC 1 Projections (",round(pooldata.pca$perc.var[1],2),"%)")) + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.9, 0.55))
dev.off()

# Projection of f3 statistics
pdf("output/figures/pop_structure/poolfstat/F3_Map.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = f3.object.sum.meta, aes(x = Long, y = Lat, fill = avg.F3), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  geom_text(data=sites, aes(long.site.labels.abrev, lat.site.labels, label=site.abrev)) +
  ggtitle("F3 Projections") + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.85, 0.55))
dev.off()

# ================================================================================== #

# Compute parameters Fs, F3, F4, D, heterozygosities

# Estimation of f-statistics on Pool-Seq data (with computation of Dstat)
pooldata.fstats <- compute.fstats(pooldata, nsnp.per.bjack.block = 1000, return.F2.blockjackknife.samples = TRUE)

# Save pooldata
save(pooldata.fstats, file="data/processed/pop_structure/pooldata.fstats.RData")
# Reload pooldata.fstats
load("data/processed/pop_structure/pooldata.fstats.RData")


# Look at stats:
head(pooldata.fstats@f2.values, 3) # 3 first f2
head(pooldata.fstats@fst.values, 3) # 3 first Fst
head(pooldata.fstats@divergence, 3) # 3 first pairwise genetic divergence
head(pooldata.fstats@f3star.values, 3) # 3 first f3*
head(pooldata.fstats@heterozygosities, 3) # 3 first D 

# Plot heterozygosities
pdf("output/figures/pop_structure/poolfstat/Heterozygosities.pdf", width = 6, height = 6)
plot(pooldata.fstats, stat.name="heterozygosities", main="Heterozygosities")
dev.off()

# Save heterozygosities as a data frame 
pooldata.fstats@heterozygosities %>% as.data.frame -> het.object
# Rename columns
het.object %>% rename(bjack_mean='bjack mean', bjack_s.e.='bjack s.e.') -> het.object
# Order sites
het.object$Site <- row.names(het.object)
het.object$Site <- factor(het.object$Site, 
levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD","KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))


# Graph heterozyogisities with ggplot
pdf("output/figures/pop_structure/poolfstat/heterozygosities_mean_se_plot.pdf", width = 8, height = 10)
ggplot(data = het.object, aes(x=bjack_mean, y=Site)) + geom_point(size=2) + 
geom_errorbar(aes(xmin=bjack_mean-bjack_s.e., xmax=bjack_mean+bjack_s.e.), width=0.3) + xlab("Heterozygosity") +
theme_classic(base_size = 30) 
dev.off()

# ================================================================================== #

# Graph F3 

# Plot F3
pdf("output/figures/pop_structure/poolfstat/F3.pdf", width = 6, height = 6)
plot(pooldata.fstats, stat.name="F3", main="F3")
dev.off()
# A negative f3 statistic is evidence that the target population is admixed. Almost all are positive. 
# Graph only negative F3
pdf("output/figures/pop_structure/poolfstat/F3_only_neg.pdf", width = 6, height = 6)
plot(pooldata.fstats, stat.name="F3", value.range=c(NA, 0), main="F3 (only populations with negative f3 statistics)")
dev.off()

# Save F3 values as a data frame
pooldata.fstats@f3.values %>% as.data.frame -> f3.object
# Rename columns
f3.object %>% rename(bjack_mean='bjack mean', bjack_s.e.='bjack s.e.') -> f3.object

# Split row name up into focal pop and the two parent/source 
f3.object %>% mutate(pair_info = row.names(.)) %>% 
separate(pair_info, into = c("Site","Parents"), sep = ";") %>% 
separate(Parents, into = c("P1","P2"), sep = ",") -> f3.object.flt

# Order sites
f3.object.flt$Site <- factor(f3.object.flt$Site, 
levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph boxplot of f3 values
pdf("output/figures/pop_structure/poolfstat/F3_box_plot.pdf", width = 12, height = 6)
ggplot(data = f3.object.flt, aes(x=Site, y =`bjack mean`)) + geom_boxplot() + theme_classic(base_size = 20) 
dev.off()

# Summarize F3 for each Site
f3.object.flt %>% group_by(Site) %>% dplyr::summarize(avg.F3=mean(bjack_mean)) -> f3.object.sum

# Combine metadata with F3 statistics
f3.object.sum.meta <- left_join(f3.object.sum, metadata, by="Site")

# Coordinates of site labels
sites <- data.frame(
  longitude = c(-124.0593, -124.0848, -124.1148, -124.4015, -124.5647, -124.2529, -124.0809, -123.7895, -123.8036, 
                -123.2551, -123.0740, -122.3976, -121.9537, -121.9290, -121.3187, -121.2868, -120.8838, -120.6399, -120.6157),
  latitude = c(44.83777, 44.50540, 44.24999, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461, 39.28090, 38.51198, 38.31900, 
               37.18506, 36.51939, 36.44750, 35.72893, 35.66549, 35.28994, 34.88117, 34.73024),
  site.abrev = c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

lat.site.labels <- c(44.83777+0.04, 44.50540, 44.24999-0.04, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461, 39.28090, 38.51198+0.04, 38.31900-0.07, 
               37.18506, 36.51939+0.2, 36.44750-0.1, 35.72893+0.14, 35.66549-0.07, 35.28994-0.02, 34.88117, 34.73024-0.14)

long.site.labels.abrev <- c(-124.0593-0.7, -124.0848-0.8, -124.1148-0.75, -124.4015-0.8, -124.5647-0.75, -124.2529-0.8, 
                      -124.0809-0.75, -123.7895-0.7, -123.8036-0.7, -123.2551-0.75, -123.0740-0.85, -122.3976-0.8, 
                      -121.9537-0.65, -121.9290-0.8, -121.3187-0.8, -121.2868-0.75, -120.8838-0.8, -120.6399-0.8, -120.6157-0.8)

# Projection of f3 statistics
pdf("output/figures/pop_structure/poolfstat/F3_Map.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = f3.object.sum.meta, aes(x = Long, y = Lat, fill = avg.F3), shape = 21, size = 6) + 
  scale_fill_gradient(low = "cyan", high = "black") + 
             coord_fixed(1.3) + theme_classic(base_size = 20) +
  xlim(c(-125.5, -114))  +
  xlab("Longitude") + ylab("Latitude") + 
  geom_text(data=sites, aes(long.site.labels.abrev, lat.site.labels, label=site.abrev)) +
  ggtitle("F3 Projections") + 
  theme(plot.title=element_text(family='', face='bold', size=25)) +
  theme(legend.position =  c(0.85, 0.55))
dev.off()

# ================================================================================== #

# Test for admixture (note: a Z-score < −1.65 provides evidence for admixture at the 95% significance threshold)
# Z-score = ratio of the block-jackknife estimated mean and standard-error
tst.sel <- pooldata.fstats@f3.values$`Z-score`< -1.65

pooldata.fstats@f3.values[tst.sel,] -> test.sel.df
test.sel.df
#                 Estimate    bjack mean   bjack s.e.    Z-score
#CBL;ARA,PSG -0.0006874833 -0.0006756751 3.032273e-05 -22.282791
#CBL;ARA,STC -0.0004121394 -0.0003928102 3.317190e-05 -11.841655
#SBR;BMR,PB  -0.0012524884 -0.0013422484 6.276391e-04  -2.138567
#SBR;BMR,PSN -0.0017790245 -0.0016602182 6.470612e-04  -2.565782
#SBR;FR,PB   -0.0013355210 -0.0013971432 6.297338e-04  -2.218625
#SBR;FR,PSN  -0.0018677455 -0.0017165174 6.497936e-04  -2.641635
#SBR;KH,PB   -0.0016284524 -0.0017639007 6.320927e-04  -2.790573
#SBR;KH,PSN  -0.0021293997 -0.0020184018 6.516912e-04  -3.097175
#SBR;VD,PB   -0.0019272877 -0.0021025096 6.302147e-04  -3.336180
#SBR;VD,PSN  -0.0024482292 -0.0023965779 6.500269e-04  -3.686890
#SBR;PB,PGP  -0.0012845715 -0.0014763223 6.287333e-04  -2.348090
#SBR;PGP,PSN -0.0018096816 -0.0017313944 6.502507e-04  -2.662657

# Evidence of an admixed origin for CBL with ancestral sources related to ARA, PSG and STC. 
# Evidence of an admixed origin for SBR with populations from both north and south of that site.

# Rename columns
test.sel.df %>% rename(bjack_mean='bjack mean', bjack_s.e.='bjack s.e.') -> test.sel.df

# Split row name up into focal pop and the two parent/source 
test.sel.df %>% mutate(pair_info = row.names(.)) %>% 
separate(pair_info, into = c("Site","Parents"), sep = ";") %>% 
separate(Parents, into = c("P1","P2"), sep = ",") -> test.sel.df.flt

# Graph
pdf("output/figures/pop_structure/poolfstat/F3_test_sel.pdf", width = 12, height = 6)
ggplot(data = test.sel.df.flt, aes(y=rownames(test.sel.df.flt), x=bjack_mean)) + geom_point() + 
geom_errorbar(aes(y=rownames(test.sel.df.flt), xmin=bjack_mean-bjack_s.e., xmax=bjack_mean+bjack_s.e.), width=0.4, position = position_dodge(.9)) +
geom_vline(xintercept=0, color="red", linetype="dashed", linewidth=1.5) + xlim(NA,0.0001) + xlab("F3 Statistic") + ylab("") +
theme_classic(base_size = 20) 
dev.off()
