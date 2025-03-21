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

# Save pairwise matrix as xls for follow-up analyses
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
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 1 Projections")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC2.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC2), shape = 21, size = 5) + 
  scale_fill_gradient(low = "cyan1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 2 Projections")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC3.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC3), shape = 21, size = 5) + 
  scale_fill_gradient(low = "orchid", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 3 Projections")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC4.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC4), shape = 21, size = 5) + 
  scale_fill_gradient(low = "gold", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 4 Projections")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC5.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC5), shape = 21, size = 5) + 
  scale_fill_gradient(low = "darkolivegreen2", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 5 Projections")
dev.off()

pdf("output/figures/pop_structure/poolfstat/Map_PC6.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = metadata, aes(x = Long, y = Lat, fill = PC6), shape = 21, size = 5) + 
  scale_fill_gradient(low = "darkorange1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() + ggtitle("PC 6 Projections")
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

# Plot F3
pdf("output/figures/pop_structure/poolfstat/F3.pdf", width = 6, height = 6)
plot(pooldata.fstats, stat.name="F3", main="F3")
dev.off()
# A negative f3 statistic is evidence that the target population is admixed. Almost all are positive. 
# Graph only negative F3
pdf("output/figures/pop_structure/poolfstat/F3_only_neg.pdf", width = 6, height = 6)
plot(pooldata.fstats, stat.name="F3", value.range=c(NA, 0), main="F3 (only populations with negative f3 statistics)")
dev.off()

# ================================================================================== #

# Test for admixture (note: a Z-score < −1.65 provides evidence for admixture)
tst.sel <- pooldata.fstats@f3.values$`Z-score`< -1.65
pooldata.fstats@f3.values[tst.sel,]
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

# Test for four-population treeness (note:  a Z-score lower than 1.96 in absolute value provides no evidence against the null-hypothesis of treeness for the tested population configuration at the 95% significance threshold)
tst.sel<-abs(pooldata.fstats@f4.values$`Z-score`)<1.96
as.data.frame(pooldata.fstats@f4.values)[tst.sel,]

# Estimating admixture proportions with f4-ratios
# Example of estimating f4 ratios
compute.f4ratio(pooldata.fstats, num.quadruplet = "HZD,PSN;PB,OCT", den.quadruplet="HZD,PSN;SBR,STR")
#  Estimate bjack mean bjack s.e.    CI95inf    CI95sup 
# 3.0296472  2.9891599  0.0992849  2.7945615  3.1837583 
# Note: Simulated value (α = 0.25) is within the confidence interval

# ================================================================================== #

# Build admixture graph from scratch

# The find.tree.popset function selects maximal sets of unadmixed populations from an fstats object
scaf.pops <- find.tree.popset(pooldata.fstats, verbose=FALSE)

# List the 15 passing the treeness test for the identified set
scaf.pops$passing.quaduplet

# State the range of variation of the passing quadruplet
scaf.pops$Z_f4.range

# The rooted.njtree.builder function builds the scaffold tree based on the set of unadmixed populations (identified above with find.tree.popset) 
scaf.tree <- rooted.njtree.builder(fstats=pooldata.fstats, pop.sel=scaf.pops$pop.sets[1,], plot.nj=FALSE)
# Score of the NJ tree: 27.19838 

# Check the fit of the neighbor-joining tree
scaf.tree$nj.tree.eval

# Note: need to graph using R GUI

# Plot best inferred rooted scaffold tree
plot(scaf.tree$best.rooted.tree)
# This tree can be used as a reference graph to construct the complete admixture graph

# Add other populations
add.pool.OCT <- add.leaf(scaf.tree$best.rooted.tree,leaf.to.add="OCT", 
                         fstats=pooldata.fstats, verbose=FALSE, drift.scaling=TRUE)

plot(add.pool.OCT$best.fitted.graph)

add.pool.OCT.HZD <- add.leaf(scaf.tree$best.rooted.tree,leaf.to.add="HZD", 
                         fstats=pooldata.fstats, verbose=FALSE, drift.scaling=TRUE)

plot(add.pool.OCT.HZD$best.fitted.graph)

add.pool.OCT.HZD.SBR <- add.leaf(scaf.tree$best.rooted.tree,leaf.to.add="SBR", 
                             fstats=pooldata.fstats, verbose=FALSE, drift.scaling=TRUE)

plot(add.pool.OCT.HZD.SBR$best.fitted.graph)