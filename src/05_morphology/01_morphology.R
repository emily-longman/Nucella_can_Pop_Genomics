# Graph shell morphometric data.

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
install.packages(c('poolfstat', 'tidyverse', 'ggplot2', 'RColorBrewer', 'pheatmap'))
library(poolfstat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(geomorph)

# ================================================================================== #

# Read in datasets

# Read in landmark file
Nucella_landmarks <- readland.tps("data/processed/morphometrics/N_canaliculata_ventral_rand_remove_missing.tps", specID="ID")

# Read in metadata
metadata <- read.csv("data/processed/morphometrics/N_canaliculata_tps_ordered_no_missing_metadata.csv", header=T)
metadata$Site.Code <- factor(metadata$Site.Code, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

metadata <- metadata %>% mutate(genetic.structure = ifelse(Site.Code == "FC" | Site.Code == "SLR" | Site.Code ==  "SH" | Site.Code ==  "ARA" | Site.Code ==  "CBL" | Site.Code ==  "PSG" | Site.Code ==  "STC" | Site.Code ==  "KH" | Site.Code ==  "VD" | Site.Code ==  "FR" | Site.Code ==  "BMR" | Site.Code ==  "PGP", "N", "S"))
metadata$genetic.structure <- factor(metadata$genetic.structure, levels=c("N", "S"))

# Read in Procrustes pair-wise distance from MorphoJ
proc_dist <- read.csv("data/processed/morphometrics/Procrustes_dist.csv", header=T)
# Re-name rownames
row.names(proc_dist) <- proc_dist$X
proc_dist <- proc_dist[,-1]

# ================================================================================== #

# Generalized Procrustes Analysis

# Graph distribution of raw landmarks
pdf("output/figures/morphology/Nucella_landmarks.pdf", width = 8, height = 8)
plot(Nucella_landmarks)
dev.off()

# Perform a Procrustes analysis of landmark data
Nucella_gpa <- gpagen(Nucella_landmarks)

# Graph Procrustes coordinates
pdf("output/figures/morphology/Procrustes_coord.pdf", width = 8, height = 8)
plot(Nucella_gpa)
dev.off()

# ================================================================================== #

# Principal component analysis

# PCA on Procrustes coordinates
Nucella_pca <- gm.prcomp(Nucella_gpa$coords)
summary(Nucella_pca)

# Graph Eigenvectors
pdf("output/figures/morphology/PCA_scree_plot.pdf", width = 6, height = 6)
par(mar=c(5,5,4,1)+.1) # Adjust margins
plot(Nucella_pca$d, xlab="Principal Component", ylab="Eigenvalue", col="black", bg="black" , pch=21, cex=1.5, cex.lab=1.75)
dev.off()

# Graph PCA
pdf("output/figures/morphology/PCA.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
plot(Nucella_pca, xlab="PC1 (21.60%)", ylab="PC2 (19.65%)", 
xlim=c(-.1, .1), ylim=c(-.1, .1), col="black", bg=metadata$Site.Code, pch=21, cex=1, cex.lab=2)
abline(h=0, lty=2, col="grey"); abline(v=0, lty=2, col="grey")
dev.off()

# Save PC scores
pc_scores <- Nucella_pca$x

# Merge PC scores with metadata
pc_scores_metadata <- cbind(metadata, pc_scores)

# Color palette 
two_colors = c("#4575b4", "#d73027")
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph PCA with ggplot - genetic structure (PC1 vs PC2)
pdf("output/figures/morphology/PCA_ggplot_PC1_PC2.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(pc_scores_metadata, aes(Comp1, Comp2)) + 
geom_hline(yintercept=0, linetype="dashed", color ="grey") + geom_vline(xintercept=0, linetype="dashed", color ="grey") +
geom_point(aes(fill=genetic.structure), size=3, shape = 21) + 
scale_fill_manual(values=two_colors) + xlim(-.1, .1) + ylim(-.1, .1) + 
xlab("PC1 (21.60%)") + ylab("PC2 (19.65%)") + 
theme_classic(base_size = 25) + guides(fill="none")
dev.off()

# Graph PCA with ggplot - genetic structure (PC1 vs PC3)
pdf("output/figures/morphology/PCA_ggplot_PC1_PC3.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(pc_scores_metadata, aes(Comp1, Comp3)) + 
geom_hline(yintercept=0, linetype="dashed", color ="grey") + geom_vline(xintercept=0, linetype="dashed", color ="grey") +
geom_point(aes(fill=genetic.structure), size=3, shape = 21) + 
scale_fill_manual(values=two_colors) + xlim(-.1, .1) + ylim(-.1, .1) + 
xlab("PC1 (21.60%)") + ylab("PC3 (10.90%)") + 
theme_classic(base_size = 25) + guides(fill="none")
dev.off()

# Graph PCA with ggplot - populations (PC1 vs PC2)
pdf("output/figures/morphology/PCA_ggplot_PC1_PC2_sites.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(pc_scores_metadata, aes(Comp1, Comp2)) + 
geom_hline(yintercept=0, linetype="dashed", color ="grey") + geom_vline(xintercept=0, linetype="dashed", color ="grey") +
geom_point(aes(fill=Site.Code), size=3, shape = 21) + 
scale_fill_manual(values=mycolors) + xlim(-.1, .1) + ylim(-.1, .1) + 
xlab("PC1 (21.60%)") + ylab("PC2 (19.65%)") + 
theme_classic(base_size = 25) + guides(fill="none")
dev.off()

# Graph PCA with ggplot - populations (PC1 vs PC3)
pdf("output/figures/morphology/PCA_ggplot_PC1_PC3_sites.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(pc_scores_metadata, aes(Comp1, Comp3)) + 
geom_hline(yintercept=0, linetype="dashed", color ="grey") + geom_vline(xintercept=0, linetype="dashed", color ="grey") +
geom_point(aes(fill=Site.Code), size=3, shape = 21) + 
scale_fill_manual(values=mycolors) + xlim(-.1, .1) + ylim(-.1, .1) + 
xlab("PC1 (21.60%)") + ylab("PC3 (10.90%)") + 
theme_classic(base_size = 25) + guides(fill="none")
dev.off()

# ================================================================================== #


# Heatmap of Procrustes distances among the populations
pdf("output/figures/morphology/pheatmap.pdf", width = 8.5, height = 8)
pheatmap(proc_dist, border_color = "black", fontsize_col = 16, fontsize_row = 16)
dev.off()


# ================================================================================== #


# Load pooldata object 
load("data/processed/pop_structure/pooldata.RData")

#PCA on the read count data (the object)
pooldata.pca = randomallele.pca(pooldata, main="Read Count data")

# Read in metadata 
metadata <- read.csv("data/processed/pop_structure/guide_files/Populations_metadata.csv", header=T)
metadata$Population <- metadata$Site_full

# Extract PC1 
metadata$PC1_gen <- round(pooldata.pca$pop.loadings[,1],3)

# Join
pc_meta <- left_join(pc_scores, metadata, by="Population")

# Graph Genetic PC1 against morphology PC1
pdf("output/figures/morphology/Genetic_morph_PC.pdf", width = 8, height = 8)
pcs <- plot(pc_meta$PC1_gen, pc_meta$PC1, 
xlab="Genetic PC1",
ylab="Morphology PC1", pch=21, cex=1, cex.lab=1.75)
dev.off()






# ================================================================================== #
# ================================================================================== #
# Read in PC scores
pc_scores <- read.csv("data/processed/morphometrics/PC_scores.csv", header=T)

