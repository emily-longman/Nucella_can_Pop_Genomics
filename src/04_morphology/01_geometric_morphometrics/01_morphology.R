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
install.packages(c('geomorph', 'poolfstat', 'tidyverse', 'ggplot2', 'RColorBrewer', 'pheatmap'))
library(geomorph)
library(poolfstat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggbiplot)

# ================================================================================== #

# Generate Folders and files

# Make output directory
output_dir_morph="output/figures/morphology"
if (!dir.exists(output_dir_morph)) {dir.create(output_dir_morph)}

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

# Perform a Procrustes analysis of landmark data
Nucella_gpa <- gpagen(Nucella_landmarks)

# ================================================================================== #

# Graph distribution of raw landmarks
pdf("output/figures/morphology/Nucella_landmarks.pdf", width = 8, height = 8)
plot(Nucella_landmarks)
dev.off()

# Graph Procrustes coordinates
pdf("output/figures/morphology/Procrustes_coord.pdf", width = 8, height = 8)
plot(Nucella_gpa)
dev.off()

# Generate consensus shape
consensus <- apply(Nucella_gpa$coords, c(1,2), mean)
# Graph consensus shape
pdf("output/figures/morphology/Consensus.pdf", width = 6, height = 6)
plot(consensus, asp=1, type="n")
for(i in 1: length(Nucella_gpa$coords[,,3]))
points(Nucella_gpa$coords[,,i], pch=21, col="black", bg="darkgrey")
points(consensus, col="black", bg="blue", pch=21, cex=1.5)
dev.off()

# Transformation grid 
# Make reference specimen
ref <- mshape(Nucella_gpa$coords)
# Target specimen (in this case all specimens)
gp1.mn <- mshape(Nucella_gpa$coords[,,1:15])
# Graph transformation grid of target specimen to reference
pdf("output/figures/morphology/Transformation_grid.pdf", width = 6, height = 6)
plotRefToTarget(ref, gp1.mn, mag=2)
dev.off()

# If want to make a transformation grid for a specific specimen
gp.370 <- mshape(Nucella_gpa$coords[,,370])
pdf("output/figures/morphology/Transformation_grid_370.pdf", width = 6, height = 6)
plotRefToTarget(ref, gp.370, mag=2)
dev.off()

# ================================================================================== #

# Principal component analysis

# PCA on Procrustes coordinates
Nucella_pca <- gm.prcomp(Nucella_gpa$coords)
summary(Nucella_pca)

# Graph Eigenvectors
pdf("output/figures/morphology/PCA_scree_plot.pdf", width = 6, height = 6)
par(mar=c(5,5,4,1)+.1) # Adjust margins
plot(Nucella_pca$d, xlab="Principal Component", ylab="Eigenvalue", col="black", bg="black" , pch=21, cex=1.5, cex.lab=1.4)
dev.off()

# Graph PCA
pdf("output/figures/morphology/PCA.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
plot(Nucella_pca, xlab="PC1 (21.60%)", ylab="PC2 (19.65%)", 
xlim=c(-.1, .1), ylim=c(-.1, .1), col="black", bg=metadata$Site.Code, pch=21, cex=1, cex.lab=2)
abline(h=0, lty=2, col="grey"); abline(v=0, lty=2, col="grey")
dev.off()

#####

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

# Picknplot of PCA - allows you to pick a point then draw a deformation grid - need to do this on R studio
PCA <- plot(Nucella_pca)
#picknplot.shape(PCA)

# ================================================================================== #

# Heatmap of Procrustes distances among the populations
pdf("output/figures/morphology/pheatmap.pdf", width = 8.5, height = 8)
pheatmap(proc_dist, border_color = "black", fontsize_col = 16, fontsize_row = 16)
dev.off()
pdf("output/figures/morphology/pheatmap_alt.pdf", width = 8.5, height = 8)
pheatmap(proc_dist, border_color = "black", fontsize_col = 16, fontsize_row = 16, color=rev(viridis(n=361, alpha=1,begin=0, end=1, option="viridis")))
dev.off()

# ================================================================================== #

# Look at association of shape with classifiers

# Create data frame with shape data and classifiers
gdf <- geomorph.data.frame(Nucella_gpa, genetic.structure=metadata$genetic.structure, Site.Code=metadata$Site.Code, latitude=metadata$Latitude)
attributes(gdf)

# ================================================================================== #

# Perform Procrustes ANOVA with permuations 
fit.size <- procD.lm(coords ~ log(Csize), data = gdf)
summary(fit.size) 
anova(fit.size) #Significant association with shape and size

# Graph Allometry
# Predictor Line
pdf("output/figures/morphology/Allometry_PredLine.pdf", width = 8, height = 8)
plotAllometry(fit.size, size=gdf$Csize, logsz=TRUE, method="PredLine", col=gdf$Site.Code)
dev.off()
# Regression Score
pdf("output/figures/morphology/Allometry_RegScore.pdf", width = 8, height = 8)
plotAllometry(fit.size, size=gdf$Csize, logsz=TRUE, method="RegScore", col=gdf$Site.Code)
dev.off()
# Alternate way to generate predictor line with regression 
#pdf("output/figures/morphology/Allometry_PredLine.pdf", width = 8, height = 8)
#plot(fit.size, type="regression", reg.type="PredLine", predictor=log(gdf$Csize))
#dev.off()

# Partial least squares analyses to see relationship between shape and size
PLS <- two.b.pls(log(gdf$Csize), gdf$coords)

# Graph PLS
pdf("output/figures/morphology/PLS.pdf", width = 8, height = 8)
plot(PLS, col=gdf$Site.Code)
dev.off()

# Graph Common allometric component (CAC)
pdf("output/figures/morphology/Allometry_CAC.pdf", width = 8, height = 8)
plot(fit.size, size=gdf$Csize, logsz=TRUE, method="CAC", col=gdf$Site.Code)
dev.off()

# ================================================================================== #

# Perform Procrustes ANOVA analysis of Site with size
fit <- procD.lm(coords ~ log(Csize) * Site.Code, data = gdf)
# ANOVA summary
summary(fit)

# Fixed effects
fit.2 <- procD.lm(coords ~ log(Csize) + Site.Code, data = gdf)
# ANOVA summary
summary(fit.2)

# Pairwise differences between sites
PW <- pairwise(fit.size, fit, group=gdf$Site.Code)
summary(PW, test.type="dist", confidence=0.95)

# Graph Common allometric component (CAC)
pdf("output/figures/morphology/Model_CAC.pdf", width = 8, height = 8)
plot(fit, size=gdf$Csize, logsz=TRUE, method="CAC", col=gdf$Site.Code)
dev.off()

# Morphological disparity of the groups from the overall mean; Procrustes variance for each group, pairwise differences and significance of differences
MD1 <- morphol.disparity(coords~1, groups=~Site.Code, data=gdf)
# Morphological disparity while accounting for allometry
MD2 <- morphol.disparity(coords~Csize, groups=NULL, data=gdf)
# Morphological disparity of the groups from the allometric mean
MD3 <- morphol.disparity(coords~Csize, groups=~Site.Code, data=gdf)

# ================================================================================== #

# Perform Procrustes regression with permuations 
lat.regression.fit <- procD.lm(coords ~ log(Csize) * latitude, data = gdf)
anova(lat.regression.fit)

# Regression
pdf("output/figures/morphology/latitude_regression.pdf", width = 8, height = 8)
plot(lat.regression.fit)
dev.off()

# Extract C size data
regression.scores.Csize <- as.matrix(lat.regression.fit$X)

# Merge PC scores with metadata
regression_scores_Csize_metadata <- cbind(metadata, regression.scores.Csize)

# Graph lat versus log Csize 
pdf("output/figures/morphology/lat_csize.pdf", width = 9, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(regression_scores_Csize_metadata, aes(x=latitude, y=regression_scores_Csize_metadata[,20])) + 
geom_point(aes(fill=Site.Code), size=3, shape = 21) + geom_vline(xintercept=36.8007, linetype="dashed", color="black") +
scale_fill_manual(values=mycolors) + xlab("Latitude") + ylab("log(Centroid Size)") +
theme_classic(base_size = 25) + guides(fill="none")
dev.off()

# ================================================================================== #

# Extract PC1 and PC2 population data for baypass

pop.data <- pc_scores_metadata %>% 
group_by(Site.Code) %>% summarise(mean.pc1=mean(Comp1), sd.pc1=sd(Comp1), CV.pc1=sd(Comp1)/mean(Comp1), mean.pc2=mean(Comp2), sd.pc2=sd(Comp2), CV.pc2=sd(Comp2)/mean(Comp2))

write.csv(pop.data, "data/processed/morphometrics/pc.morphology.csv")

# ================================================================================== #

# Graph PC1 versus latitude
pdf("output/figures/morphology/PCA_ggplot_PC1_lat.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(pc_scores_metadata, aes(Latitude, Comp1)) + 
geom_point(aes(fill=Site.Code), size=3, shape = 21) + 
scale_fill_manual(values=mycolors) + geom_vline(xintercept=36.8007, linetype="dashed", color="black") +
xlab("Latitude") + ylab("PC1 (21.60%)") + 
theme_classic(base_size = 25) + guides(fill="none")
dev.off()

# Graph PC2 versus latitude
pdf("output/figures/morphology/PCA_ggplot_PC2_lat.pdf", width = 8, height = 8)
par(mar=c(5,5,4,1)+.1) # Adjust margins
ggplot(pc_scores_metadata, aes(Latitude, Comp2)) + 
geom_point(aes(fill=Site.Code), size=3, shape = 21) + 
scale_fill_manual(values=mycolors) + geom_vline(xintercept=36.8007, linetype="dashed", color="black") +
xlab("Latitude") + ylab("PC2 (19.65%)") + 
theme_classic(base_size = 25) + guides(fill="none")
dev.off()
