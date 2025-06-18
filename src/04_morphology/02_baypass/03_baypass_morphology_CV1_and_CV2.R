# Graph Baypass morphology output.

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
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages(c('data.table', 'dplyr', 'ggplot2', 'mvtnorm', 'geigen'))
library(data.table)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(geigen)

# ================================================================================== #

# Generate Folders and files

# Make output directory
output_dir_morph_baypass="output/figures/morphology/baypass"
if (!dir.exists(output_dir_morph_baypass)) {dir.create(output_dir_morph_baypass)}

# ================================================================================== #

# Read in Baypass results 
baypass_morph_CV1 <- read.table("data/processed/morphometrics/baypass/CV_1/NC_baypass_morphology_CV_1_summary_pi_xtx.out", header=T)
baypass_morph_CV2 <- read.table("data/processed/morphometrics/baypass/CV_2/NC_baypass_morphology_CV_2_summary_pi_xtx.out", header=T)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Join baypass morphology CV results with snp metadata
baypass_morph_CV1_pos <- cbind(snp.meta, baypass_morph_CV1)
baypass_morph_CV2_pos <- cbind(snp.meta, baypass_morph_CV2)

# ================================================================================== #

# Graph XtX

# Graph XtX for CV1
pdf("output/figures/morphology/baypass/Baypass_xtx_CV_1.pdf", width = 10, height = 5)
plot(baypass_morph_CV1$XtXst)
dev.off()

# Graph XtX for CV2
pdf("output/figures/morphology/baypass/Baypass_xtx_CV_2.pdf", width = 10, height = 5)
plot(baypass_morph_CV2$XtXst)
dev.off()

# ================================================================================== #

# Identify bonferroni significant SNPs 
baypass_morph_CV1_pos_bonf_sig_SNPs <- baypass_morph_CV1_pos[which(baypass_morph_CV1_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV1_pos)[1])),]
baypass_morph_CV2_pos_bonf_sig_SNPs <- baypass_morph_CV2_pos[which(baypass_morph_CV2_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV2_pos)[1])),]

# Write file of bonferroni outlier SNPs
write.csv(baypass_morph_CV1_pos_bonf_sig_SNPs, "data/processed/morphometrics/baypass/CV_1/baypass_morph_CV1_pos_bonf_sig_SNPs.csv", row.names=FALSE)
write.csv(baypass_morph_CV2_pos_bonf_sig_SNPs, "data/processed/morphometrics/baypass/CV_2/baypass_morph_CV2_pos_bonf_sig_SNPs.csv", row.names=FALSE)

# ================================================================================== #

# Graph outliers for CV1
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_1.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV1$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV1)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()

# Graph outliers for CV2
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_2.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV2$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV2)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()

# ================================================================================== #

# Highlight outliers in graphs

# Set colors for outliers
baypass_morph_CV1_pos$color <- ifelse(baypass_morph_CV1_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV1_pos)[1]), "sig", "not.sig")
baypass_morph_CV2_pos$color <- ifelse(baypass_morph_CV2_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV2_pos)[1]), "sig", "not.sig")

# CV1
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_1_alt.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV1$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position", 
col=ifelse(baypass_morph_CV1_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV1_pos)[1]), "black", "#bebebe93"), pch=19)
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV1_pos)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()

# CV2
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_2_alt.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV2_pos$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position", 
col=ifelse(baypass_morph_CV2_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV2_pos)[1]), "black", "#bebebe93"), pch=19)
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV2_pos)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()


# graph with ggplot
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_1_ggplot.pdf", width = 10, height = 5)
ggplot(baypass_morph_CV1_pos, aes(x = chr, y = log10.1.pval., color=color)) +
geom_jitter(size = 1.2, alpha = 0.7, show.legend = FALSE, height=NULL) + 
scale_color_manual(values=c("#bebebecc", "black")) +
geom_hline(yintercept = -log10(0.05/dim(baypass_morph_CV1_pos)[1]), linetype = "dashed", color = "red", linewidth = 0.8)  +
theme_classic(base_size = 20) +  
theme(panel.spacing = unit(0.5, "lines"),
axis.text.x = element_blank(), axis.ticks.x = element_blank(),
axis.text.y = element_text(size = 12)) +
labs(x = "Position", y = expression(-log[10](italic(p))))
dev.off()