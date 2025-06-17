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

# Read in Baypass results
baypass_morph_PC <- read.table("data/processed/morphometrics/baypass/PC/NC_baypass_morphology_PC_summary_pi_xtx.out", header=T)
baypass_morph_CV <- read.table("data/processed/morphometrics/baypass/CV/NC_baypass_morphology_CV_summary_pi_xtx.out", header=T)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
 
# ================================================================================== #
# ================================================================================== #

# Analyze and Graph PC morphology baypass results

# Look at the behavior of the p-values associated to the XtX estimator
pdf("output/figures/morphology/baypass/Baypass_xtx_hist_PC.pdf", width = 5, height = 5)
hist(10**(-1*baypass_morph_PC$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/morphology/baypass/Baypass_xtx_PC.pdf", width = 10, height = 5)
plot(baypass_morph_PC$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_PC.pdf", width = 10, height = 5)\
par(mar=c(5,5,4,1)+.1) # Adjust margins
plot(baypass_morph_PC$log10.1.pval., ylab="-log10(P-value)", xlab="Position")
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_PC)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()

# ================================================================================== #

# Fancier graphing with ggplot

# Join baypass morphology PC results with snp metadata
baypass_morph_PC_pos <- cbind(snp.meta, baypass_morph_PC)

# Graph with ggplot
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_PC_ggplot.pdf", width = 10, height = 5)
ggplot(baypass_morph_PC_pos, aes(x = chr, y = log10.1.pval.)) +
geom_jitter(size = 1, alpha = 0.7, show.legend = FALSE, height=NULL) + 
geom_hline(yintercept = -log10(0.05/dim(baypass_morph_PC)[1]), linetype = "dashed", color = "red", linewidth = 0.8)  +  
theme_classic(base_size = 20) +  
theme(panel.spacing = unit(0.5, "lines"),
axis.text.x = element_blank(), axis.ticks.x = element_blank(),
axis.text.y = element_text(size = 12)) +
labs(x = "Position", y = expression(-log[10](italic(p))))
dev.off()

# ================================================================================== #

# Identify SNP list for PC morphology baypass results

# Identify bonferroni significant SNPs -- 10 SNPs
baypass_morph_PC_pos_bonf_sig_SNPs <- baypass_morph_PC_pos[which(baypass_morph_PC_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_PC_pos)[1])),]

# Write file of bonferroni outlier SNPs
write.csv(baypass_morph_PC_pos_bonf_sig_SNPs, "data/processed/morphometrics/baypass/baypass_morph_PC_pos_bonf_sig_SNPs")

# ================================================================================== #
# ================================================================================== #

# Analyze and Graph CV morphology baypass results

# Look at the behavior of the p-values associated to the XtX estimator
pdf("output/figures/morphology/baypass/Baypass_xtx_hist_CV.pdf", width = 5, height = 5)
hist(10**(-1*baypass_morph_CV$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/morphology/baypass/Baypass_xtx_CV.pdf", width = 10, height = 5)
plot(baypass_morph_CV$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()

# ================================================================================== #

# Fancier graphing with ggplot

# Join baypass morphology PC results with snp metadata
baypass_morph_CV_pos <- cbind(snp.meta, baypass_morph_CV)

# Graph with ggplot
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_ggplot.pdf", width = 10, height = 5)
ggplot(baypass_morph_CV_pos, aes(x = chr, y = log10.1.pval.)) +
geom_jitter(size = 1, alpha = 0.7, show.legend = FALSE, height=NULL) + 
geom_hline(yintercept = -log10(0.05/dim(baypass_morph_CV)[1]), linetype = "dashed", color = "red", linewidth = 0.8)  +  
theme_classic(base_size = 20) +  
theme(panel.spacing = unit(0.5, "lines"),
axis.text.x = element_blank(), axis.ticks.x = element_blank(),
axis.text.y = element_text(size = 12)) +
labs(x = "Position", y = expression(-log[10](italic(p))))
dev.off()

# ================================================================================== #

# Identify SNP list for CV morphology baypass results

# Identify bonferroni significant SNPs -- 6 SNPS
baypass_morph_CV_pos_bonf_sig_SNPs <- baypass_morph_CV_pos[which(baypass_morph_CV_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV_pos)[1])),]

# Write file of bonferroni outlier SNPs
write.csv(baypass_morph_CV_pos_bonf_sig_SNPs, "data/processed/morphometrics/baypass/baypass_morph_CV_pos_bonf_sig_SNPs")
