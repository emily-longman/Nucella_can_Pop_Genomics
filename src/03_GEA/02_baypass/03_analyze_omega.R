# Analyze omega matrix

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
install.packages(c('corrplot', 'ggplot2', 'RColorBrewer', 'WriteXLS'))
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(WriteXLS)

# ================================================================================== #

# Load data
NC.omega <- as.matrix(read.table("data/processed/GEA/baypass/omega/NC_baypass_mat_omega.out"))
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

# ================================================================================== #

# Create a correlation matrix of the omega values -- assess genomic differentiation between pools
cor.mat <- cov2cor(NC.omega)

pdf("output/figures/GEA/Baypass_omega_cor_matrix.pdf", width = 5, height = 5)
corrplot(cor.mat, method = "color", mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))
dev.off()
