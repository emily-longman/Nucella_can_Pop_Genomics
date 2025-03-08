# Graph Baypass XtX output.

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
install.packages(c('WriteXLS'))
library(WriteXLS)

# ================================================================================== #

# Read in Baypass results
NC.omega <- as.matrix(read.table("data/processed/GEA/baypass/omega/NC_baypass_mat_omega.out"))
XtX <- read.table("data/processed/GEA/baypass/xtx/NC_baypass_core_summary_pi_xtx.out", header=T)
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

# ================================================================================== #

# Check the behavior of the p-values associated to the XtXst estimator
pdf("output/figures/GEA/Baypass_xtx_hist.pdf", width = 5, height = 5)
hist(10**(-1*XtX$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/GEA/Baypass_xtx.pdf", width = 5, height = 5)
plot(XtX$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/GEA/Baypass_xtx_outliers.pdf", width = 5, height = 5)
plot(XtX$log10.1.pval., ylab="XtX P-value (-log10 scale)" )
abline(h=3, lty=2) #0.001 p-value threshold
dev.off()

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Use pseudo-observed data (POD) to calibrate the XtX

# Load packages
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages('mvtnorm')
library(mvtnorm)


# Get estimates (posterior mean) for both the a_pi and b_pi parameters or the Pi Beta distribution
pi.beta.coef=read.table("data/processed/GEA/baypass/xtx/NC_baypass_core_summary_beta_params.out", h=T)$Mean

# Update the original data to obtain total allele count
NC.genobaypass=geno2YN("data/processed/GEA/baypass/genobaypass")

# Create the POD
simu.NC<-simulate.baypass(omega.mat=NC.omega, nsnp=8277206, sample.size=NC.genobaypass$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="data/processed/GEA/baypass/NC.baypass.sim")

