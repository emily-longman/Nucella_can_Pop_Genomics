# Use pseudo-observed data (POD) to calibrate the XtX

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
install.packages('mvtnorm')
library(mvtnorm)

# ================================================================================== #

# Read in Baypass results

# Population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# Get omega matrix
NC.omega <- as.matrix(read.table("data/processed/GEA/baypass/omega/NC_baypass_mat_omega.out"))
# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

# Get estimates (posterior mean) for both the a_pi and b_pi parameters or the Pi Beta distribution
pi.beta.coef=read.table("data/processed/GEA/baypass/xtx/NC_baypass_core_summary_beta_params.out", h=T)$Mean

# Update the original data to obtain total allele count
NC.genobaypass=geno2YN("data/processed/GEA/baypass/genobaypass")

# Create the POD
sim.NC<-simulate.baypass(omega.mat=NC.omega, nsnp=8277206, sample.size=NC.genobaypass$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="NC.baypass.sim")

