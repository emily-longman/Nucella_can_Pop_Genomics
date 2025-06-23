# Filter pooldata object for only the SNPs within the protein of interest

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
install.packages(c('poolfstat', 'tidyverse'))
library(poolfstat)
library(tidyverse)

# ================================================================================== #

# Read in population names
pops <- read.table("guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
# Load pooldata
load("data/processed/pop_structure/pooldata.RData")

# ================================================================================== #


# Get the snp.info matrix
snp_info <- pooldata@snp.info

# Filter for a specific chromosome (e.g., "chr1")
selected_snps_ntLink_3633 <- which(snp_info$Chromosome == "ntLink_3633")
selected_snps_prot_g26813 <- which(snp_info$Chromosome == "ntLink_3633" & snp_info$Position >= 35566 & snp_info$Position <= 48839)

# Subset the pooldata object using the selected SNP indices
pooldata_ntLink_3633 <- pooldata.subset(pooldata, snp.index = selected_snps_ntLink_3633)
pooldata_prot_g26813 <- pooldata.subset(pooldata, snp.index = selected_snps_prot_g26813)

# ================================================================================== #

# Save pooldata
save(pooldata_ntLink_3633, file="data/processed/pop_structure/pooldata_subset_ntLink_3633.RData")
save(pooldata_prot_g26813, file="data/processed/pop_structure/pooldata_subset_prot_g26813.RData")