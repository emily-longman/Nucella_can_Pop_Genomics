# Create GDS object from the vcf 

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
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Create a gds object from the vcf

# Load the vcf file (needs to be a .gz file)
vcf.fn="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf"

# Parse the header
#seqVCF_Header(vcf.fn)

# ================================================================================== #

gds.fn=gsub(".vcf", ".gds", vcf.fn)
vcf.fn=paste(vcf.fn, ".gz", sep="")

# ================================================================================== #

# Convert VCF to GDS
# NOTE: make sure you make a directory called gds first
seqVCF2GDS(vcf.fn, storage.option="ZIP_RA", "data/processed/fastq_to_vcf/gds/N.canaliculata_pops.gds")

