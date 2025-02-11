# Generate PCA of all SNPs

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move gwas output to results dir 

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Set relative path of results directory from root
dir(find_root_file("data", "processed",  criterion = has_file("README.md")))
data_path_from_root <- find_root_file("data", "processed", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(data_path_from_root)

# Set working directory as path from root
setwd(data_path_from_root)

# ================================================================================== #

# Load packages
install.packages(c('poolfstat'))
require(poolfstat)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools )
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools
pooldata <-vcf2pooldata(vcf.file="fastq_to_vcf/vcf_freebayes/N_can_pops.vcf",poolsizes=rep(40,19))
pooldata <-vcf2pooldata(vcf.file="fastq_to_vcf/vcf_clean/N_can_pops_output_snps-only.vcf.recode.vcf",poolsizes=rep(40,19))

sim6p.readcount30X <-vcf2pooldata(vcf.file="sim6p.poolseq30X.vcf.gz",poolsizes=rep(50,6))