# Use poolfstat to convert VCF to Baypass input file

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
install.packages(c('poolfstat'))
library(poolfstat)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/pop_structure/guide_files/Nucella_pops.list", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf.gz", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 2, max.cov.per.pool = 200, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 13,350,988 SNPs for 19 Pools

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

# ================================================================================== #

# Convert to BayPass input file 
pooldata2genobaypass(pooldata, writing.dir = "data/processed/GEA/baypass", subsamplesize = -1)
# Three output files = genobaypass (allele counts), poolsize (haploid size per pool), & snpdet (snp info matrix). 
# Subsample size can be used to sample to a smaller number of SNPs. If the subsample size is <0, then all SNPs are included in the BayPass files.