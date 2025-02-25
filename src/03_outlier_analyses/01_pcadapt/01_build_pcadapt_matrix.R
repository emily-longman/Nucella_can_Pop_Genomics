# Build pcadapt matrix

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
install.packages(c('poolfstat', 'WriteXLS'))
library(poolfstat)
library(WriteXLS)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh - # variants: 26,206,958)
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 15, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 11,656,080 SNPs for 19 Pools

# ================================================================================== #

ref_STR <- pooldata@refallele.readcount[,1]
ref_SLR <- pooldata@refallele.readcount[,2]
ref_SH <- pooldata@refallele.readcount[,3]
ref_ARA <- pooldata@refallele.readcount[,4]
ref_BMR <- pooldata@refallele.readcount[,5]
ref_CBL <- pooldata@refallele.readcount[,6]
ref_FC <- pooldata@refallele.readcount[,7]
ref_FR <- pooldata@refallele.readcount[,8]
ref_HZD <- pooldata@refallele.readcount[,9]
ref_SBR <- pooldata@refallele.readcount[,10]
ref_PSG <- pooldata@refallele.readcount[,11]
ref_KH <- pooldata@refallele.readcount[,12]
ref_STC <- pooldata@refallele.readcount[,13]
ref_PL <- pooldata@refallele.readcount[,14]
ref_VD <- pooldata@refallele.readcount[,15]
ref_OCT <- pooldata@refallele.readcount[,16]
ref_PB <- pooldata@refallele.readcount[,17]
ref_PGP <- pooldata@refallele.readcount[,18]
ref_PSN <- pooldata@refallele.readcount[,19]

alt_STR <- pooldata@readcoverage[,1] - pooldata@refallele.readcount[,1]
alt_SLR <- pooldata@readcoverage[,2] - pooldata@refallele.readcount[,2]
alt_SH <- pooldata@readcoverage[,3] - pooldata@refallele.readcount[,3]
alt_ARA <- pooldata@readcoverage[,4] - pooldata@refallele.readcount[,4]
alt_BMR <- pooldata@readcoverage[,5] - pooldata@refallele.readcount[,5]
alt_CBL <- pooldata@readcoverage[,6] - pooldata@refallele.readcount[,6]
alt_FC <- pooldata@readcoverage[,7] - pooldata@refallele.readcount[,7]
alt_FR <- pooldata@readcoverage[,8] - pooldata@refallele.readcount[,8]
alt_HZD <- pooldata@readcoverage[,9] - pooldata@refallele.readcount[,9]
alt_SBR <- pooldata@readcoverage[,10] - pooldata@refallele.readcount[,10]
alt_PSG <- pooldata@readcoverage[,11] - pooldata@refallele.readcount[,11]
alt_KH <- pooldata@readcoverage[,12] - pooldata@refallele.readcount[,12]
alt_STC <- pooldata@readcoverage[,13] - pooldata@refallele.readcount[,13]
alt_PL <- pooldata@readcoverage[,14] - pooldata@refallele.readcount[,14]
alt_VD <- pooldata@readcoverage[,15] - pooldata@refallele.readcount[,15]
alt_OCT <- pooldata@readcoverage[,16] - pooldata@refallele.readcount[,16]
alt_PB <- pooldata@readcoverage[,17] - pooldata@refallele.readcount[,17]
alt_PGP <- pooldata@readcoverage[,18] - pooldata@refallele.readcount[,18]
alt_PSN <- pooldata@readcoverage[,19] - pooldata@refallele.readcount[,19]

fq_STR <- ref_STR/pooldata@readcoverage[,1]
fq_SLR <- ref_SLR/pooldata@readcoverage[,2]
fq_SH <- ref_SH/pooldata@readcoverage[,3]
fq_ARA <- ref_ARA/pooldata@readcoverage[,4]
fq_BMR <- ref_BMR/pooldata@readcoverage[,5]
fq_CBL <- ref_CBL/pooldata@readcoverage[,6]
fq_FC <- ref_FC/pooldata@readcoverage[,7]
fq_FR <- ref_FR/pooldata@readcoverage[,8]
fq_HZD <- ref_HZD/pooldata@readcoverage[,9]
fq_SBR <- ref_SBR/pooldata@readcoverage[,10]
fq_PSG <- ref_PSG/pooldata@readcoverage[,11]
fq_KH <- ref_KH/pooldata@readcoverage[,12]
fq_STC <- ref_STC/pooldata@readcoverage[,13]
fq_PL <- ref_PL/pooldata@readcoverage[,14]
fq_VD <- ref_VD/pooldata@readcoverage[,15]
fq_OCT <- ref_OCT/pooldata@readcoverage[,16]
fq_PB <- ref_PB/pooldata@readcoverage[,17]
fq_PGP <- ref_PGP/pooldata@readcoverage[,18]
fq_PSN <- ref_PSN/pooldata@readcoverage[,19]


pooldata@snp.info[,1] <- substring(pooldata@snp.info[,1],1,12)
SNP <-paste(pooldata@snp.info[,1],pooldata@snp.info[,2] ,sep="_")


poolfstat_matrix = matrix(nrow=19, ncol=length(SNP))
colnames(poolfstat_matrix) <- SNP
rownames(poolfstat_matrix)=c("STR","SLR","SH","ARA", "BMR", "CBL", "FC", "FR", "HZD", "SBR", "PSG", "KH", "STC", "PL", "VD", "OCT", "PB", "PGP", "PSN")


poolfstat_matrix[1,]=fq_STR
poolfstat_matrix[2,]=fq_SLR
poolfstat_matrix[3,]=fq_SH
poolfstat_matrix[4,]=fq_ARA
poolfstat_matrix[5,]=fq_BMR
poolfstat_matrix[6,]=fq_CBL
poolfstat_matrix[7,]=fq_FC
poolfstat_matrix[8,]=fq_FR
poolfstat_matrix[9,]=fq_HZD
poolfstat_matrix[10,]=fq_SBR
poolfstat_matrix[11,]=fq_PSG
poolfstat_matrix[12,]=fq_KH
poolfstat_matrix[13,]=fq_STC
poolfstat_matrix[14,]=fq_PL
poolfstat_matrix[15,]=fq_VD
poolfstat_matrix[16,]=fq_OCT
poolfstat_matrix[17,]=fq_PB
poolfstat_matrix[18,]=fq_PGP
poolfstat_matrix[19,]=fq_PSN




# Save relative allele frequencies data 
poolfstat_matrix <- as.data.frame(poolfstat_matrix)
pooldata <- read.pcadapt(poolfstat_matrix, type = "pool")


poolfstat_matrix[1:19,]