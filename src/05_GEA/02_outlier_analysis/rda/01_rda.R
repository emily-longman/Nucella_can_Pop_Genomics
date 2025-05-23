# Redundancy analysis (RDA) as a genotype-environment association (GEA)
# Note: prior to running the R script, need to load R module on the VACC (module load R/4.4.1)
# Useful tutorial: https://popgen.nescent.org/2018-03-27_RDA_GEA.html

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('poolfstat', 'data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'gameofthrones', 'vegan', 'psych'))
library(poolfstat)
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(vegan) # Used to run the RDA
library(psych)

# ================================================================================== #

# Load pooldata object
load("data/processed/pop_structure/pooldata.RData")
str(pooldata)


# Load environmental data

# ================================================================================== #

# Note: RDA requires complete data frames (i.e., no missing data)

# Turn the ref allele read count data object into a dataframe 
readcount <- as.data.frame(pooldata@refallele.readcount)
# Transpose dataframe
readcount_t <- t(readcount)

# Check to make sure no missing data 
sum(is.na(readcount_t))

# Subsample SNP list for testing - 10,000 SNPs
readcount_t_sub <- readcount_t[,sample(1:ncol(readcount_t), 10000)]