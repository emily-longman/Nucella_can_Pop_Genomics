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

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Format read count data and subsample
# Note: RDA requires complete data frames (i.e., no missing data)

# Check structure
str(pooldata)

# Turn the ref allele read count data object into a dataframe 
readcount <- as.data.frame(pooldata@refallele.readcount)
# Transpose dataframe
readcount_t <- t(readcount)

# Set rownames as sites
rownames(readcount_t) <- pooldata@poolnames

# Check to make sure no missing data 
sum(is.na(readcount_t))

# Subsample SNP list for testing - 10,000 SNPs
readcount_t_sub <- readcount_t[,sample(1:ncol(readcount_t), 10000)]

# ================================================================================== #

# Check structure 
str(bio_oracle_sites_2010)

# Re-order bio-oracle data so in same order as read count data
sites <- rownames(readcount_t)
bio_oracle_sites_2010$location = factor(bio_oracle_sites_2010$location, levels = sites)
bio_oracle_sites_2010_ordered <- bio_oracle_sites_2010[order(bio_oracle_sites_2010$location), ]
# Make site names characters, not factors 
bio_oracle_sites_2010_ordered$location <- as.character(bio_oracle_sites_2010_ordered$location)

# Confirm that read count data and environmental data are in the same order
identical(rownames(readcount_t), bio_oracle_sites_2010_ordered[,13]) 

# ================================================================================== #

# Assess correlations among the Bio-oracle environmental variables
# Note: typically remove predictors with |r| > 0.7

# Bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal
pdf("output/figures/GEA/enviro/Bio-oracle/Bio-oracle_correlations.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_ordered[,4:12], scale=T)
dev.off()

# Many of the variables are correlated
# Remove temp max, temp min, temp range, o2, 
bio_oracle_sites_2010_ordered_sub <- bio_oracle_sites_2010_ordered[,-c(4,5,6,11,13)]

pdf("output/figures/GEA/enviro/Bio-oracle/Bio-oracle_correlations_sub.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_ordered_sub[,4:7], scale=T)
dev.off()

# Everything is super correlated! I don't know which to keep
# One option is to keep: thetao_mean, chl_mean, O2_mean, and ph_min

