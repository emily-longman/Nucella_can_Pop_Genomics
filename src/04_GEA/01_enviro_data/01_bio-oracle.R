# Get Bio-oracle data

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Create new library for package

# Check library paths
.libPaths() 
# Create new path for Biooracle libraries
.libPaths(c("/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Rpackages_biooracle", .libPaths()))
# Set library path to biooracle path
libpath <- .libPaths()[1]

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
raw_data_path_from_root <- find_root_file("data", "raw", "Bio-oracle", criterion = has_file("README.md"))
# Set working directory as path from root
setwd(raw_data_path_from_root)
setwd("/Users/emilylongman/Documents/GitHub/Nucella_can_Pop_Genomics/data/raw/Bio-oracle")

# ================================================================================== #

# Load packages
install.packages('devtools')
library(devtools)


install.packages("Rcpp")
#library(Rcpp)
install.packages('terra')
library(terra)


devtools::install_github("bio-oracle/biooracler")
library(biooracler)
library(tidyverse)

# ================================================================================== #

# 