# Graph Treemix output

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
install.packages(c('ggplot2', 'RColorBrewer'))
library(ggplot2)
library(RColorBrewer)

# Get Treemix R plotting script (in src folder of Treemix source code)
source("src/02_pop_structure/11_treemix_plotting_functions.R")

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
