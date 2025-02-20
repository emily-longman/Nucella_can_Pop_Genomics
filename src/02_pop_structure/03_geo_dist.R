# Calculate geographic distance between sites

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
install.packages(c('gdistance', 'ggplot2', 'tidyr', 'dplyr', 'raster', 'sf', 'rnaturalearth', 'nngeo', 'fasterize', 'reshape2', 'patchwork'))
library(gdistance)
library(ggplot2)
library(tidyr)
library(dplyr)
library(raster)
library(sf)
library(rnaturalearth)
library(nngeo)
library(fasterize)
library(reshape2)
library(patchwork)

# ================================================================================== #

# Projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"


