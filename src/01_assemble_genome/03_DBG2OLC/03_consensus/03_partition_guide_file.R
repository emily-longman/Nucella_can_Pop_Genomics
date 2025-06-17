# Create partition file. This file will be used in the next step to break the contigs file and backbone into 50 contig chunks. 

# ================================================================================== #

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo

# Install and load rprojroot
install.packages(c('rprojroot'))
library(rprojroot)
# Create relative path from root
rel_path_from_root <- find_root_file("data", "processed", "genome_assembly", "consensus", criterion = has_file("README.md"))
# Set working directory as path from root
setwd(rel_path_from_root)

# ================================================================================== #

# Load packages
install.packages("foreach")
library(foreach)

# ================================================================================== #

# Create partition file to break the contigs file and backbone into 50 contig chunks.
# Break up the number of contigs (i.e. 27922) into 50 contig chunks

# Create sequence
seq(from=1, to=27922, by= 50) -> vec.num

# Create partitions
dat.win=foreach(i=vec.num, .combine = "rbind")%do%{
  data.frame(start=i, end=i+49)
  }

# Write table
write.table(dat.win, file = "dat.win.partitions.txt", append = FALSE, quote = FALSE, sep = "\t",
eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
