# Compare gtf for cDNA and ab initio methods

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move the covmatrix to the pcangsd folder in results > stats

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Set relative path of results directory from root
dir(find_root_file("results", criterion = has_file("README.md")))
rel_path_from_root <- find_root_file("results", "genome_annotation", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(rel_path_from_root)

# Set working directory as path from root
setwd(rel_path_from_root)

# ================================================================================== #

# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Load data
ab_initio <- read.table("braker_ab_initio.gtf", header = F, sep = '\t')
cDNA <- read.table("braker_cDNA.gtf", header = F, sep = '\t')

str(ab_initio)
str(cDNA)

ab_initio$V6 <- as.numeric(ab_initio$V6) # NAs will be generated for some elements
cDNA$V6 <- as.numeric(cDNA$V6) # NAs will be generated for some elements

# ================================================================================== #

# Calculate length of element
ab_initio$length <- ab_initio$V5 - ab_initio$V4
cDNA$length <- cDNA$V5 - cDNA$V4

# Calculate mid-point of element
ab_initio$mid.point <- (ab_initio$length/2) + ab_initio$V4
cDNA$mid.point <- (cDNA$length/2) + cDNA$V4

# Calculate number of each element
ab_initio_sum <- ab_initio %>%
  group_by(Element = V3) %>% 
  summarise(n=n(), mean.length = mean(length), tot.length = sum(length))

cDNA_sum <- cDNA %>%
  group_by(Element = V3) %>% 
  summarise(n=n(), mean.length = mean(length), tot.length = sum(length))

# ================================================================================== #

# Density graphs of introns and exons

# Subset data for only introns
ab_initio_introns <- ab_initio[which(ab_initio$V3 == "intron"),]
cDNA_introns <- cDNA[which(cDNA$V3 == "intron"),]

# Density graph of intron length
ab_initio_introns_graph <- ggplot(ab_initio_introns, aes(x=length)) +
  geom_density(col="black", alpha=0.8) + 
  xlab("Intron Length (bp)") + xlim(0,1500) + 
  ylab("Density") +
  theme_bw()
cDNA_introns_graph <- ggplot(cDNA_introns, aes(x=length)) + 
  geom_density(col="black", alpha=0.8, linetype = "dashed") + 
  xlab("Intron Length (bp)") + xlim(0,1500) + 
  ylab("Density") +
  theme_bw() 


####

# Subset data for only exons
ab_initio_exons <- ab_initio[which(ab_initio$V3 == "exon"),]
cDNA_exons <- cDNA[which(cDNA$V3 == "exon"),]

# Density graph of exon length
ab_initio_exons_graph <- ggplot(ab_initio_exons, aes(x=length)) +
  geom_density(col="black", alpha=0.8) + 
  xlab("Exon Length (bp)") + xlim(0,500) + 
  ylab("Density") +
  theme_bw()
cDNA_exons_graph <- ggplot(cDNA_exons, aes(x=length)) + 
  geom_density(col="black", alpha=0.8, linetype = "dashed") + 
  xlab("Exon Length (bp)") +  xlim(0,500) +
  ylab("Density") +
  theme_bw() 

plot_grid(ab_initio_introns_graph, cDNA_introns_graph, ab_initio_exons_graph, cDNA_exons_graph, ncol=2)

# ================================================================================== #

# Compare posterior for both methods

# Subset data for only CDS
ab_initio_CDS <- ab_initio[which(ab_initio$V3 == "CDS"),]
cDNA_CDS <- cDNA[which(cDNA$V3 == "CDS"),]

# Drop NA rows
ab_initio_CDS_subset <- ab_initio_CDS %>% drop_na()
cDNA_CDS_subset <- cDNA_CDS %>% drop_na()

# Compare 

ggplot(ab_initio_CDS_subset, aes(x=mid.point, y=V6)) + 
  geom_point(col="black",alpha=0.8, size=1.3) + 
  ylab("Posterior") + xlab("Position (bp)") +
  theme_bw() 
  
ggplot(cDNA_subset, aes(x=mid.point, y=V6)) + 
  geom_point(col="black",alpha=0.8, size=1.3) + 
  ylab("Posterior") + xlab("Position (bp)") +
  theme_bw() 

# ================================================================================== #

# Drop NA rows
ab_initio_subset <- ab_initio %>% drop_na()
cDNA_subset <- cDNA %>% drop_na()


ab_initio_posterior <- ggplot(ab_initio_subset, aes(x = V3, y = V6)) + 
  geom_boxplot() +
  ylab("Posterior") + xlab("Element") +
  theme_bw()

cDNA_posterior <- ggplot(cDNA_subset, aes(x = V3, y = V6)) + 
  geom_boxplot() +
  ylab("Posterior") + xlab("Element") +
  theme_bw()

plot_grid(ab_initio_posterior, cDNA_posterior, ncol=2)
