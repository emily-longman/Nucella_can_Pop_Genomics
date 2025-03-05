# Analyze grenedalf output pop gen statistics

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

# ================================================================================== #

# Load data
meta <- read.csv("data/processed/pop_gen/guide_files/Populations_metadata.csv", header=T)
grenedalf <- read.csv("data/processed/pop_gen/grenedalf_filt/diversity.csv", header=F)

# ================================================================================== #

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
colors.alphabetical <- mycolors[c(4,11,5,1,10,17,8,18,16,12,13,6,15,14,3,2,7,19,9)]

# ================================================================================== #

# Extract and reformat data for theta pi
theta_pi <- t(grenedalf[, grep("theta_pi", grenedalf[1,])])
theta_pi_names <- data.frame(do.call('rbind', strsplit(as.character(theta_pi[,1]),'.',fixed=TRUE)))
theta_pi <- data.frame(theta_pi_names[,1], as.numeric(theta_pi[,2]))
colnames(theta_pi) <- c("Site", "theta_pi")
theta_pi <- merge(theta_pi, meta, by="Site")

pdf("output/figures/pop_structure/Theta_pi_filt.pdf", width = 5, height = 5)
ggplot(data = theta_pi, aes(x = Lat, y = theta_pi)) + 
  geom_point(shape = 21, size = 3, fill = colors.alphabetical) + 
  xlab("Latitude") + ylab("Theta pi") + theme_classic() 
dev.off()

# ================================================================================== #

# Extract and reformat data for theta watterson
theta_watterson <- t(grenedalf[, grep("theta_watterson", grenedalf[1,])])
theta_watterson_names <- data.frame(do.call('rbind', strsplit(as.character(theta_watterson[,1]),'.',fixed=TRUE)))
theta_watterson <- data.frame(theta_watterson_names[,1], as.numeric(theta_watterson[,2]))
colnames(theta_watterson) <- c("Site", "theta_watterson")
theta_watterson <- merge(theta_watterson, meta, by="Site")

pdf("output/figures/pop_structure/Theta_watterson_filt.pdf", width = 5, height = 5)
ggplot(data = theta_watterson, aes(x = Lat, y = theta_watterson)) + 
  geom_point(shape = 21, size = 3, fill = colors.alphabetical) + 
  xlab("Latitude") + ylab("Theta Watterson") + theme_classic() 
dev.off()


# ================================================================================== #

# Extract and reformat data for Tajima D
tajimas_d <- t(grenedalf[, grep("tajimas_d", grenedalf[1,])])
tajimas_d_names <- data.frame(do.call('rbind', strsplit(as.character(tajimas_d[,1]),'.',fixed=TRUE)))
tajimas_d <- data.frame(tajimas_d_names[,1], as.numeric(tajimas_d[,2]))
colnames(tajimas_d) <- c("Site", "tajimas_d")
tajimas_d <- merge(tajimas_d, meta, by="Site")

pdf("output/figures/pop_structure/Tajimas_d_filt.pdf", width = 5, height = 5)
ggplot(data = tajimas_d, aes(x = Lat, y = tajimas_d)) + 
  geom_point(shape = 21, size = 3, fill = colors.alphabetical) + 
  xlab("Latitude") + ylab("Tajima's d") + theme_classic() 
dev.off()

# A negative Tajima's D signifies an excess of low frequency polymorphisms relative to expectation, indicating population size expansion
