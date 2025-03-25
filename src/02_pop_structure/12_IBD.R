# Mantel test for isolation by distance (IBD)

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
install.packages(c('ecodist', 'ggplot2', 'RColorBrewer'))
library(ecodist)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Read in geographic distance matrix
geo_dist <- read.csv("data/processed/pop_structure/guide_files/Geo_distance.csv", header=T, row.names=1)

# Read in Fst matrix
fsts <- read.csv("data/processed/pop_structure/Fst/pooldata.pairwisefst.csv", header=T, row.names=1)

# ================================================================================== #

# Make geo dist a matrix
geo_mat <- data.matrix(geo_dist, rownames.force=NA)
# Log transform the geographic distance matrix
geo_mat_log <- log(geo_mat)

# ================================================================================== #

# Run Mantel Test
fst_mantel <- mantel(as.dist(fsts) ~ as.dist(geo_mat_log), nperm=999)
fst_mantel

# ================================================================================== #

# Linearize fsts

# Function to linearize fst
linearize <- function(x) {
  return((x/(1-x)))
} 

# Linearize fsts
fst_lin = apply(fsts, 2, linearize)

# ================================================================================== #

# Run Mantel Test
fst_mantel_lin <- mantel(as.dist(fst_lin) ~ as.dist(geo_mat_log), nperm=999)
fst_mantel_lin

# ================================================================================== #

# Make IBD plot 

e_Dist <- as.matrix(geo_mat_log)
g_Dist <- as.matrix(fst_lin)

df <- data.frame(Genetic_Distance = g_Dist[lower.tri(g_Dist)], Slope_Distance = e_Dist[lower.tri(e_Dist)])
df <- df[ !is.infinite(df$Slope_Distance),]

# Write df
write.csv(df, "output/tables/IBD.plot.csv")
# Add columns classifying pairwise contrast (S, N, X)
df.clust <- read.csv("output/tables/IBD.plot.mod.csv", header=T)

# Specify colors
colors = c("blue", "red", "purple")
colors = c("#0000FF", "#CC0033", "#990099")

pdf("output/figures/pop_structure/IBD/IBD_plot.pdf")
ggplot(df.clust, aes(x=Slope_Distance, y=Genetic_Distance)) + geom_point(aes(fill=Cluster), size=3, shape = 21) + 
scale_fill_manual(values=colors) + stat_smooth(method=lm, formula = y ~ x, colour="black") + 
xlab("log(Geographic Distance)") + ylab("FST/(1-FST)")+ theme_classic(base_size = 16) + guides(fill="none")
dev.off()

