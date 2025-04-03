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

# Make geographic distance a matrix
geo_mat <- data.matrix(geo_dist, rownames.force=NA)
# Log transform the geographic distance matrix
geo_mat_log <- log(geo_mat)

# ================================================================================== #

# Run Mantel Test
fst_mantel <- mantel(as.dist(fsts) ~ as.dist(geo_mat_log), nperm=999)
fst_mantel
#    mantelr       pval1       pval2       pval3   llim.2.5%  ulim.97.5% 
#0.478995128 0.001001001 1.000000000 0.001001001 0.428394132 0.564644807 

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
#    mantelr       pval1       pval2       pval3   llim.2.5%  ulim.97.5% 
#0.480088794 0.001001001 1.000000000 0.001001001 0.436936062 0.558648777 

# ================================================================================== #

# Make Isolation by Distance (IBD) plot 

e_Dist <- as.matrix(geo_mat_log)
g_Dist <- as.matrix(fst_lin)

df <- data.frame(Genetic_Distance = g_Dist[lower.tri(g_Dist)], Slope_Distance = e_Dist[lower.tri(e_Dist)])
df <- df[ !is.infinite(df$Slope_Distance),]

# Write df
write.csv(df, "output/tables/IBD.plot.csv")
# Add columns classifying pairwise contrast (S, N, X)
df.clust <- read.csv("output/tables/IBD.plot.mod.csv", header=T)

# Specify colors
colors = c("#4575b4", "#d73027", "#fee090")

# Graph IBD
pdf("output/figures/pop_structure/IBD/IBD_plot.pdf")
ggplot(df.clust, aes(x=Slope_Distance, y=Genetic_Distance)) + geom_point(aes(fill=Cluster), size=4, shape = 21) + 
scale_fill_manual(values=colors) + stat_smooth(method=lm, formula = y ~ x, colour="black") + 
xlab("log(Geographic Distance)") + ylab("FST/(1-FST)")+ theme_classic(base_size = 25) + guides(fill="none")
dev.off()

