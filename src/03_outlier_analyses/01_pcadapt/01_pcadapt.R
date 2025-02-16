# Use pcadapt to identify outliers
# https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

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
install.packages("pcadapt")
library(pcadapt)
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)

# ================================================================================== #

# Problems either loading the data if in poolfstat pooldata object or can load the full vcf but then can't run pcadapt

# Read in VCF data 
path_to_file <- "data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf"
pooldata <- read.pcadapt(path_to_file, type = "pool")

# Read in metadata
meta_path <- "data/processed/outlier_analyses/guide_files/Populations_metadata.csv"
meta <- read.csv(meta_path, header=T)

# ================================================================================== #

#With Pool-Seq data, the package computes again a Mahalanobis distance based on PCA loadings.

# Note: To choose K, principal component analysis should first be performed with a large enough number of principal components. 
# Use scree plot and PCA to identify optimal K. 

# Run pcadapt on data with large number of K
# You can also set the parameter min.maf that corresponds to a threshold of minor allele frequency. 
# By default, the parameter min.maf is set to 5%
x <- pcadapt(input = pooldata, K = 18)

# Graph screeplot (The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. 
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure 
# lie on a steep curve). Keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule)
pdf("output/figures/pcadapt_screeplot.pdf", width = 10, height = 10)
plot(x, option = "screeplot")
dev.off()

# Graph PCA - PC1 and PC2
pdf("output/figures/pcadapt_pca_PC1_PC2.pdf", width = 10, height = 10)
plot(x, option = "scores", pop = meta$Site)
dev.off()

# Graph PCA - PC3 and PC4
pdf("output/figures/pcadapt_pca_PC3_PC4.pdf", width = 10, height = 10)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
dev.off()

# ================================================================================== #

# Run pcadapt on data with optimal K
x <- pcadapt(input = pooldata, K = 15)

summary(x)

# Identify which PCs SNPs are associate with
get.pc(x, 1:150)->aux
print(aux[,2])

# ================================================================================== #

# Identify outliers 

# Manhattan plot
pdf("output/figures/pcadapt_screeplot.pdf", width = 10, height = 10)
plot(x, option = "manhattan")
dev.off()

# QQ plot
pdf("output/figures/pcadapt_qqplot.pdf", width = 10, height = 10)
plot(x, option = "qqplot")
dev.off()

# Histogram of test statistic and p-val 
pdf("output/figures/pcadapt_hist.pdf", width = 10, height = 10)
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()

# Histogram of the test statistic 𝐷𝑗
pdf("output/figures/pcadapt_stat.dist.pdf", width = 10, height = 10)
plot(x, option = "stat.distribution")
dev.off()

# ================================================================================== #

# Choose cut-off for outlier detection

# q-value method 
# SNPs with q-values less than 𝛼 (10%) will be considered as outliers with an expected false discovery rate bounded by 𝛼
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)

# Benjamini-Hochberg Procedue
padj <- p.adjust(x$pvalues, method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

# Bonferroni Correction
padj <- p.adjust(x$pvalues, method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)


# Get outlier snps associated with PCs
snp_pc <- get.pc(x, outliers)
