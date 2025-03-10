# Graph Baypass XtX output.

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
install.packages(c('data.table', 'dplyr', 'ggplot2'))
library(data.table)
library(dplyr)
library(ggplot2)

# ================================================================================== #

# Read in Baypass results
NC.omega <- as.matrix(read.table("data/processed/GEA/baypass/omega/NC_baypass_mat_omega.out"))
XtX <- read.table("data/processed/GEA/baypass/xtx/NC_baypass_core_summary_pi_xtx.out", header=T)
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
snp.meta <- read.table("data/processed/GEA/baypass/snpdet", header=F)

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
 
# ================================================================================== #

# Check the behavior of the p-values associated to the XtXst estimator
pdf("output/figures/GEA/Baypass_xtx_hist.pdf", width = 5, height = 5)
hist(10**(-1*XtX$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/GEA/Baypass_xtx.pdf", width = 5, height = 5)
plot(XtX$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/GEA/Baypass_xtx_outliers.pdf", width = 5, height = 5)
plot(XtX$log10.1.pval., ylab="XtX P-value (-log10 scale)" )
abline(h=3, lty=2) #0.001 p-value threshold
dev.off()

# ================================================================================== #

# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Import local score functions (https://forge-dga.jouy.inra.fr/documents/809).
source('data/processed/GEA/baypass/xtx/localscore/scorelocalfunctions.R')

# Create key based on the chromosome
setkey(SNP.XtX.dt, chr)
# Total number of chromosomes (18,369)
Nchr=length(SNP.XtX.dt[,unique(chr)])

# Compute the absolute position in the genome. This is useful for doing genome-wide plots.
chrInfo=SNP.XtX.dt[,.(L=.N, cor=autocor(log10.1.pval.)), chr]
setkey(chrInfo, chr)
tmp=data.table(chr=SNP.XtX.dt[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr])))
setkey(tmp,chr)
SNP.XtX.dt[tmp,posT:=pos+S]


# Choice of $\xi$ (1,2,3 or 4)

# To choose the apropiate threshold ($\xi = 1$ or 2) we look at the distribution of −log10(p − value). Then the score function will be $−log10(p − value) − \xi$. 

pdf("output/figures/GEA/Baypass_xtx_lindley_qplot.pdf", width = 5, height = 5)
qplot(log10.1.pval., data=SNP.XtX.dt, geom='histogram', binwidth=0.1, main='P-values histogram')
dev.off()

mean(-log10(SNP.XtX.dt$log10.1.pval.)) #0.6148583

# Remember that we have to choose some value between mean(-log10(mydata$pval)) and max(-log10(mydata$pval)).
mean(SNP.XtX.dt$log10.1.pval.)  #0.4177791
max(SNP.XtX.dt$log10.1.pval.) #10.44488

pdf("output/figures/GEA/Baypass_xtx_lindley_plot.pdf", width = 5, height = 5)
p<-ggplot(data=SNP.XtX.dt, aes(posT,log10.1.pval.)) + geom_point()
p + geom_abline(intercept=c(1,2) , slope=0) 
dev.off()

# Computation of the score and the Lindley Process

xi=2
SNP.XtX.dt[,score:= log10.1.pval.-xi]
# The score mean must be negative
mean(SNP.XtX.dt$score) # -1.582221
SNP.XtX.dt[,lindley:=lindley(score),chr]

# Compute significance threshold for each chromosome

# Uniform distribution of p-values

# If the distribution of the p-values is uniform and if $\xi=1$ or 2, 
# it is possible to compute the thresholds of significancy for each chromosome directely, given the length and the autocorrelation of the chromosome .

chrInfo[,th:=thresUnif(L, cor, 2),chr]
SNP.XtX.dt=SNP.XtX.dt[chrInfo]
sigZones=SNP.XtX.dt[chrInfo, sig_sl(lindley, pos, unique(th)),chr]



