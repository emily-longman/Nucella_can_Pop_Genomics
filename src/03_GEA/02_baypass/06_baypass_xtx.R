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
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages(c('data.table', 'dplyr', 'ggplot2', 'mvtnorm', 'geigen', 'tidyverse', 'foreach'))
library(data.table)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(geigen)
library(tidyverse)
library(foreach)
library(WriteXLS)

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
abline(h=3, lty=2, col="red") #0.001 p-value threshold
dev.off()


# ================================================================================== #

# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Rank normalize p-values

L = dim(SNP.XtX.dt)[1] 
    
SNP.XtX.dt %>%
arrange(-log10.1.pval.) %>% 
as.data.frame() %>%
mutate(rank = seq(from=1,  to = L, by = 1)) %>% 
mutate(rank_norm = rank/L) %>% 
group_by(chr) %>%
arrange(as.numeric(pos))->
inner.rnf


# Window and step size
win.bp = 100000
step.bp = 50000

# Generate windows
wins <- foreach(chr.i=unique(SNP.XtX.dt$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  tmp <-  inner.rnf %>%
                    filter(chr == chr.i)

                    if(max(tmp$pos) <= step.bp){
                      data.table(chr=chr.i,
                             start=min(tmp$pos),
                             end=max(tmp$pos))  
                    } else if(max(tmp$pos) > step.bp){
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                } }

wins[,i:=1:dim(wins)[1]]

# Save windows
save(wins, file="data/processed/GEA/baypass/baypass_windows.RData")
save.image("data/processed/GEA/baypass/N.canaliculata_baypass_windows.RData")

# Reload windows
load("data/processed/GEA/baypass/baypass_windows.RData")

# ================================================================================== #

# Start the summarization process

win.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- inner.rnf %>%
    filter(chr == wins[win.i]$chr) %>%
    filter(pos >= wins[win.i]$start & pos <= wins[win.i]$end)
  
  pr.i <- c(0.01)
  
  win.tmp %>% 
    filter(!is.na(rank_norm)) %>%
    summarise(chr = wins[win.i]$chr,
              pos_mean = mean(pos),
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              win=win.i,
              pr=pr.i,
              rnp.pr=c(mean(rank_norm<=pr.i)),
              rnp.binom.p=c(binom.test(sum(rank_norm<=pr.i), 
                                       length(rank_norm), pr.i)$p.value),
              max.p=max(log10.1.pval.),
              nSNPs = n(),
              sum.rnp=sum(rank_norm<=pr.i),
    )  -> win.out
}

# Save window out
save(win.out, file="data/processed/GEA/baypass/baypass_win.out.RData")

# Reload windows
load("data/processed/GEA/baypass/baypass_win.out.RData")

# ================================================================================== #

# Graph 

# Create unique Chromosome number 
chr.unique <- unique(win.out$chr)
win.out$chr.unique <- as.numeric(factor(win.out$chr, levels = chr.unique))

# Graph based on chr
pdf("output/figures/GEA/Baypass_rnp.pdf", width = 5, height = 5)
ggplot(win.out, aes(y=rnp.binom.p, x=chr.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()
dev.off()

# Graph based on position
pdf("output/figures/GEA/Baypass_rnp_pos.pdf", width = 5, height = 5)
ggplot(win.out, aes(y=rnp.binom.p, x=pos_mean/1e6)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()
dev.off()

##### NOTE: THESE GRAPHS LOOK OFF - should I not use -log10(0.01) or are all of the 1s because I am hitting limits becuase contigs are short?

# ================================================================================== #

# Identify contigs with highly significant rnp p
win.out.sig <- win.out[which(-log10(win.out$rnp.binom.p)>-log10(0.01)),]

# ================================================================================== #

# Create outlier SNP list 
SNPs.Interest <- foreach(i=1:dim(win.out.sig)[1], .combine = "rbind")%do%{
  tmp.snps <- data.binary.SNP.filt %>%
  filter(Chromosome == win.out.sig[i,]$Chromosome) %>%
  filter(Position >= win.out.sig[i,]$pos_min & Position <= win.out.sig[i,]$pos_max)
}

# Write file of outlier SNPs
write.csv(SNPs.Interest, "Nucella_GWAS_outlier_SNPs.csv")

