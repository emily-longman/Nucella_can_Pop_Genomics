# Analyze omega matrix

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
install.packages(c('corrplot', 'ggplot2', 'RColorBrewer', 'WriteXLS', 'ape'))
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(WriteXLS)
library(ape)

# ================================================================================== #

# Load data
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
NC.omega <- as.matrix(read.table("data/processed/GEA/baypass/omega/NC_baypass_mat_omega.out"))

# ================================================================================== #

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

write.table(NC.omega, "output/tables/NC.omega.txt", sep='\t')

# Re-order populations so in latitudinal order
NC.omega_reorder <- as.matrix(read.table("output/tables/NC.omega_reorder.txt", header=T))

# ================================================================================== #

# Create a correlation matrix of the omega values -- assess genomic differentiation between pools
cor.mat <- cov2cor(NC.omega_reorder)

pdf("output/figures/GEA/Baypass_omega_cor_matrix.pdf", width = 5, height = 5)
corrplot(cor.mat, method = "color", mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))
dev.off()

# ================================================================================== #

# Assess population differentiation with hierarchical clustering
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))

pdf("output/figures/GEA/Baypass_hier_clustering.pdf", width = 5, height = 5)
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

# ================================================================================== #

# Read the xtx BayPass output
NC.snp.res <- read.table("data/processed/GEA/baypass/omega/NC_baypass_summary_pi_xtx.out")

# Get the Pi Beta distribution for POD generation
NC.pi.beta.coef <- read.table("data/processed/GEA/baypass/omega/NC_baypass_summary_beta_params.out",h=T)$Mean

# Upload original data to get read counts
NC.data <- geno2YN("SG.genobaypass")

# Simulate POD dataset to use for outlier SNP detection
simu.SG <-simulate.baypass(omega.mat=SG.omega,nsnp=10000,sample.size=SG.data$NN,
                           beta.pi=SG.pi.beta.coef,pi.maf=0,suffix="SG.BP.sim")