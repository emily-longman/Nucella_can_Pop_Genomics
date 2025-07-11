# Analyze the SNPs within the protein of interest
# Note: prior to running the R script, need to load R and GDAL module on the VACC
# module load R/4.4.1
# module load gdal

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
install.packages(c('RColorBrewer', 'data.table', 'ggplot2', 'poolfstat', 'tidyverse', 'foreach', 'magrittr', 'rnaturalearth', 'rnaturalearthdata'))
library(RColorBrewer)
library(data.table)
library(ggplot2)
require(poolfstat)
library(tidyverse)
library(foreach)
library(magrittr)
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# ================================================================================== #

# Read in datasets

# Metadata
meta <- fread("guide_files/Populations_metadata.csv")

# Morphology data
morph <- fread("data/processed/morphometrics/pc.morphology.csv")
morph_all_data <- fread("data/processed/morphometrics/morphology_geomorph_data.csv")
# Rename sites column for morph
names(morph)[2] = "Site"

# Load subset of SNPS
load("data/processed/pop_structure/pooldata_subset_prot_g26813.RData")

# ================================================================================== #

# Extract and manipulate snp info

# Extract SNP info for all SNPs
pooldata_prot_g26813@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  snp.info.annot

# Rename columns
names(snp.info.annot)[1:2] = c("chr","pos")

# Make snp_id column
snp.info.annot <- snp.info.annot %>%
  mutate(snp_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Extract and manipulate coverage

# Extract read count and coverage data for SNPs
ref_count <- pooldata_prot_g26813@refallele.readcount
coverage <- pooldata_prot_g26813@readcoverage

# Extract coverage for SNPs of interest
coverage %>%
  as.data.frame %>%
  mutate(snp_id = snp.info.annot$snp_id ) ->
  covs.id

# Rename columns (19 sites and snp_id)
names(covs.id) = c(pooldata_prot_g26813@poolnames, "snp_id")

# Restructure data so long format
reshape2::melt(covs.id, id = "snp_id", variable.name = "Site", value.name = "COV") -> covs.id.melt

# ================================================================================== #

# Calculate allele frequencies

# Calculate allele frequency for SNPs
afs <- ref_count/coverage

# Extract allele freq for SNPs of interest
afs %>%
  as.data.frame %>%
  mutate(snp_id = snp.info.annot$snp_id ) ->
  afs.id

# Rename columns (19 sites and snp_id)
names(afs.id) = c(pooldata_prot_g26813@poolnames, "snp_id")

# Restructure data so long format
reshape2::melt(afs.id, id = "snp_id", variable.name = "Site", value.name = "AF") %>%
  separate(snp_id, remove = F, into = c("chr_Or", "chr_id", "pos"), sep = "_") -> afs.id.melt

# ================================================================================== #

# Join datasets

# Join datasets
left_join(covs.id.melt, afs.id.melt, by = join_by(snp_id, Site)) %>%
  left_join(meta, by = join_by(Site)) %>%
  left_join(morph, by = join_by(Site)) ->
  afs.id.mapped

setDT(afs.id.mapped)

# Calculate mean effective coverage ('nEff')
afs.id.mapped %>% mutate(nEff:=round((COV*2*40 )/(COV+2*40- 1))) %>%
  mutate(af_nEff:=round(AF*nEff)/nEff) ->
  afs.id.mapped

# Add demography information - North (N), Admixed (ADMX), South (S)
afs.id.mapped %<>%
  mutate(cluster = case_when(
    Site %in% c("FC","SLR","SH","ARA","CBL","PSG","STC","KH","VD","FR","BMR","PGP") ~ "N",
    Site %in% c("PL","SBR") ~ "ADMX",
    Site %in% c("PSN","PB","HZD","OCT","STR") ~ "S"
  )) 

# ================================================================================== #

# Filter dataset for only candidate SNPs
afs.id.mapped %>%
  filter(snp_id %in% c(
    "ntLink_3633_41113",
    "ntLink_3633_41601",
    "ntLink_3633_41553",
    "ntLink_3633_42901"
  )) -> afs.id.mapped.targets

# Make site a factor 
afs.id.mapped.targets$Site <- factor(afs.id.mapped.targets$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Graph allele frequencies

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# AF versus latitude
pdf("output/figures/morphology/Allele_freq_lat.pdf", width = 18, height = 4)
ggplot(data = afs.id.mapped.targets, aes(x = Lat, y = AF)) +
  geom_point(aes(fill = Site), size = 5, shape = 21) + 
  geom_vline(xintercept=36.8007, linetype="solid", color="black") +
  ylab("Allele Frequency") + xlab("Latitude") +
  scale_fill_manual(values = mycolors) +
  facet_grid(~pos) + theme_bw(base_size = 22) + theme(legend.position="none")
dev.off()


# Graph map
pdf("output/figures/morphology/Allele_freq_maps.pdf", width = 8, height = 3)
ggplot(data = world) +
  geom_sf(fill= "grey70") +
  coord_sf(xlim = c(-127, -115), ylim = c(32, 47), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "white")) +
  geom_jitter(data = afs.id.mapped.targets,
              color = "black",
              aes(
                x=Long,
                y=Lat,
                fill = AF,
                shape = cluster
              ), alpha = 0.9, size = 2.5) + 
  scale_fill_gradient2(low = "steelblue", ,high = "firebrick",
                       midpoint = 0.5) +
  scale_shape_manual(values = 21:23) +
  facet_grid(~snp_id) + theme_bw()
dev.off()

# Graph allele freq on morphology PC 
pdf("output/figures/morphology/Allele_freq_pcplots.pdf", width = 18, height = 4)
afs.id.mapped.targets %>%
  ggplot(aes(y=mean.pc2, x=mean.pc1)) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_point(aes(fill = AF), shape = 21, size = 5) + 
  ylab("Mean PC2") + xlab("Mean PC1") +
  scale_fill_gradient2(low = "#23d458", ,high = "black", midpoint = 0.5)  +
  facet_grid(~pos) + theme_bw(base_size = 25)
dev.off()

# Graph allele freq on morphology PC alternative
pdf("output/figures/morphology/Allele_freq_pcplots_2.pdf", width = 8, height = 3)
afs.id.mapped.targets %>%
  ggplot(aes(y=mean.pc2, x=mean.pc1, fill = AF, shape = cluster)) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  geom_point(size = 3.1) + 
  scale_shape_manual(values = 21:23) +
  scale_fill_gradient2(low = "steelblue", ,high = "firebrick",
                       midpoint = 0.5)  +
  facet_grid(~snp_id) + theme_bw() 
dev.off()


# ================================================================================== #
# ================================================================================== #

# Graph Baypass outliers and highlight top SNPs

# ================================================================================== #

# Read in Baypass results 
baypass_morph_CV1 <- read.table("data/processed/morphometrics/baypass/CV_1/NC_baypass_morphology_CV_1_summary_pi_xtx.out", header=T)
baypass_morph_CV2 <- read.table("data/processed/morphometrics/baypass/CV_2/NC_baypass_morphology_CV_2_summary_pi_xtx.out", header=T)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Join baypass morphology CV results with snp metadata
baypass_morph_CV1_pos <- cbind(snp.meta, baypass_morph_CV1)
baypass_morph_CV2_pos <- cbind(snp.meta, baypass_morph_CV2)

# Make snp_id column
baypass_morph_CV1_pos <- baypass_morph_CV1_pos %>%
  mutate(snp_id = paste(chr, pos, sep = "_"))
baypass_morph_CV2_pos <- baypass_morph_CV2_pos %>%
  mutate(snp_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Highlight outliers in graphs

# Set colors for outliers
baypass_morph_CV1_pos$color <- ifelse(baypass_morph_CV1_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV1_pos)[1]), "sig", "not.sig")
baypass_morph_CV2_pos$color <- ifelse(baypass_morph_CV2_pos$log10.1.pval. >= -log10(0.05/dim(baypass_morph_CV2_pos)[1]), "sig", "not.sig")

# Color top SNPS
baypass_morph_CV1_pos$color[which(baypass_morph_CV1_pos$snp_id == "ntLink_3633_41113" | baypass_morph_CV1_pos$snp_id == "ntLink_3633_41601" )] <- "top"
baypass_morph_CV2_pos$color[which(baypass_morph_CV2_pos$snp_id == "ntLink_3633_41553" | baypass_morph_CV2_pos$snp_id == "ntLink_3633_41601" | 
baypass_morph_CV2_pos$snp_id == "ntLink_3633_42901")] <- "top"

# Change to factor
baypass_morph_CV1_pos$color <- factor(baypass_morph_CV1_pos$color)
baypass_morph_CV2_pos$color <- factor(baypass_morph_CV2_pos$color)

# CV1
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_1_alt_top_SNPs.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV1_pos$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position", 
col=c("#bebebe93", "black", "red")[baypass_morph_CV1_pos$color], pch=19)
abline(h=-log10(0.001), lty=2, col="#1515ed") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV1_pos)[1]), lty=1, col="#1515ed") # Bonferroni threshold
dev.off()

# CV2
pdf("output/figures/morphology/baypass/Baypass_xtx_outliers_CV_2_alt_top_SNPs.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(baypass_morph_CV2_pos$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position", 
col=c("#bebebe93", "black", "red")[baypass_morph_CV2_pos$color], pch=19)
abline(h=-log10(0.001), lty=2, col="#1515ed") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass_morph_CV2_pos)[1]), lty=1, col="#1515ed") # Bonferroni threshold
dev.off()
