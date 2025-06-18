# Analyze npstat output pop gen statistics

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
install.packages(c('ggplot2', 'RColorBrewer', 'dplyr', 'WriteXLS', 'lme4', 'emmeans', 'lmerTest', 'car'))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(WriteXLS)
library(lme4)
library(emmeans)
library(lmerTest)
library(car)

# ================================================================================== #

# Generate Folders and files

# Make output directory
output_dir_npstat="output/figures/pop_structure/npstat"
if (!dir.exists(output_dir_npstat)) {dir.create(output_dir_npstat)}

# ================================================================================== #

# Load data
meta <- read.csv("guide_files/Populations_metadata_alphabetical.csv", header=T)
pops <- meta$Site

# Load data for each population
for (i in pops){
    print(i)
    filename <- paste0("npstat.", i)
    wd <- paste0("data/processed/pop_structure/npstat/npstat.", i, ".txt")
    assign(filename, read.table(wd, header=F))
}

# Rename columns for each population
column.names <- c("window", "length", "length_outgroup", "read_depth", "S", "Watterson", "Pi", "Tajima_D", "var_S", 
"var_Watterson", "unnorm_FayWu_H", "FayWu_H", "div", "nonsyn_pol", "syn_pol", "nonsyn_div", "syn_div", "alpha", "chr_name")  

colnames(npstat.ARA) <- column.names
colnames(npstat.BMR) <- column.names
colnames(npstat.CBL) <- column.names
colnames(npstat.FC) <- column.names
colnames(npstat.FR) <- column.names
colnames(npstat.HZD) <- column.names
colnames(npstat.KH) <- column.names
colnames(npstat.OCT) <- column.names
colnames(npstat.PB) <- column.names
colnames(npstat.PGP) <- column.names
colnames(npstat.PL) <- column.names
colnames(npstat.PSG) <- column.names
colnames(npstat.PSN) <- column.names
colnames(npstat.SBR) <- column.names
colnames(npstat.SH) <- column.names
colnames(npstat.SLR) <- column.names
colnames(npstat.STC) <- column.names
colnames(npstat.STR) <- column.names
colnames(npstat.VD) <- column.names

# Note: becuase these statistics are calculated based on the bam files, each file contains a slightly different number of rows.

# ================================================================================== #

# Format dataframe

# Merge all datasets into one long dataset
npstat <- bind_rows(npstat.ARA, npstat.BMR, npstat.CBL, npstat.FC, npstat.FR, npstat.HZD, npstat.KH, npstat.OCT, npstat.PB, 
npstat.PGP, npstat.PL, npstat.PSG, npstat.PSN, npstat.SBR, npstat.SH, npstat.SLR, npstat.STC, npstat.STR, npstat.VD)

# Split chr_name column into multiple columns (the first column contains the scaffold name)
npstat.chr.names <- data.frame(do.call('rbind', strsplit(as.character(npstat$chr_name),'.',fixed=TRUE)))
# Remove last 2 columns (they just say "pileup" "stats")
npstat.chr.names <- npstat.chr.names[,1:2]
# Rename columns
colnames(npstat.chr.names) <- c("chr", "pop")
# Join data frames
npstat <- cbind(npstat, npstat.chr.names)

# ================================================================================== #

# Graph Statistics

# Make population an ordered factor for graphing
npstat$pop <- factor(npstat$pop, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph statistics and express as log(pi), log(wattersons), and tajima's D 

# Graph Pi
pdf("output/figures/pop_structure/npstat/Npstat_boxplot_logPi.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(x = pop, y = log(Pi), fill=pop)) + 
  geom_boxplot() + 
  scale_color_manual(values=mycolors) +
  xlab("Population") + ylab("log(Pi)") + theme_classic(base_size=20) +
  guides(fill="none")
dev.off()

# Graph Pi - density
pdf("output/figures/pop_structure/npstat/Npstat_density_Pi.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(log(Pi), color=pop)) + 
  geom_density(size=2) + 
  scale_color_manual(values=mycolors) + ylim(NA,0.955) +
  xlab("log(Pi)") + ylab("Density") + theme_classic(base_size=20) +
  theme(legend.position="none")
dev.off()

# Graph Watterson Theta
pdf("output/figures/pop_structure/npstat/Npstat_boxplot_Watterson.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(x = pop, y = log(Watterson), fill=pop)) + 
  geom_boxplot() + 
  scale_color_manual(values=mycolors) +
  xlab("Population") + ylab("log(Watterson)") + theme_classic(base_size=20) +
  guides(fill="none")
dev.off()

# Graph Watterson Theta - density
pdf("output/figures/pop_structure/npstat/Npstat_density_Watterson.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(log(Watterson), color=pop)) + 
  geom_density(size=2) + 
  scale_color_manual(values=mycolors) + ylim(NA,0.955) +
  xlab("log(Watterson)") + ylab("Density") + theme_classic(base_size=20) +
  theme(legend.position="none")
dev.off()

# Graph Tajima's D - boxplot
pdf("output/figures/pop_structure/npstat/Npstat_boxplot_TajimaD.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(x = pop, y = Tajima_D, fill=pop)) + 
  geom_boxplot() + 
  scale_color_manual(values=mycolors) +
  xlab("Population") + ylab("Tajima's D") + theme_classic(base_size=20) +
  guides(fill="none")
dev.off()

# Graph Tajima's D - density
pdf("output/figures/pop_structure/npstat/Npstat_density_TajimaD.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(Tajima_D, color=pop)) + 
  geom_density(size=2) + 
  scale_color_manual(values=mycolors) +
  xlab("Tajima's D") + ylab("Density") + theme_classic(base_size=20) +
  geom_vline(xintercept=0, linetype="dashed") + theme(legend.position="none")
dev.off()

# Note: data table has NAs and inf values - these will generate a warning due to taking the log

# ================================================================================== #

# Generate summary table

# Create summary table for each population
npstat.summary <- npstat %>% 
  dplyr::group_by(pop) %>% 
  dplyr::summarise(pi_mean = mean(Pi, na.rm = TRUE), pi_median = median(Pi, na.rm = TRUE), 
  Tajima_D_mean = mean(Tajima_D, na.rm = TRUE), Tajima_D_median = median(Tajima_D, na.rm = TRUE))

# Write table
write.table(npstat.summary, "output/tables/npstat.summary.txt", sep='\t')

