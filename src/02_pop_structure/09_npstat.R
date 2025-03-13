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

# Load data
meta <- read.csv("data/processed/pop_gen/guide_files/Populations_metadata_alphabetical.csv", header=T)
pops <- meta$Site

# Load data for each population
for (i in pops){
    print(i)
    filename <- paste0("npstat.", i)
    wd <- paste0("data/processed/pop_gen/npstat/npstat.", i, ".txt")
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

# Graph Pi
pdf("output/figures/pop_structure/diversity_stats/npstat/Pi_npstat_boxplot_logPi.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(x = pop, y = log(Pi), fill=pop)) + 
  geom_boxplot() + 
  scale_color_manual(values=mycolors) +
  xlab("Population") + ylab("log(Pi)") + theme_classic(base_size=20) +
  guides(fill="none")
dev.off()

# Graph Pi - density
pdf("output/figures/pop_structure/diversity_stats/npstat/Pi_npstat_density_Pi.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(log(Pi), color=pop)) + 
  geom_density() + 
  scale_color_manual(values=mycolors) +
  xlab("log(Pi)") + ylab("Density") + theme_classic(base_size=20) +
  theme(legend.position="none")
dev.off()

# Graph Watterson Theta
pdf("output/figures/pop_structure/diversity_stats/npstat/Pi_npstat_boxplot_Watterson.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(x = pop, y = log(Watterson), fill=pop)) + 
  geom_boxplot() + 
  scale_color_manual(values=mycolors) +
  xlab("Population") + ylab("log(Watterson)") + theme_classic(base_size=20) +
  guides(fill="none")
dev.off()

# Graph Watterson Theta - density
pdf("output/figures/pop_structure/diversity_stats/npstat/Pi_npstat_density_Watterson.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(log(Watterson), color=pop)) + 
  geom_density() + 
  scale_color_manual(values=mycolors) +
  xlab("log(Watterson)") + ylab("Density") + theme_classic(base_size=20) +
  theme(legend.position="none")
dev.off()

# Graph Tajima's D - boxplot
pdf("output/figures/pop_structure/diversity_stats/npstat/Pi_npstat_boxplot_TajimaD.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(x = pop, y = Tajima_D, fill=pop)) + 
  geom_boxplot() + 
  scale_color_manual(values=mycolors) +
  xlab("Population") + ylab("Tajima's D") + theme_classic(base_size=20) +
  guides(fill="none")
dev.off()

# Graph Tajima's D - density
pdf("output/figures/pop_structure/diversity_stats/npstat/Pi_npstat_density_TajimaD.pdf", width = 8, height = 5)
ggplot(data = npstat, aes(Tajima_D, color=pop)) + 
  geom_density() + 
  scale_color_manual(values=mycolors) +
  xlab("Tajima's D") + ylab("Density") + theme_classic(base_size=20) +
  geom_vline(xintercept=0, linetype="dashed") + theme(legend.position="none")
dev.off()

# ================================================================================== #

# Generate summary table

# Create summary table for each population
npstat.summary <- npstat %>% 
  dplyr::group_by(pop) %>% 
  dplyr::summarise(pi_mean = mean(Pi, na.rm = TRUE), pi_median = median(Pi, na.rm = TRUE), 
  Tajima_D_mean = mean(Tajima_D, na.rm = TRUE), Tajima_D_median = median(Tajima_D, na.rm = TRUE))

# Write table
write.table(npstat.summary, "output/tables/npstat.summary.txt", sep='\t')

# ================================================================================== #

# Model statistics

# Add region column to dataframe
npstat$region <- ifelse(npstat$pop == "FC"| npstat$pop == "SLR" | npstat$pop == "SH" | npstat$pop == "ARA" | npstat$pop == "CBL" | 
npstat$pop == "PSG" | npstat$pop == "STC" | npstat$pop == "KH" | npstat$pop == "VD" | npstat$pop == "FR" | npstat$pop == "BMR" | 
npstat$pop == "PGP", "N", "S")

# Model 
npstat.mod <- lmer(Tajima_D ~ region + (1|pop), npstat)

# Check assumptions
# qq plot ---- DOESN"T FIT ASSUMPTIONS
pdf("output/figures/pop_structure/diversity_stats/npstat/Tajima_D_npstat_qqplot.pdf", width = 5, height = 5)
qqPlot(resid(npstat.mod))
dev.off()
# approximate S/L plot
#pdf("output/figures/pop_structure/diversity_stats/npstat/Tajima_D_npstat_SL.pdf", width = 5, height = 5)
#plot(sqrt(abs(resid(npstat.mod, scaled=T)))~fitted(npstat.mod))
#lines(fitted(loess(sqrt(abs(resid(npstat.mod, scaled=T)))~fitted(npstat.mod))))
#dev.off()
# qqplot of random effects
#qqPlot(ranef(npstat.mod)$Region[,1])
shapiro.test(resid(npstat.mod))

# Anova 
anova(npstat.mod)
#Type III Analysis of Variance Table with Satterthwaite's method
#       Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#region 14.342  14.342     1 17.002  27.023 7.248e-05 ***


# ================================================================================== #
# ================================================================================== #

# Analyze just one population at a time (e.g., PB)

# Rename columns
colnames(npstat.PB) <- c("window", "length", "length_outgroup", "read_depth", "S", "Watterson", "Pi", "Tajima_D", "var_S", 
"var_Watterson", "unnorm_FayWu_H", "FayWu_H", "div", "nonsyn_pol", "syn_pol", "nonsyn_div", "syn_div", "alpha", "chr_name")

# Split chr_name column into multiple columns (the first column contains the scaffold name)
npstat.PB.chr.names <- data.frame(do.call('rbind', strsplit(as.character(npstat.PB$chr_name),'.',fixed=TRUE)))
# Remove last 2 columns (they just say "pileup" "stats")
npstat.PB.chr.names <- npstat.PB.chr.names[,1:2]
# Rename columsn
colnames(npstat.PB.chr.names) <- c("chr", "pop")
# Join data frames
npstat.PB.chr <- cbind(npstat.PB, npstat.PB.chr.names)

# ================================================================================== #

# Plot sliding window of Pi for PB
pdf("output/figures/pop_structure/Pi_window_npstat_PB.pdf", width = 5, height = 5)
ggplot(data = npstat.PB.chr, aes(x = chr, y = Pi)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Pi") + theme_classic() 
dev.off()

# Plot sliding window of Watterson for PB
pdf("output/figures/pop_structure/Watterson_window_npstat_PB.pdf", width = 5, height = 5)
ggplot(data = npstat.PB.chr, aes(x = chr, y = Watterson)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Watterson") + theme_classic() 
dev.off()

# Plot sliding window of Tajima D for PB
pdf("output/figures/pop_structure/Tajimas_d_window_npstat_PB.pdf", width = 5, height = 5)
ggplot(data = npstat.PB.chr, aes(x = chr, y = Tajima_D)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Tajima's d") + theme_classic() 
dev.off()
