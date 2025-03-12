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
install.packages(c('ggplot2', 'RColorBrewer'))
library(ggplot2)
library(RColorBrewer)

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

# ================================================================================== #


# Split chr_name column into multiple columns (the first column contains the scaffold name)
npstat.ARA.chr.names <- data.frame(do.call('rbind', strsplit(as.character(npstat.ARA$chr_name),'.',fixed=TRUE)))
# Remove last 2 columns (they just say "pileup" "stats")
npstat.ARA.chr.names <- npstat.ARA.chr.names[,1:2]
# Rename columsn
colnames(npstat.ARA.chr.names) <- c("chr", "pop")
# Join data frames
npstat.ARA.chr <- cbind(npstat.ARA, npstat.ARA.chr.names)




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




# Rename columns
colnames(npstat.SBR) <- c("window", "length", "length_outgroup", "read_depth", "S", "Watterson", "Pi", "Tajima_D", "var_S", 
"var_Watterson", "unnorm_FayWu_H", "FayWu_H", "div", "nonsyn_pol", "syn_pol", "nonsyn_div", "syn_div", "alpha", "chr_name")

# Split chr_name column into multiple columns (the first column contains the scaffold name)
npstat.SBR.chr.names <- data.frame(do.call('rbind', strsplit(as.character(npstat.SBR$chr_name),'.',fixed=TRUE)))
# Remove last 2 columns (they just say "pileup" "stats")
npstat.SBR.chr.names <- npstat.SBR.chr.names[,1:2]
# Rename columsn
colnames(npstat.SBR.chr.names) <- c("chr", "pop")
# Join data frames
npstat.SBR.chr <- cbind(npstat.SBR, npstat.SBR.chr.names)

# ================================================================================== #

# Plot sliding window of Pi for ARA
pdf("output/figures/pop_structure/Pi_window_npstat_ARA.pdf", width = 5, height = 5)
ggplot(data = npstat.ARA.chr, aes(x = chr, y = Pi)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Pi") + theme_classic() 
dev.off()

# Plot sliding window of Pi for PB
pdf("output/figures/pop_structure/Pi_window_npstat_PB.pdf", width = 5, height = 5)
ggplot(data = npstat.PB.chr, aes(x = chr, y = Pi)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Pi") + theme_classic() 
dev.off()

# Plot sliding window of Pi for SBR
pdf("output/figures/pop_structure/Pi_window_npstat_SBR.pdf", width = 5, height = 5)
ggplot(data = npstat.SBR.chr, aes(x = chr, y = Pi)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Pi") + theme_classic() 
dev.off()

# ================================================================================== #

# Plot sliding window of Watterson for ARA
pdf("output/figures/pop_structure/Watterson_window_npstat_ARA.pdf", width = 5, height = 5)
ggplot(data = npstat.ARA.chr, aes(x = chr, y = Watterson)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Watterson") + theme_classic() 
dev.off()

# Plot sliding window of Watterson for PB
pdf("output/figures/pop_structure/Watterson_window_npstat_PB.pdf", width = 5, height = 5)
ggplot(data = npstat.PB.chr, aes(x = chr, y = Watterson)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Watterson") + theme_classic() 
dev.off()

# Plot sliding window of Watterson for SBR
pdf("output/figures/pop_structure/Watterson_window_npstat_SBR.pdf", width = 5, height = 5)
ggplot(data = npstat.SBR.chr, aes(x = chr, y = Watterson)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Watterson") + theme_classic() 
dev.off()

# ================================================================================== #

# Plot sliding window of Tajima D for ARA
pdf("output/figures/pop_structure/Tajimas_d_window_npstat_ARA.pdf", width = 5, height = 5)
ggplot(data = npstat.ARA.chr, aes(x = chr, y = Tajima_D)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Tajima's d") + theme_classic() 
dev.off()

# Plot sliding window of Tajima D for PB
pdf("output/figures/pop_structure/Tajimas_d_window_npstat_PB.pdf", width = 5, height = 5)
ggplot(data = npstat.PB.chr, aes(x = chr, y = Tajima_D)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Tajima's d") + theme_classic() 
dev.off()

# Plot sliding window of Tajima D for SBR
pdf("output/figures/pop_structure/Tajimas_d_window_npstat_SBR.pdf", width = 5, height = 5)
ggplot(data = npstat.SBR.chr, aes(x = chr, y = Tajima_D)) + 
  geom_point(shape = 21, size = 1) + 
  xlab("Chr") + ylab("Tajima's d") + theme_classic() 
dev.off()