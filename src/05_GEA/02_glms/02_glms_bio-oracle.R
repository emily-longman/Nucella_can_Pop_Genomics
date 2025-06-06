# Generalized linear models to assess relationship between SNPs and environmental variables

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
#install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
# install.packages(c('data.table', 'tidyverse', 'foreach', 'poolfstat', 'magrittr', 'reshape2', 'broom', 'SNPRelate'))
library(data.table)
library(tidyverse)
library(foreach)
library(poolfstat)
library(magrittr)
library(reshape2)
library(broom)
library(stats)

# Install and load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
y = as.numeric(args[1]) # Chunk: this is the chunk of 1000 windows  (total of 32 chunks)
k = as.numeric(args[2]) # Array: this is a specific window within a given chunk (1000 windows in each chunk)

# ================================================================================== #

# Load data

# Load SNPs of interest (baypass POD outlier SNPs - 320,742 SNPs)
baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs", header=T)

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Reload windows (generated on previous script - 01_windows.R)
load("data/processed/GEA/glms/glms_windows.RData")

# ================================================================================== #

# Format data

# Rename "location" column as "sampleId"
names(bio_oracle_sites_2010)[names(bio_oracle_sites_2010) == "location"] <- "sampleId"

# Create SNP_id column for outlier SNP list
baypass_POD_sig_SNPs <- baypass_POD_sig_SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# Filter for only outlier SNPs (NOTE: somehow getting one less SNP than in baypass_POD_sig_SNPs)
snp.dt.sig <- snp.dt %>% filter(snp.dt$SNP_id %in% baypass_POD_sig_SNPs$SNP_id)

# ================================================================================== #
# ================================================================================== #

# Split windows up into chunks

# Define window and step size
win.bp = 5e4
step.bp = win.bp+1

# Specify chunk size
chunk_size <- 1000

# Create list with chunks of 1000 windows
split(wins$i, ceiling(seq_along(wins$i) / chunk_size)) -> subdivision_list

# Check number of elements/chunks in list (total number of chunks: 32)
length(subdivision_list)

# ================================================================================== #
# ================================================================================== #

# Get data for a specified window
# To do so, specify a given chunk ('y'), then a given window ('k' array)
message(paste("I am doing chunk number y =", y, "and array number k =", k,  sep = " "))

# Get window list for chunk y
chunk_of_choice <- subdivision_list[[y]]

# Check number of elements/windows in chunk
message(paste("Chunk", y, "has", length(chunk_of_choice), "windows",  sep = " "))

# Get window information for chunk 'y'
wins %>%
filter(i %in% chunk_of_choice) -> 
wins.y

# Filter snp.dt for a given window "k" in chunk "y"
snp.dt %>%
filter(chr == wins.y$chr[k]) %>%
filter(pos > wins.y$start[k] & pos < wins.y$end[k]) -> 
data_win

# ================================================================================== #

model.output =
  # For each SNP in a given window extract allele freq then run GLM with bio-oracle data
  foreach(i=1:dim(data_win)[1], .combine = "rbind")%do%{
    
    seqResetFilter(genofile)

    ###############################################################

    # Calculate allele frequency for SNP i in window k in chunk y
    seqSetFilter(genofile, variant.id=data_win$variant.id[i], verbose=T)

      # Extract allele depth
      ad_i <- seqGetData(genofile, "annotation/format/AD")
      # Extract total depth
      dp_i <- seqGetData(genofile, "annotation/format/DP")

        # Create af data table with ad, dp, sample ID and variant id
      af_i <- data.table(ad=expand.grid(ad_i$data)[,1], dp=dp_i[,1], 
      sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad_i$data)[2]),
      variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad_i$data)[1]))
      
      # Merge allele freq table af_i and snp.dt
      af_i_snp <- merge(af_i, snp.dt, by="variant.id")
      
      # Calculate allele frequency ('af')
      af_i_snp[,af:=ad/dp]
      # Calculate the mean effective coverage ('nEff') (note: each pool consists of 20 dogwhelks)
      af_i_snp[,nEff:=round((dp*2*20 - 1)/(dp+2*20))]
      # Calculate the effective read-depth
      af_i_snp[,af_nEff:=round(af*nEff)/nEff]
      
    ###############################################################

    # Join with bio-oracle environmental data
    left_join(af_i_snp, bio_oracle_sites_2010, by ="sampleId") -> af_i_snp_enviro
    
    # Create long format data table with the enviro data in column "value" and the specific variable identified in column "column"
    af_i_snp_enviro %>% as_tibble %>% gather(key = "enviro_var", value = "value", `thetao_max`:`so_mean`) ->
      gathered_data
    
    # Get names of environmental variables
    unique(gathered_data$enviro_var) -> enviro_vars_names
    
    ###############################################################

    # Null t0 model
    t0 <- lm(af~1, data=af_i_snp_enviro)

    # Run model for each variable
    real_estimates =
      foreach(j=enviro_vars_names, .combine = "rbind", .errorhandling = "remove")%do%{
        
        # Extract data for 'j' environmental variable
        gathered_data %>% filter(enviro_var == j) -> inner.tmp
        
        # Model allele freq for 'j' enviro variable
        t1 <- lm(af~(value), data=inner.tmp)
        t1.sum <- (summary(t1))
        # Likelihood ratio test (LRT) between enviro model and null model
        LRT.aov <- anova(t0, t1, test="Chisq")
        
        # Generate output table with model comparison information for each variable
        data.frame(
          test_code = 0,
          chr = unique(af_i_snp$chr),
          pos = unique(af_i_snp$pos),
          variable = j,
          missing=seqMissing(genofile),
          data = "real",
          lrt_p = LRT.aov$`Pr(>Chi)`[2],
          AIC_enviro = AIC(t1),
          AIC_null = AIC(t0),
          Beta = last(t1$coef))
          
      } # End run for all j enviro var

    ###############################################################

    # Permutations to generate null expectation of association between af and enviro var
    permutation_estimates = 
      foreach(j=enviro_vars_names_sub, .combine = "rbind", .errorhandling = "remove")%do%{
                
        # Extract data for 'j' environmental variable
        gathered_data %>% filter(enviro_var == j) -> inner.tmp
                
                # Do 100 permutations; 
                foreach(l=1:100, .combine = "rbind")%do%{
                  message(paste("Permutation number:", l))
                  
                  # Shuffle enviro data for 'j' enviro variable
                  inner.tmp %>% mutate(suffle_value = sample(value)) -> inner.tmp.shuffle
                  
                  # Model allele freq for 'j' enviro variable with shuffled data
                  t1 <- lm(af~(suffle_value), data=inner.tmp.shuffle)
                  t1.sum <- (summary(t1))
                  LRT.aov <- anova(t0, t1, test="Chisq")
                  
                  # Generate output table with model comparison information for each variable for the shuffled data
                  data.frame(
                    test_code = l,
                    chr = unique(af_i_snp$chr),
                    pos = unique(af_i_snp$pos),
                    variable = j,
                    missing = seqMissing(genofile),
                    data = "permutation",
                    lrt_p = LRT.aov$`Pr(>Chi)`[2],
                    AIC_wea = AIC(t1),
                    AIC_null = AIC(t0),
                    Beta = last(t1$coef))
                  
                } # End for i permutation
              } # End run for all j enviro var

    # Combine real estimates and permutations
    rbind(real_estimates, permutation_estimates) -> all_data
    return(all_data)
  }

# ================================================================================== #

# Generate folders and save output

# Make folder for model outputs
folder_name <- paste("data/processed/GEA/glms/glms_window_analysis")
system(paste("mkdir", folder_name, sep = " "))

# Make folder for chunk y
folder_name <- paste("data/processed/GEA/glms/glms_window_analysis/LM_100perm_Bio-Oracle_window_", y, sep = "")
system(paste("mkdir", folder_name, sep = " "))

# Save file for window y
file_name <- paste("LM_100perm_Bio-Oracle_y", y, "chr", wins$chr[k], "start", wins$start[k], "stop", wins$end[k], sep = "_")
save(model.output.test, file = paste(folder_name, "/" , file_name, ".Rdata", sep = "") )
