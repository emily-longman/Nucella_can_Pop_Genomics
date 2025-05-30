# Generalized linear models to assess relationship between SNPs and environmental variables

# Clear memory
rm(list=ls())

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
#library(SNPRelate) # not available for this version of R
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
k = as.numeric(args[1]) ### this is the array argument -- this is 1-1900
y = as.numeric(args[2]) ### this is the chunk argument  1-X

# ================================================================================== #

# Load data

# Load pooldata object
load("data/processed/pop_structure/pooldata.RData")
# Load SNPs of interest (baypass rnp outlier SNPs: p=0.001)
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/SNPs.Interest.pval.001.RData")
# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)
# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# ================================================================================== #

# Create SNP_id column for outlier SNP list
SNPs.Interest.pval.001 <- SNPs.Interest.pval.001 %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Calculate allele frequency

# Turn the ref allele read count data object into a dataframe and add SNP name (Chr and pos)
readcount <- as.data.frame(pooldata@refallele.readcount)
readcov <- as.data.frame(pooldata@readcoverage)
allele_freq <- as.data.frame(readcount/readcov)

# Set rownames as sites
colnames(allele_freq) <- pooldata@poolnames

# Extract SNP info and set as rownames
snp_info <- as.data.frame(pooldata@snp.info)
# Create SNP_id column
snp_info_full <- snp_info %>% mutate(SNP_id = paste(Chromosome, Position, sep = "_"))

# Set SNP_id as rownames
rownames(allele_freq) <- snp_info_full$SNP_id
allele_freq$SNP_id <- snp_info_full$SNP_id

# ================================================================================== #

# Filter for only outlier SNPs (total: 70,205)

# Filter allele freq for only outlier SNPs (NOTE: somehow getting only 61,897 SNPs)
allele_freq_sig <- allele_freq %>% filter(rownames(allele_freq) %in% SNPs.Interest.pval.001$SNP_id)


# ================================================================================== #
# ================================================================================== #

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# Filter for only outlier SNPs (NOTE: somehow getting only 61,896 SNPs rather than 70,205)
snp.dt.sig <- snp.dt %>% filter(snp.dt$SNP_id %in% SNPs.Interest.pval.001$SNP_id)

# ================================================================================== #
# ================================================================================== #

# Define window and step size
win.bp = 5e4
step.bp = win.bp+1

# Generate windows (if a contig is less that the specified window size, the if statement, includes it)
wins <- foreach(chr.i=unique(snp.dt$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)

                  tmp <-  snp.dt %>%
                    filter(chr == chr.i)

                    if(max(tmp$pos) <= step.bp){
                      data.table(chr=chr.i, start=min(tmp$pos), end=max(tmp$pos))
                    } else if(max(tmp$pos) > step.bp){
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                } }
# Add column with window number
wins[,i:=1:dim(wins)[1]]

# Save windows
save(wins, file="data/processed/GEA/glms/glms_windows.RData")
# Reload windows
load("data/processed/GEA/glms/glms_windows.RData")

# ================================================================================== #

# Split windows up into chunks

# Specify chunk size
chunk_size = 1000

# Create list with chunks of 1000 windows
split(wins$i, ceiling(seq_along(wins$i) / chunk_size)) -> subdivision_list

# Check number of elements/chunks in list (total number of chunks: 32)
length(subdivision_list)

# ================================================================================== #

# Get data for a specified window
# To do so, specify a given chunk ('y'), then a given window ('k' array)
message(paste("I am doing chunk number y = ", y, " and array number k = ", k,  sep = " "))

# Get window list for chunk y
chunk_of_choice = subdivision_list[[y]]

# Get window information for chunk 'y'
wins %>%
filter(i %in% chunk_of_choice) -> 
wins.y

# Filter snp.dt for a given window, specified as array "k" in window "y"
snp.dt %>%
filter(chr == wins.y$chr[k]) %>%
filter(pos > wins.y$start[k] & pos < wins.y$end[k]) -> 
data_win

# ================================================================================== #

# Calculate allele frequencies

# Extract allele depth
ad <- seqGetData(genofile, "annotation/format/AD")
# Extract total depth
dp <- seqGetData(genofile, "annotation/format/DP")

# Change dp to list
if(class(dp)[1]!="SeqVarDataList") {
          dp.list <- list()
          dp.list$data <- dp
          dp <- dp.list
        }

# Create data table with ad, dp, sample ID and variant id (NOT WORKING PROPERLY!)
af <- data.table(ad=expand.grid(ad$data)[,1],
                 dp=expand.grid(dp$data)[,1],
                 sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                 variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
#Warning messages:
#1: In as.data.table.list(x, keep.rownames = keep.rownames, check.names = check.names,  :
#  Item 2 has 283051892 rows but longest item has 570077748; recycled with remainder.

# Merge af and snp.dt
af.i <- merge(af, snp.dt, by="variant.id")

# Calculate allele frequency
af.i[,af:=ad/dp]

# Calculate effective read-depth (note: each pool consists of 20 dogwhelks)
af.is <- merge(af.i, samps, by="sampleId") #STOPPED HERE - DON"T NEED SAMPLE info bc poolsizes the same

af.is[chr!="X", nEff:=round((dp*2*20 - 1)/(dp+2*20))]
af.is[,af_nEff:=round(af*nEff)/nEff]

### return
af.is[,c("sampleId", "af_nEff", "nEff"), with=F]


# ================================================================================== #

model.output =
  # For each SNP in a given window extract 
  foreach(i=1:dim(data_win)[1], .combine = "rbind")%do%{
    
    seqResetFilter(genofile)


    ######

    # genotype info
    seqSetFilter(genofile, variant.id=data_win$variant.id[i], verbose=T)
    seqGetData(genofile, "genotype") -> genotypes_i
    seqGetData(genofile, "sample.id") -> sampsids_raw
    gsub("/users/c/p/cpetak/WGS/BWA_out/", "", sampsids_raw) -> sampsids_f1
    gsub("_200925_A00421_0244_AHKML5DSXY_", "_", sampsids_f1) -> sampsids_f2
    gsub("_L002_R1_001.rmdup.bam", "", sampsids_f2) -> sampsids_f3
    gsub("18170X", "", sampsids_f3) -> sampsids_clean
    
    colnames(genotypes_i) = sampsids_clean
    
    genotypes_i %>% colSums() %>%
      data.frame(genotype = .) %>%
      mutate(SampleId = row.names(.)) %>%
      separate(SampleId, into = c("id","numid","Sid"), sep = "_") ->
      genotable_v1
    
    genotable_v1 %>%
      group_by(id) %>%
      summarise(AF = mean(genotype, na.rm = T),
                MissDat = sum(is.na(genotype))
      ) %>%
      mutate(chr =seqGetData(genofile, "chromosome"),
             pos=seqGetData(genofile, "position"),
             variant.id=seqGetData(genofile, "variant.id"),
             nAlleles=seqNumAllele(genofile),
             missing=seqMissing(genofile),
             allele=seqGetData(genofile, "allele")
      ) ->
      dat_for_analysis
    #####
    left_join(dat_for_analysis, weather_dat, by = "id") -> dat_for_analysis.wea
    
    t0 <- lm(AF~1, data=dat_for_analysis.wea)
    
    dat_for_analysis.wea %>% 
      as_tibble %>% 
      gather(key = "column", value = "value", `thetao_max`:`so_mean`) ->
      gathered_data
    
    unique(gathered_data$column) -> ecovars_names
    
    real_ests =
      foreach(j=ecovars_names, .combine = "rbind", .errorhandling = "remove")%do%{
        
        gathered_data %>%
          filter(column == j) ->
          inner.tmp
        
        t1 <- lm(AF~(value), data=inner.tmp)
        t1.sum <- (summary(t1))
        LRT.aov <- anova(t0, t1, test="Chisq")
        
        data.frame(
          test_code = 0,
          chr = unique(dat_for_analysis$chr),
          pos = unique(dat_for_analysis$pos), ### k - 1-999
          missing=seqMissing(genofile),
          data= "real",
          variable = j,
          lrt_p = LRT.aov$`Pr(>Chi)`[2],
          AIC_wea = AIC(t1),
          AIC_null = AIC(t0),
          Beta = last(t1$coef)
          
        )
      }
    
    perms_dat = 
      foreach(j=ecovars_names, 
              .combine = "rbind", 
              .errorhandling = "remove")%do%{
                
                gathered_data %>%
                  filter(column == j) ->
                  inner.tmp
                
                foreach(l=1:100, .combine = "rbind")%do%{
                  
                  inner.tmp %>%
                    mutate(suffle_value = sample(value)) ->
                    inner.tmp.shuffle
                  
                  t1 <- lm(AF~(suffle_value), data=inner.tmp.shuffle)
                  t1.sum <- (summary(t1))
                  LRT.aov <- anova(t0, t1, test="Chisq")
                  message(l)
                  
                  data.frame(
                    test_code = l,
                    chr = unique(dat_for_analysis$chr),
                    pos = unique(dat_for_analysis$pos), ### k - 1-999
                    variable = j,
                    missing=seqMissing(genofile),
                    data= "perm",
                    lrt_p = LRT.aov$`Pr(>Chi)`[2],
                    AIC_wea = AIC(t1),
                    AIC_null = AIC(t0),
                    Beta = last(t1$coef)
                  )
                  
                  
                  
                } ## j
              } ## i
    rbind(real_ests, perms_dat) -> all_dat 
    return(all_dat)
  }

name_file <- paste("GLMout_multiarray_100perm_allmodels_BioOrac", y, wins$chr[k], wins$start[k] , wins$end[k], sep = "_")
folname = paste("./GLM_outs_100perm_arrayno_allmodels_BioOrac", "_", y, sep = "")
system( paste("mkdir",  folname, sep = " ") )
save(model.output, file = paste( folname, "/" ,name_file, ".Rdata", sep = "") )
