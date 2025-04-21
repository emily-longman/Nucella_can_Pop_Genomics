# Format Bio-Oracle data (https://www.bio-oracle.org/index.php )

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
install.packages(c('purr', 'ggplot2', 'RColorBrewer'))
library(purr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Load site metadata
meta <- read.csv("data/processed/GEA/guide_files/Populations_metadata.csv", header=T)

# Read in Bio-oracle data
sst <- read.csv("data/raw/Bio-oracle/thetao_baseline_2000_2019_depthsurf.csv", header=T)
chl <- read.csv("data/raw/Bio-oracle/chl_baseline_2000_2018_depthsurf.csv", header=T)
o2 <- read.csv("data/raw/Bio-oracle/o2_baseline_2000_2018_depthsurf.csv", header=T)
ph <- read.csv("data/raw/Bio-oracle/ph_baseline_2000_2018_depthsurf.csv", header=T)
so <- read.csv("data/raw/Bio-oracle/so_baseline_2000_2019_depthsurf.csv", header=T)

# ================================================================================== #

# Population names
pops <- meta$Site
sites <- cbind("longitude" = meta$long, "latitude" =  meta$Lat)

# Combine datasets into one large dataset
bio_oracle <- reduce(list(sst, chl, o2, ph, so), full_join, by = c("time", "latitude", "longitude"))


# Extract site specific data
sites_environ <- data.frame(Name=sites, extract(bio_oracle, sites))

allpresentdata <- combined_data %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    location = case_when(
      latitude == 43.304 & longitude == -124.402 ~ "ARA",
      latitude == 38.319 & longitude == -123.074 ~ "BMR",
      latitude == 42.841 & longitude == -124.565 ~ "CBL",
      latitude == 44.838 & longitude == -124.059 ~ "FC",
      latitude == 38.512 & longitude == -123.255 ~ "FR",
      latitude == 35.290 & longitude == -120.884 ~ "HZD",
      latitude == 39.605 & longitude == -123.789 ~ "KH",
      latitude == 34.881 & longitude == -120.640 ~ "OCT",
      latitude == 35.665 & longitude == -121.287 ~ "PB",
      latitude == 37.185 & longitude == -122.398 ~ "PGP",
      latitude == 36.519 & longitude == -121.954 ~ "PL",
      latitude == 41.771 & longitude == -124.253 ~ "PSG",
      latitude == 35.729 & longitude == -121.319 ~ "PSN",
      latitude == 36.448 & longitude == -121.929 ~ "SBR",
      latitude == 44.250 & longitude == -124.115 ~ "SH",
      latitude == 44.505 & longitude == -124.085 ~ "SLR",
      latitude == 40.030 & longitude == -124.081 ~ "STC",
      latitude == 34.730 & longitude == -120.616 ~ "STR",
      latitude == 39.281 & longitude == -123.804 ~ "VD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(location))

allpresentdata

