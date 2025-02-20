# Get Bio-Oracle data (https://www.bio-oracle.org/index.php )
# Use R package biooracler to access the Bio-Oracle dataset via ERDDAP server (https://github.com/bio-oracle/biooracler)

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Create new library for package

# Check library paths
#.libPaths() 
# Create new path for Biooracle libraries
#.libPaths(c("/Library/Frameworks/R.framework/Versions/4.1/Resources/library/Rpackages_biooracle", .libPaths()))
# Set library path to biooracle path
#libpath <- .libPaths()[1]

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
raw_data_path_from_root <- find_root_file("data", "raw", "Bio-oracle", criterion = has_file("README.md"))
# Set working directory as path from root
setwd(raw_data_path_from_root)

# ================================================================================== #

# Load packages
install.packages('devtools')
library(devtools) 
devtools::install_github("bio-oracle/biooracler")
library(biooracler)
library(tidyverse)

#install.packages("Rcpp")
#library(Rcpp)
#install.packages('terra')
#library(terra)

#install.packages(c("tidyverse", "sdmpredictors", "raster", "sp", "dismo"))
#library(tidyverse)
#library(sdmpredictors)
#library(raster)
#library(sp)
#library(dismo)

# ================================================================================== #

# See available layers
list_layers()
# See all layers for a specific variable (e.g., temp)
print(list_layers("temp"), n=100)

# Define time, lat, and long  
time = c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
latitude = c(34, 45)
longitude = c(-120, -125)

# Set constraints
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

# Example online (https://github.com/bio-oracle/biooracler)
#dataset_id <- "tas_baseline_2000_2020_depthsurf"
#variables = c("tas_max", "tas_min")
#layers <- download_layers(dataset_id, variables, constraints)
#dir <- "/Users/emilylongman/Documents/GitHub/Nucella_can_Pop_Genomics/data/raw/Bio-oracle"
#download_layers(dataset_id, variables, constraints, fmt = "csv", directory = dir)

# Set dataset IDs
SST_Present <-"thetao_baseline_2000_2019_depthsurf"
Chl_Present <-"chl_baseline_2000_2018_depthsurf"
Fe_Present<-"dfe_baseline_2000_2018_depthsurf"
No3_Present<-"no3_baseline_2000_2018_depthsurf"
O2_Present<-"o2_baseline_2000_2018_depthsurf"
pH_Present<-"ph_baseline_2000_2018_depthsurf"
phyto_Present<-"phyc_baseline_2000_2020_depthsurf"
po4_Present<-"po4_baseline_2000_2018_depthsurf"
sal_Present<-"so_baseline_2000_2019_depthsurf"
si_Present<-"si_baseline_2000_2018_depthmax"
sws_Present<-"sws_baseline_2000_2019_depthsurf"
swd_Present<-"swd_baseline_2000_2019_depthsurf"





datasets <- list(
  list(dataset_id = SST_Present,
       variables = c("thetao_max","thetao_min", "thetao_range", "thetao_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = Chl_Present,
       variables = c("chl_max","chl_min", "chl_range", "chl_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = Fe_Present,
       variables = c("dfe_max","dfe_min", "dfe_range", "dfe_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = No3_Present,
       variables = c("no3_max","no3_min", "no3_range", "no3_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = O2_Present,
       variables = c("o2_max","o2_min", "o2_range", "o2_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = pH_Present,
       variables = c("ph_max","ph_min", "ph_range", "ph_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = phyto_Present,
       variables = c("phyc_max","phyc_min", "phyc_range", "phyc_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = po4_Present,
       variables = c("po4_max","po4_min", "po4_range", "po4_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = sal_Present,
       variables = c("so_max","so_min", "so_range", "so_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = si_Present,
       variables = c("si_max","si_min", "si_range", "si_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = sws_Present,
       variables = c("sws_max","sws_min", "sws_range", "sws_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = swd_Present,
       variables = c("swd_max","swd_min", "swd_range", "swd_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir)
)



# Download datasets 
for (dataset in datasets) {
  
  dataset_id <- dataset$dataset_id
  variables <- dataset$variables
  constraints <- dataset$constraints
  
  # List files before downloading
  files_before <- list.files(dir, full.names = TRUE)
  
  # Download dataset
  download_layers(dataset_id, variables = variables, constraints = constraints, fmt = "csv", directory = dir)
  
  # List files after downloading
  files_after <- list.files(dir, full.names = TRUE)
  
  # Identify the newly downloaded file
  new_file <- setdiff(files_after, files_before)
  
  if (length(new_file) == 1) {
    new_name <- file.path(dir, paste0(dataset_id, ".csv"))
    file.rename(new_file, new_name)
    message("Renamed ", new_file, " to ", new_name)
  } else {
    message("No new file found for ", dataset_id, " or multiple new files detected.")
  }
}



list.files(dir)


files <- list.files(dir, pattern = "*.csv")
files

setwd(dir)
filtered_files <- files[str_detect(files, "baseline")]
filtered_files




temp<-read_csv("thetao_baseline_2000_2019_depthsurf.csv")
chl<-read_csv("chl_baseline_2000_2018_depthsurf.csv")
dfe<-read_csv("dfe_baseline_2000_2018_depthsurf.csv")
no3<-read_csv("no3_baseline_2000_2018_depthsurf.csv")
o2<-read_csv("o2_baseline_2000_2018_depthsurf.csv")
ph<-read_csv("ph_baseline_2000_2018_depthsurf.csv")
phyc<-read_csv("phyc_baseline_2000_2020_depthsurf.csv")
po4<-read_csv("po4_baseline_2000_2018_depthsurf.csv")
si<-read_csv("si_baseline_2000_2018_depthmax.csv")
so<-read_csv("so_baseline_2000_2019_depthsurf.csv")
swd<-read_csv("swd_baseline_2000_2019_depthsurf.csv")
sws<-read_csv("sws_baseline_2000_2019_depthsurf.csv")


temp<-temp[-1,]
chl<-chl[-1,]
dfe<-dfe[-1,]
no3<-no3[-1,]
o2<-o2[-1,]
ph<-ph[-1,]
phyc<-phyc[-1,]
po4<-po4[-1,]
si<-si[-1,]
so<-so[-1,]
swd<-swd[-1,]
sws<-sws[-1,]

head(temp)

combined_data <- reduce(list(temp, chl, dfe, no3, o2, ph, phyc, po4, si, so, swd, sws), 
                        full_join, by = c("time", "latitude", "longitude"))


combined_data



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

