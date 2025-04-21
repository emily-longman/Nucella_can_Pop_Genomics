# Download Bio-Oracle data (https://www.bio-oracle.org/index.php )
# Use R package biooracler to access the Bio-Oracle dataset via ERDDAP server (https://github.com/bio-oracle/biooracler)

# Note: need to run this script on my computer rather than the VACC becuase incompatibilities with R version and some of the programs

# Clear memory
rm(list=ls()) 

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
library(tidyverse)
#install.packages('devtools')
library(devtools) 
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)

# ================================================================================== #

# Example online (https://github.com/bio-oracle/biooracler)

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

# Specify output directory
dir <- "/Users/emilylongman/Documents/GitHub/Nucella_can_Pop_Genomics/data/raw/Bio-oracle"

# Set dataset IDs
SST_Present <-"thetao_baseline_2000_2019_depthsurf"
pH_Present<-"ph_baseline_2000_2018_depthsurf"
Chl_Present <-"chl_baseline_2000_2018_depthsurf"
O2_Present<-"o2_baseline_2000_2018_depthsurf"
sal_Present<-"so_baseline_2000_2019_depthsurf"


# Specify datasets and summary statistics to download 
datasets <- list(
  list(dataset_id = SST_Present,
       variables = c("thetao_max","thetao_min", "thetao_range", "thetao_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = pH_Present,
       variables = c("ph_min", "ph_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = Chl_Present,
       variables = c("chl_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = O2_Present,
       variables = c("o2_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir),
  list(dataset_id = sal_Present,
       variables = c("so_mean"),
       constraints = constraints,
       fmt = "csv", directory = dir)
)


# Download datasets using for loop
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
    message("Rename ", new_file, " to ", new_name)
  } else {
    message("No new file found for ", dataset_id, " or multiple new files detected.")
  }
}

# Check list of files
list.files(dir)

# Note: upload these files to the appropriate folder on the VACC
