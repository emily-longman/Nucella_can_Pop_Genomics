# Graph Bio-Oracle data (https://www.bio-oracle.org/index.php)
# Note: prior to running the R script, need to load GDAL module on the VACC
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

#dir(find_root_file(criterion = has_file("README.md")))
#processed_data_path_from_root <- find_root_file("data", "processed", "GEA", "data", criterion = has_file("README.md"))
# Set working directory as path from root
#setwd(processed_data_path_from_root)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'terra', 'raster'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(terra)
library(raster)

# ================================================================================== #

# Read in Bio-oracle data
bio_oracle <- read.csv("data/processed/GEA/data/bio_oracle.csv", header=T)

# Extract the most recent present data (i.e., 2010-01-01T00:00:00Z, which represents 2010-2020)
bio_oracle_present <- bio_oracle %>% 
  filter(time == "2010-01-01T00:00:00Z")

# ================================================================================== #

# Specify parameters

# Set geographic constraints
latitude_range <- c(34, 45)
longitude_range <- c(-125, -120)

# Set study extent
study_extent <- extent(longitude_range[1], longitude_range[2], latitude_range[1], latitude_range[2])

# Raster resolution
raster_resolution <- 0.05

# Create raster layer object
study_raster <- raster(study_extent, res=raster_resolution)

# Specify coordinate reference system
crs(study_raster) <- CRS("+proj=longlat +datum=WGS84")

# Graph empty template
values(study_raster) <- NA
pdf("output/figures/GEA/Raster_template.pdf", width = 5, height = 5)
plot(study_raster, main = "Raster Template for West Coast")
dev.off()

# ================================================================================== #

# Test graphing with just one environmental variable

# Set coordinates 
coordinates <- cbind(bio_oracle_present$longitude, bio_oracle_present$latitude)

# Extract one variable - temperature 
bio_oracle_present_temp_mean <- bio_oracle_present %>% 
  dplyr::select(longitude, latitude, thetao_mean)

# Rasterize data
temp_raster <- rasterize(coordinates, study_raster, bio_oracle_present_temp_mean$thetao_mean, fun = mean, na.rm = TRUE)

# Graph temperature raster
pdf("output/figures/GEA/Test_Raster_bio-oracle_temp_mean.pdf", width = 5, height = 5)
plot(temp_raster, main = "Rasterized Bio-ORACLE Data (thetao_mean)",
     xlim = c(-125, -120), ylim = c(34, 45))
dev.off()

# Write raster tif file
writeRaster(temp_raster, filename = "output/figures/GEA/Test_biooracle_thetao_mean_test_raster.tif", format = "GTiff")

# Look at structure
str(temp_raster)

# Graph with ggplot
# Change to data frame
raster_df <- as.data.frame(temp_raster, xy = TRUE, na.rm = TRUE)
# Graph 
pdf("output/figures/GEA/Test_Raster_bio-oracle_temp_mean_ggplot.pdf", width = 5, height = 5)
ggplot(raster_df, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  scale_fill_viridis_c() +  # Color scale for the raster values
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  labs(title = "Rasterized Bio-ORACLE (thetao_mean)", x = "Longitude", y = "Latitude") +
  theme_void()
dev.off()

# ================================================================================== #

# Environmental variables
env_variables <- c("thetao_max",   "thetao_min",   "thetao_range", "thetao_mean", "chl_mean", "o2_mean", "ph_min", "ph_mean", "so_mean")  # List all env variables here

# Loop through environmental variables and rasterize each dataset
for (var in env_variables) {
  
  # Ensure that each env variable column exists
  if (var %in% colnames(bio_oracle_present)) {
    
    # Select the relevant column from bio_oracle_sites_present
    env_data <- bio_oracle_present %>% dplyr::select(longitude, latitude, all_of(var))
    
    # Convert the data to matrix of coordinates and values
    coords <- cbind(env_data$longitude, env_data$latitude)
    
    # Rasterize the data onto the study_raster template
    env_raster <- rasterize(coords, study_raster, env_data[[var]], fun = mean, na.rm = TRUE)
    common_zlim <- c(10, 26)
    
    # Plot the raster for the current environmental variable
    pdf(paste("output/figures/GEA/Raster_bio-oracle_",var,".pdf", sep = ""), width = 5, height = 5)
    plot(env_raster,
         xlim = c(-125, -112), ylim = c(30, 45), 
         col = colorRampPalette(c("darkblue","lightblue", "yellow","orange", "red"))(100)) 
    dev.off()
 
    # Save the raster to a file
    raster_filename <- paste0("data/processed/GEA/data/tif_files/bio-oracle_", var, "_raster.tif")
    writeRaster(env_raster, filename = raster_filename, format = "GTiff", overwrite = TRUE)
  } else {
    warning(paste("Environmental variable", var, "not found in biooracle_data"))
  }
}

# ================================================================================== #

# Save full path of tif files as envtif
env_tif <- list.files("data/processed/GEA/data/tif_files/", pattern ="bio-oracle_", full.names = TRUE)
env_tif

# Stack tif files
env_stack <- stack(env_tif)
env_stack

# ================================================================================== #

# Convert the raster stack to a dataframe
env_values <- as.data.frame(raster::extract(env_stack, 1:ncell(env_stack), df = TRUE))
env_values$cell <- 1:nrow(env_values)

# ================================================================================== #

# Save the raster stack
writeRaster(env_stack, "data/processed/GEA/data/tif_files/bio-oracle_env_trns_output.tif", format = "GTiff")