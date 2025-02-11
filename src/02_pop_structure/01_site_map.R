# Generate Site Map

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file("output",  criterion = has_file("README.md")))
output_path_from_root <- find_root_file("output", criterion = has_file("README.md"))
# Set working directory as path from root
setwd(output_path_from_root)

# ================================================================================== #

# Load packages
install.packages(c('ggplot2', 'maps', 'mapdata', 'ggrepel'))
library(ggplot2)
library(maps) #Contains a lot of outlines of continents, countries, states, and counties
library(mapdata) #Supplement to maps package
library(ggrepel)
#library(ggsn) 

# ================================================================================== #

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon")) 

# ================================================================================== #

# Coordinates of sites 
sites <- data.frame(
  longitude = c(-123.074, -123.25506, -122.397576, -121.286767, -121.928989, -121.95367, -123.803572,
                -124.114767, -124.084799, -124.0593, -124.401546, -124.564709, -124.252929, 
                -124.080911, -123.78945, -121.31867, -120.63994, -120.61569, -120.88384),
  latitude = c(38.319, 38.51198, 37.185064, 35.66549, 36.447503, 36.51939, 39.280903, 44.249988, 
               44.505402, 44.83777, 43.304022, 42.840972, 41.771209, 40.030114, 39.604613, 35.72893,
               34.88117, 34.73024, 35.28994),
  site = c("BMR", "FR", "PGP", "PB", "SBR", "PL", "VD", "SH", "SLR", "FC", "ARA", "CBL", "PSG", "STC",
           "KH", "PSN", "OCT", "STR", "HZD"), 
  site.full = c("Bodega Marine Reserve", "Fort Ross", "Pigeon Point", "Piedras Blancas", 
                "Soberanes Point", "Point Lobos", "Van Damme", "Strawberry Hill", 
                "Seal Rock", "Fogarty Creek", "Cape Arago", "Cape Blanco", "Point Saint George", 
                "Shelter Cove", "Kibesillah Hill", "Point Sierra Nevada", "Occulto", "Stairs", "Hazards")
  )

# Coordinates of site labels
long.site.labels <- c(-123.074-0.75, -123.25506-0.6, -122.397576-0.75, -121.286767-0.75, 
                      -121.928989-0.75, -121.95367-0.5, -123.803572-0.75, -124.114767-0.75,
                      -124.084799-0.75, -124.0593-0.7, -124.401546-0.75, -124.564709-0.75,
                      -124.252929-0.75, -124.080911-0.75, -123.78945-0.5, -121.31867-0.75, 
                      -120.63994-0.75, -120.61569-0.6, -120.88384-0.7)
long.site.labels.full <- c(-123.074-2.75, -123.25506-1.6, -122.397576-1.75, -121.286767-1.75, 
                      -121.928989-1.9, -121.95367-1.5, -123.803572-1.75, -124.114767-1.75,
                      -124.084799-1.75, -124.0593-1.7, -124.401546-1.5, -124.564709-1.5,
                      -124.252929-2.25, -124.080911-1.75, -123.78945-1.5, -121.31867-2.75, 
                      -120.63994-1.2, -120.61569-1.2, -120.88384-1.2)
lat.site.labels <- c(38.319-0.05, 38.51198+0.04, 37.185064, 35.66549-0.07, 36.447503-0.09, 36.51939+0.1,
                     39.280903, 44.249988-0.04, 44.505402, 44.83777+0.04, 43.304022, 42.840972, 
                     41.771209, 40.030114, 39.604613-0.06, 35.72893+0.12, 34.88117, 
                     34.73024-0.12, 35.2899409-0.02)

# ================================================================================== #

# Create site map (site codes)
pdf("figures/Site_map.pdf", width = 10, height = 10)
g<- ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
             shape = 20, fill = "black") + coord_fixed(1.3) +
  geom_text(data=sites, aes(long.site.labels, lat.site.labels, label=site))+ 
  xlim(c(-134, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() 
g + annotate(geom = 'text', size = 5,
         x = -127, y = 34.2, label = 'Pacific Ocean', fontface = 'italic')
dev.off()

# Create site map (full site names)
pdf("figures/Site_map_full_names.pdf", width = 10, height = 10)
g<- ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
             shape = 20, fill = "black") + coord_fixed(1.3) +
  geom_text(data=sites, aes(long.site.labels.full, lat.site.labels, label=site.full))+ 
  xlim(c(-134, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic()
g + annotate(geom = 'text', size = 5,
         x = -127, y = 34.2, label = 'Pacific Ocean', fontface = 'italic')
dev.off()
