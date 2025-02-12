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
install.packages(c('ggplot2', 'maps', 'mapdata', 'ggrepel', 'RColorBrewer'))
library(ggplot2)
library(maps) 
library(mapdata)
library(ggrepel)
library(RColorBrewer)
#library(ggsn) 

# ================================================================================== #

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon")) 

# ================================================================================== #

# Coordinates of sites 
sites <- data.frame(
  longitude = c(-124.0593, -124.0848, -124.1148, -124.4015, -124.5647, -124.2529, -124.0809, -123.7895, -123.8036, 
                -123.2551, -123.0740, -122.3976, -121.9537, -121.9290, -121.3187, -121.2868, -120.8838, -120.6399, -120.6157),
  latitude = c(44.83777, 44.50540, 44.24999, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461, 39.28090, 38.51198, 38.31900, 
               37.18506, 36.51939, 36.44750, 35.72893, 35.66549, 35.28994, 34.88117, 34.73024),
  site = c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"), 
  site.full = c("Fogarty Creek"," Seal Rock",  "Strawberry Hill", "Cape Arago", "Cape Blanco",  "Point Saint George", "Shelter Cove", 
                "Kibesillah Hill", "Van Damme", "Fort Ross", "Bodega Marine Reserve", "Pigeon Point", "Point Lobos", "Soberanes Point", 
                "Point Sierra Nevada", "Piedras Blancas", "Hazards" , "Occulto", "Stairs")  
  )

# Coordinates of site labels
long.site.labels <- c(-124.0593-0.7, -124.0848-0.75, -124.1148-0.75, -124.4015-0.75, -124.5647-0.75, -124.2529-0.75, 
                      -124.0809-0.75, -123.7895-0.5, -123.8036-0.75, -123.2551-0.6, -123.0740-0.75, -122.3976-0.75, 
                      -121.9537-0.5, -121.9290-0.75, -121.3187-0.75, -121.2868-0.75, -120.8838-0.7, -120.6399-0.75, -120.6157-0.6)

long.site.labels.full <- c(-124.0593-1.7, -124.0848-1.75, -124.1148-1.75, -124.4015-1.5, -124.5647-1.5, -124.2529-2.25, 
                      -124.0809-1.75, -123.7895-1.5, -123.8036-1.75, -123.2551-1.6, -123.0740-2.75, -122.3976-1.75, 
                      -121.9537-1.5, -121.9290-1.9, -121.3187-2.75, -121.2868-1.75, -120.8838-1.2, -120.6399-1.2, -120.6157-1.2)

lat.site.labels <- c(44.83777+0.04, 44.50540, 44.24999-0.04, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461-0.06, 39.28090, 38.51198+0.04, 38.31900-0.05, 
               37.18506, 36.51939+0.1, 36.44750-0.09, 35.72893+0.12, 35.66549-0.07, 35.28994-0.02, 34.88117, 34.73024-0.12)


# ================================================================================== #

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# ================================================================================== #

# Create site map (site codes)
pdf("figures/Site_map.pdf", width = 8, height = 8)
g <- ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 5, 
             shape = 21, fill = mycolors) + coord_fixed(1.3) +
  geom_text(data=sites, aes(long.site.labels, lat.site.labels, label=site))+ 
  xlim(c(-134, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic() 
g + annotate(geom = 'text', size = 5,
         x = -127, y = 34.2, label = 'Pacific Ocean', fontface = 'italic')
dev.off()

# Create site map (full site names)
pdf("figures/Site_map_full_names.pdf", width = 8, height = 8)
g <- ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
             shape = 21, fill = mycolors) + coord_fixed(1.3) +
  geom_text(data=sites, aes(long.site.labels.full, lat.site.labels, label=site.full))+ 
  xlim(c(-134, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic()
g + annotate(geom = 'text', size = 5,
         x = -127, y = 34.2, label = 'Pacific Ocean', fontface = 'italic')
dev.off()
