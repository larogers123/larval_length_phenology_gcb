#Prep map data for plots
#
library(rnaturalearth)
library(sf)

map_data <- rnaturalearth::ne_countries(
  returnclass = "sf", scale = "large", ##NEED scale = "large" in hires data
  country = "united states of america")

ak <- subset(map_data, name == "Alaska")

# Crop the polygon for plotting and efficiency:
ak_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -179.5, ymin = 51.5, xmax = -145, ymax = 63))))

ak_coast_goa_proj <- sf::st_transform(ak_coast, crs = 32604)
ak_coast_ebs_proj <- sf::st_transform(ak_coast, crs = 32603)


#read in goa prediction grid (based on IERP grid - rotated, AEA original projection)
goagrid<-read.csv("Data/goagrid.csv") %>%
  add_utm_columns(ll_names = c("LON", "LAT"),utm_crs=32604) %>%
  mutate(across(c(X, Y), \(x) round(x, digits = 1)))

#read in goa projection grid (based on UTM, will plot using raster in UTM)
goagrid2<-read.csv("Data/goagrid_utm.csv") %>%
  mutate(X=X/1000,Y=Y/1000)