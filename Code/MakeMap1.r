## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Map of sampling regions (GOA and EBS) showing sampling intensity
## 
## 

library(ggplot2)
source("Code/MyFunctions.r")
library(akmarineareas2)
library(marmap)
#theme_update(panel.grid.major = element_line(color="gray"), panel.grid.minor = element_line(color="gray"))

goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")
goa_ebs<-bind_rows(goa,ebs)

## This geographic filtering now in SubsetLengthData... 
#goa_ebs<-goa_ebs %>% 
#  filter(LON > -175) %>%
#  filter(!(LON < -170 & LAT< 54)) %>%
#  filter(LAT<59.8)

goa_ebs_proj<- goa_ebs %>% 
  st_as_sf(coords = c("LON","LAT"),crs="EPSG:4326") %>%
  sf::st_transform(crs = 3338) #code for EA Albers

goa_ebs_xy<-as.data.frame(st_coordinates(goa_ebs_proj)) 

goa_ebs_xy<-bind_cols(goa_ebs,goa_ebs_xy)


#map_data <- rnaturalearth::ne_countries(
#  returnclass = "sf", scale = "large", ##NEED scale = "large" in hires data
#  country = "united states of america")

#ak <- subset(map_data, name == "Alaska")

# Crop the polygon for plotting and efficiency:
#ak_coast <- suppressWarnings(suppressMessages(
#  st_crop(map_data,
#          c(xmin = -179.5, ymin = 51.5, xmax = -145, ymax = 63))))

#CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")  ##Equal area Albers"

ak_coast_albers_proj <- sf::st_transform(ak_coast, crs = 3338) #code for EA Albers
st_bbox(ak_coast_albers_proj)
#xmin       ymin       xmax       ymax 
#-1684245.1   423767.1   490718.3  1646476.1 

# pull bathymetry contours coordinates
bat <- getNOAA.bathy(-179, -140, 50, 70, res = 4, keep = TRUE,path="Data/")
bat_xyz <- as.xyz(bat)
# this projects the contours
bat_proj <- as.raster(bat) %>% 
  raster::projectRaster(crs = 3338) %>% 
  as.bathy() %>% 
  as.xyz() %>%
  as_tibble() %>%
  drop_na(V3)


dat_to_plot<-goa_ebs_xy %>% group_by(HAUL_ID) %>% summarize(nspec=n(),X=mean(X),Y=mean(Y))
scale_name<-"Number \nof hauls"
plot_type<-"hauls"

dat_to_plot<-goa_ebs_xy
scale_name<-"Number of \nspecimens"
plot_type<-"specimens"

ggplot() +
  theme_light()  +  
  geom_hex(data=dat_to_plot,aes(x=X,y=Y),bins=30) +
  scale_fill_viridis_c(trans = "log10",name=scale_name,alpha=0.8) +
  geom_contour(data = bat_proj, # light gray contours for 1000m
               aes(x = V1, y = V2, z = V3),
               breaks = c( -200), color = "grey40", linewidth = 0.2,
               linetype = 2, show.legend=TRUE) +
  geom_sf(data=ak_coast_albers_proj,color=NA,fill="gray20")  +
  scale_x_continuous(breaks = seq(-180, -145, 5),limits=c(-1280000,350000)) +
  #  xlim(-1450000,350000) +
  ylim(350000,1300000) +
  ylab("Latitude") +
  xlab("Longitude") +
  coord_sf(expand = FALSE) +
  annotate("text",x=-850000, y=950000,label="Eastern Bering Sea",size=6) +
  annotate("text",x=100000, y=550000,label="Gulf of Alaska",size=6)




ggsave(here("Figures","Maps",paste0("Map1_sampling_",plot_type,"_",Sys.Date(),".png")),width=10,height=6)





