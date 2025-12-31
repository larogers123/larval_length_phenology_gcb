## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
##
## Plot summary maps of survey stations by year, with date of sampling
## Maps included in supplementary materials.
## 


# Read in data prep and mapping functions, packages, graphics defaults
source("Code/MyFunctions.r")
library(marmap)
library(ggforce)
today<-Sys.Date()


######################
### PREP DATA
######################

goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")
goa_ebs<-bind_rows(goa,ebs)

goa_ebs_proj<- goa_ebs %>% 
  st_as_sf(coords = c("LON","LAT"),crs="EPSG:4326") %>%
  sf::st_transform(crs = 3338) #code for EA Albers

goa_ebs_xy<-as.data.frame(st_coordinates(goa_ebs_proj)) 

goa_ebs_xy<-bind_cols(goa_ebs,goa_ebs_xy)


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


dat_to_plot<-goa_ebs_xy %>% 
  group_by(HAUL_ID) %>% 
  summarize(YDAY=mean(YDAY), YEAR=mean(YEAR), X=mean(X),Y=mean(Y))

scale_name<-"Day of Year"
plot_type<-"hauls by doy"

## Plot across two pages
for(pages in 1:2){

ggplot() +
  theme_light()  +  
  geom_point(data=dat_to_plot,aes(x=X,y=Y,color=YDAY)) +
#  scale_fill_viridis_c(name=scale_name) +
  labs(color=scale_name) +
  geom_contour(data = bat_proj, # light gray contours for 1000m
               aes(x = V1, y = V2, z = V3),
               breaks = c( -200), color = "grey40", linewidth = 0.2,
               linetype = 2, show.legend=TRUE) +
  geom_sf(data=ak_coast_albers_proj,color=NA,fill="gray20")  +
  scale_x_continuous(breaks = seq(-178, -146, 8),limits=c(-1280000,350000)) +
  #  xlim(-1450000,350000) +
  ylim(350000,1300000) +
  ylab("Latitude") +
  xlab("Longitude") +
  coord_sf(expand = FALSE) +
  facet_wrap_paginate(~YEAR,nrow=6,ncol=4,page=pages)

ggsave(here("Figures","Maps",paste0("MapS1_",plot_type,pages,".png")),width=8,height=9)
}



## Plot one large landscape page (43 panels) (5x9, 6x8)

ggplot() +
  theme_light()  +  
  geom_point(data=dat_to_plot,aes(x=X,y=Y,color=YDAY),size=0.8) +
  #  scale_fill_viridis_c(name=scale_name) +
  labs(color=scale_name) +
  geom_contour(data = bat_proj, # light gray contours for 1000m
               aes(x = V1, y = V2, z = V3),
               breaks = c( -200), color = "grey40", linewidth = 0.2,
               linetype = 2, show.legend=TRUE) +
  geom_sf(data=ak_coast_albers_proj,color=NA,fill="gray20")  +
  scale_x_continuous(breaks = seq(-178, -146, 8),limits=c(-1280000,350000)) +
  #  xlim(-1450000,350000) +
  ylim(350000,1300000) +
  ylab("Latitude") +
  xlab("Longitude") +
  coord_sf(expand = FALSE) +
  facet_wrap(~YEAR,nrow=7,ncol=7) #+
#  theme(legend.position = "bottom") #,
#    legend.position.inside = c(0.7,0.0))

ggsave(here("Figures","Maps",paste0("MapS1_",plot_type,"_landscape.png")),width=11,height=8.5)

