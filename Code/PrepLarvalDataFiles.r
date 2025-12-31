## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Prep larval data files by adding mean monthly temperatures.
## 
## Larval data files were previously queried from NOAA AFSC EcoFOCI database,
## and trimmed to include only observations within a particular seasonal time window
## and geographic area. 
## 
## For the EBS: 
# EBS <- IchL %>%
# filter(GEOGRAPHIC_AREA=="BS") %>%
#   filter(LAT < 60 & LAT > 50) %>%  
#   filter(LON > -175) %>%  #to remove stations not adjacent to others
#   filter(!(LON < -170 & LAT< 54.05)) %>%
#   filter(!(LON < -173 & LAT< 55.2))
#   
## For the GOA:
# GOA <- IchL %>%
#   filter(GEOGRAPHIC_AREA=="GOA") %>%
#   filter(LON < -150) %>% 
#   filter(LON > -165.2) %>%
#   filter(LAT < 59.15)   
#   
## Taxa were then filtered out if they are not identified to species, or if they occur in fewer
## than 10 years.
## 
source("Code/MyFunctions.r")

# Larval data
GOA<-read.csv(here("Data","Larval","GOA_spring_d117to160_32spp_Jul2024_larval.csv"))
EBS<-read.csv(here("Data","Larval","EBS_spring_d98to162_11spp_Apr2024_larval.csv"))

# Read in temperature data to merge with larval data.
SST_GOA<-read.csv(here("Data","Environmental","NCEP_reanalysisGOAmonthlySST_thru2022.csv"))
colnames(SST_GOA)[1]<-"YEAR"

Temps_EBS<-read.csv(here("Data","Environmental","EBStemps.csv")) 

# Merge temperature files with larval files
EBS2<-merge(EBS,Temps_EBS,all.x=TRUE)
GOA2<-merge(GOA,SST_GOA,all.x=TRUE)

# Save in "processed" folder
# 

write.csv(GOA2,here("Data","Processed","GOA_spring_d117to160_32spp_Jul2024.csv"),row.names=F)
write.csv(EBS2,here("Data","Processed","EBS_spring_d98to162_11spp_Apr2024.csv"),row.names=F)

# These are the data files used in further code.