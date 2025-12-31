## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Summarize species included in study -
## years sampled, number of specimens, date ranges, size ranges, etc
## Table 1 in manuscript 

# Read in data prep and mapping functions, packages, graphics defaults
source("Code/MyFunctions.r")
today<-Sys.Date()


######################
### PREP DATA
######################


goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")


myspp<-unique(regiondat$SPECIES_NAME)

#Filter out species that are consistently not converging in GOA or EBS.
#Also leave out Capelin as they are previous year's spawn

if(myregion=="ebs"){ 
  myspp<-myspp[!myspp %in% c("Bathylagus pacificus",#"Hippoglossoides elassodon",
                             "Leuroglossus schmidti","Podothecus acipenserinus",
                             "Poroclinus rothrocki","Stenobrachius leucopsarus")]
}else if(myregion=="goa"){ 
  myspp<-myspp[!myspp %in% c("Anoplopoma fimbria","Lumpenella longirostris","Ronquilus jordani","Mallotus villosus")]
}


######################
### Summary tables for EBS, GOA
######################

#NOTE - make sure to update subsetting below if changed in model prep code
makeSppTable_fun<-function(dat=goa){
  my_table<-dat %>%
    mutate(DATETIME = ymd_hm(GMT_DATE_TIME_TXT)) %>%
    filter(!is.na(LARVALCATCHPER10M2)) %>%
    group_by(SPECIES_NAME,YEAR) %>%
    filter(n()>=10) %>%
    ungroup() %>%
    group_by(SPECIES_NAME) %>%
    summarize(CommonName=first(COMMON_NAME),
              FirstYear=min(YEAR),LastYear=max(YEAR),
              nYears=length(unique(YEAR)),
              nSpecimens=n(),
              meanLength=round(mean(CORRECTED_LENGTH),digits=1),
              Length5pct=round(quantile(CORRECTED_LENGTH,0.05),1),
              Length95pct=round(quantile(CORRECTED_LENGTH,0.95),1)) #%>%
 #   mutate(meanLength_5to95=paste0(meanLength," (",Length5pct,"-",Length95pct,")"))
  return(my_table)
  }


goa_table<-makeSppTable_fun(goa) %>%
  filter(!SPECIES_NAME %in% 
           c("Anoplopoma fimbria","Lumpenella longirostris",
             "Ronquilus jordani","Mallotus villosus"))

ebs_table<-makeSppTable_fun(ebs) %>%
  filter(!SPECIES_NAME %in% c("Bathylagus pacificus",
                              "Leuroglossus schmidti","Podothecus acipenserinus",
                              "Poroclinus rothrocki","Stenobrachius leucopsarus"))

write.csv(goa_table,paste0("Results/SppSummaryTable_goa",today,".csv"),row.names=F)
write.csv(ebs_table,paste0("Results/SppSummaryTable_ebs",today,".csv"),row.names=F)
