## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Create plots of conditional effects for all species, 
## for temperature and year factor models.
## 

# Read in data prep and mapping functions, packages, graphics defaults
source("Code/MyFunctions.r")
library(ggpubr)
today<-Sys.Date()

goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")

goa_ebs<-bind_rows(goa,ebs)

spnames<-read.csv("Data/SpeciesNamesCodes.csv")

######################
######################
### LOAD MODEL RESULTS
######################
######################

myregion<-"goa"

date<-"2025-09-30"
#myfoldername="BaseModelFits_k5_noSP_gauss"
#mod_base_lists_goa<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myfoldername="YearFModelFits_k5_noSP_gauss"
myYearFfoldername<-myfoldername

mod_YearF_lists_goa<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myfoldername="TempModelFits_k5_noSP_gauss_JFMAM"
myTempfoldername<-myfoldername

mod_temp_lists_goa<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myfoldername="LinearYearModelFits_k5_noSP_gauss"
mod_year_lists_goa<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myregion<-"ebs"

date<-"2025-12-03"

#myfoldername="BaseModelFits_k5_noSP_gauss"
#mod_base_lists_ebs<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myfoldername="YearFModelFits_k5_noSP_gauss"
mod_YearF_lists_ebs<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myfoldername="TempModelFits_k5_noSP_gauss_JFMAM"
mod_temp_lists_ebs<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

myfoldername="LinearYearModelFits_k5_noSP_gauss"
mod_year_lists_ebs<-readRDS(here("Results","ModelResults",paste0(myregion,"_",myfoldername,"_",date),"AllSppLists.rds"))

################################
################################
## PLOTTING YEAR FACTOR MODELS
################################
################################
#need to exclude lists with errors
#then merge into one dataframe for ggplot
#inds<-which(mod_YearF_lists %>% map(1)!="error")
#allcoefs_YearF <- mod_YearF_lists[inds] %>% map(1) %>% dplyr::bind_rows()
#allpartials_YearF <- mod_YearF_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_YearF_goa <- mod_YearF_lists_goa %>% map(6) %>% 
  dplyr::bind_rows() %>% mutate(region="goa")
allfits_YearF_ebs <- mod_YearF_lists_ebs %>% map(6) %>% 
  dplyr::bind_rows() %>% mutate(region="ebs")

allfits_YearF<-bind_rows(allfits_YearF_goa,allfits_YearF_ebs) %>%
  left_join(goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],
           by=join_by("SPECIES"=="SPECIES_NAME"), multiple="any")

#create species/common names d.f., and determine y plot location for printing
spdf<-data.frame("SPECIES"=unique(allfits_YearF$SPECIES))
spdf<-left_join(spdf,goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],by=join_by("SPECIES"=="SPECIES_NAME"),multiple="any")
altmax<-allfits_YearF %>% filter(YDAYSC==0) %>% group_by(SPECIES,region) %>% summarize(maxest= max(visregUpr))
spdf<-full_join(spdf,altmax) %>% group_by(SPECIES,region) %>% mutate(maxy=maxest+1)

#convert year to numeric for plotting
allfits_YearF <-allfits_YearF %>%
  mutate(YEAR=as.numeric(levels(YEARF))[YEARF]) 

minY<-min(allfits_YearF$YEAR)

p1<-ggplot(allfits_YearF %>% filter(region=="goa"),aes(YEAR,visregFit)) +
  geom_point() +
  geom_errorbar(aes(x=YEAR,ymin=visregLwr,ymax=visregUpr)) +
  facet_wrap(~COMMON_NAME,scales="free_y") +
  geom_text(data=spdf%>% filter(region=="goa"),aes(label=SPECIES,y=maxy),x=minY+1,hjust=0,size=3,color="gray30",fontface = "italic") +
#  theme(strip.text = element_text(face = "italic")) +
  ylab("Larval Length (mm, May 25)") +
  xlab("Year") +
  xlim(1972,2022) + 
  theme(plot.margin = margin(1,0.1,0.2,1, 'cm'))


p2<-ggplot(allfits_YearF %>% filter(region=="ebs"),aes(YEAR,visregFit)) +
  geom_point() +
  geom_errorbar(aes(x=YEAR,ymin=visregLwr,ymax=visregUpr)) +
  facet_wrap(~COMMON_NAME,scales="free_y",nrow=1) +
  geom_text(data=spdf%>% filter(region=="ebs"),aes(label=SPECIES,y=maxy),x=minY+1,hjust=0,size=3,color="gray30",fontface = "italic") +
  #  theme(strip.text = element_text(face = "italic")) +
  ylab("Larval Length (mm, May 25)") +
  xlab("Year") +
  xlim(1972,2022) + 
  theme(plot.margin = margin(1,0.1,0.1,1, 'cm'))

ggarrange(p1,p2,ncol=1,heights = c(4,1),align="v",  labels = c("(a) GOA", "(b) EBS"))

ggsave(file=here("Figures",paste0(myYearFfoldername,"_combined_",today,".png")),width=14,height=14)
ggsave(file=here("Figures",paste0(myYearFfoldername,"_combined_",today,".pdf")),width=14,height=14)

 
# # Single species plot
# ggplot(allfits_YearF %>% filter(SPECIES=="Gadus chalcogrammus"),aes(YEAR,visregFit)) +
#   geom_point() +
#   geom_errorbar(aes(x=YEAR,ymin=visregLwr,ymax=visregUpr)) +
#   facet_wrap(~SPECIES,scales="free_y") +
#   geom_text(data=spdf %>% filter(SPECIES=="Gadus chalcogrammus"),aes(label=COMMON_NAME,y=maxy),x=1980,hjust=0,size=3,color="gray30") +
#   theme(strip.text = element_text(face = "italic")) +
#   ylab("Larval Length (mm)") +
#   xlab("Year") +
#   theme(panel.grid.major = element_line(),panel.grid.minor.x = element_line())
# 
# ggsave(file=here("Figures",paste0(myfoldername,"_pollock.png")),width=6,height=4)


################################
################################
## PLOTTING TEMPERATURE MODELS
################################
################################

# color palette
col1<-'#F87505'
col2<-'#048F96'

## need to exclude species with errors
#inds<-which(mod_temp_lists %>% map(3)!="error")
allcoefs_temp_goa <- mod_temp_lists_goa %>% map(1) %>%
  dplyr::bind_rows() %>% mutate(region="goa")
allcoefs_temp_ebs <- mod_temp_lists_ebs %>% map(1) %>%
  dplyr::bind_rows() %>% mutate(region="ebs")#allpreds_temp <- mod_temp_lists[inds] %>% map(3) %>% dplyr::bind_rows()

allfits_temp_goa <- mod_temp_lists_goa %>% map(6) %>%
  dplyr::bind_rows() %>% mutate(region="goa")
allfits_temp_ebs <- mod_temp_lists_ebs %>% map(6) %>%
  dplyr::bind_rows() %>% mutate(region="ebs")
allfits_temp<-bind_rows(allfits_temp_goa,allfits_temp_ebs) %>%
  left_join(goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],
            by=join_by("SPECIES"=="SPECIES_NAME"), multiple="any")

allpartials_temp_goa <- mod_temp_lists_goa %>% map(5) %>%
  dplyr::bind_rows() %>% mutate(region="goa")
allpartials_temp_ebs <- mod_temp_lists_ebs %>% map(5) %>%
  dplyr::bind_rows() %>% mutate(region="ebs")
allpartials_temp<-bind_rows(allpartials_temp_goa,allpartials_temp_ebs) %>%
  left_join(goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],
            by=join_by("SPECIES"=="SPECIES_NAME"), multiple="any")

#allpreds_temp <- mod_temp_lists[inds] %>% map(3) %>% dplyr::bind_rows()
#allpartials_temp <- mod_temp_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
#allfits_temp <- mod_temp_lists[inds] %>% map(6) %>% dplyr::bind_rows()
#allmesh_temp <- mod_temp_lists[inds] %>% map(8) %>% dplyr::bind_rows()

tempvar="JFMAM"
TEMPeffects<-bind_rows(allcoefs_temp_goa,allcoefs_temp_ebs) %>%
  left_join(goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],
            by=join_by("SPECIES"=="SPECIES_NAME"), multiple="any") %>%
  filter(term==tempvar) %>% arrange(SPECIES) %>%
  mutate(ispos = estimate>0, 
         signif = sign(conf.low)==sign(conf.high)) %>%
  mutate(slope = case_when(ispos == TRUE ~ 'positive',
                           ispos == FALSE ~ 'negative')) %>%
  mutate(significance = case_when(signif ==TRUE ~ 'p < 0.05',
                                  signif ==FALSE ~ 'NS'))

fitted_wPosNeg<- TEMPeffects %>%  select(SPECIES,slope,significance,region) %>%
  full_join(allfits_temp)

#create species/common names d.f. and 
spdf<-data.frame("SPECIES"=unique(TEMPeffects$SPECIES))
spdf<-left_join(spdf,goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],by=join_by("SPECIES"=="SPECIES_NAME"),multiple="any")

#Find max y for each species for labels
yloc<-allpartials_temp %>% filter(!is.infinite(visregRes)) %>% group_by(SPECIES,region)  %>% summarize(yloc=max(visregRes)) %>% select(SPECIES,yloc,region)
spdf<-full_join(spdf,yloc)

altmax<-fitted_wPosNeg %>% filter(YDAYSC==0) %>% group_by(SPECIES,region) %>% summarize(maxest= max(visregUpr))
spdf<-full_join(spdf,altmax) %>% group_by(SPECIES,region) %>% mutate(maxy=max(yloc,maxest))

#Find x location for label (minimum temp on x scale)
xloc_goa<-allfits_temp %>% filter(region=="goa") %>% select(all_of(tempvar)) %>% min()
xloc_ebs<-allfits_temp %>% filter(region=="ebs") %>% select(all_of(tempvar)) %>% min()

myregion<-"goa"
p3 <- ggplot() +
  geom_point(data=allpartials_temp %>% filter(region==myregion),aes(x=.data[[tempvar]],y=visregRes),alpha=0.2,size=0.2,color="gray") +
  geom_ribbon(data=fitted_wPosNeg %>% filter(region==myregion) %>% filter(YDAYSC==0),
              aes(.data[[tempvar]], visregFit, ymin = visregLwr,
                  ymax = visregUpr, fill=slope, alpha=significance)) +
  geom_line(data=fitted_wPosNeg %>% filter(region==myregion) %>% filter(YDAYSC==0), aes(.data[[tempvar]], visregFit,linetype=significance)) +
  facet_wrap(~COMMON_NAME,scales="free_y") +
  geom_text(data=spdf %>% filter(region==myregion) ,aes(label=SPECIES,y=maxy),x=xloc_goa,hjust=0,size=3,color="gray30",fontface = "italic") +
#  theme(strip.text = element_text(face = "italic")) +
  xlab(paste0(tempvar," Temperature (°C)")) +
  ylab("Larval Length (mm, May 25)") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_linetype_manual(values=c('NS' = "dashed", "p < 0.05" = "solid")) +
  scale_fill_manual(values = c("negative" = col1, "positive" = col2)) + 
  theme(plot.margin = margin(1,0.1,0.2,1, 'cm'))


myregion<-"ebs"
p4 <- ggplot() +
  geom_point(data=allpartials_temp %>% filter(region==myregion),aes(x=.data[[tempvar]],y=visregRes),alpha=0.2,size=0.2,color="gray") +
  geom_ribbon(data=fitted_wPosNeg %>% filter(region==myregion) %>% filter(YDAYSC==0),
              aes(.data[[tempvar]], visregFit, ymin = visregLwr,
                  ymax = visregUpr, fill=slope, alpha=significance)) +
  geom_line(data=fitted_wPosNeg %>% filter(region==myregion) %>% filter(YDAYSC==0), aes(.data[[tempvar]], visregFit,linetype=significance)) +
  facet_wrap(~COMMON_NAME,scales="free_y", nrow=1) +
  geom_text(data=spdf %>% filter(region==myregion) ,aes(label=SPECIES,y=maxy),x=xloc_ebs,hjust=0,size=3,color="gray30",fontface = "italic") +
  #  theme(strip.text = element_text(face = "italic")) +
  xlab(paste0(tempvar," Temperature (°C)")) +
  ylab("Larval Length (mm, May 25)") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_linetype_manual(values=c('NS' = "dashed", "p < 0.05" = "solid")) +
  scale_fill_manual(values = c("negative" = col1, "positive" = col2)) + 
  theme(plot.margin = margin(1,0.1,0.1,1, 'cm'),
        legend.position = "none")

ggarrange(p3,p4,ncol=1,heights = c(4,1),align="v",  labels = c("(a) GOA", "(b) EBS"))

ggsave(file=here("Figures",paste0(myTempfoldername,"_combined_",today,".png")),width=14,height=14)
ggsave(file=here("Figures",paste0(myTempfoldername,"_combined_",today,".pdf")),width=14,height=14)









### Single Species temperature plot
# ### 
# ggplot() +
# geom_point(data=allpartials_temp %>% filter(SPECIES=="Gadus chalcogrammus"),aes(x=.data[[tempvar]],y=visregRes),alpha=0.2,size=0.2,color="gray") +
#   geom_ribbon(data=toplot %>% filter(YDAYSC==0),
#               aes(.data[[tempvar]], visregFit, ymin = visregLwr,
#                   ymax = visregUpr, fill=slope, alpha=significance)) +
#   geom_line(data=toplot %>% filter(YDAYSC==0), aes(.data[[tempvar]], visregFit,linetype=significance)) +
#   facet_wrap(~SPECIES,scales="free_y") +
#   geom_text(data=spdf,aes(label=COMMON_NAME,y=maxy),x=xloc,hjust=0,size=3,color="gray30") +
#   theme(strip.text = element_text(face = "italic")) +
#   xlab(paste0(tempvar," Temperature (C)")) +
#   ylab("Larval Length (mm, May 25)") +
#   scale_alpha_discrete(range=c(0.1, 0.6)) +
#   scale_linetype_manual(values=c('NS' = "dashed", "p < 0.05" = "solid")) +
#   scale_fill_manual(values = c("negative" = col1, "positive" = col2))
# # 
