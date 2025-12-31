## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
# 
# Plot temperature data and trends
#
#
library(dplyr);library(ggplot2);library(gridExtra);library(here)

theme_set(theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major.y = element_blank()))


#GOA SST from NCAR/NCEP reanalysis
SST_GOA<-read.csv(here("Data","Environmental","NCEP_reanalysisGOAmonthlySST_thru2022.csv"))
colnames(SST_GOA)[1]<-"YEAR"

#EBS temps from Bering10K ROMS hindcasts. SST is ice-influenced - use depth-integrated.
#Depth-integrated (DepthAvg) temperatures, SEBS, monthly. 1970-2022
#Called SST here for coding simplicity only
#
SST_EBS<-read.csv(here("Data","Environmental","EBStemps.csv"))

## Identify years with larval data included in this study

goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")

goayrs<-unique(goa$YEAR)
ebsyrs<-unique(ebs$YEAR)

SST_GOA <- SST_GOA %>%
  mutate(inData = YEAR %in% goayrs) %>%
  rowwise() %>%
  mutate(JFMAM=mean(c(JAN,FEB,MAR,APR,MAY))) %>%
  mutate(ANNUAL=mean(c(JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC))) %>%
  filter(YEAR >= 1972)

SST_EBS <- SST_EBS %>%
  mutate(inData = YEAR %in% ebsyrs) %>%
  rowwise() %>%
  mutate(JFMAM=mean(c(JAN,FEB,MAR,APR,MAY))) %>%
  mutate(ANNUAL=mean(c(JAN,FEB,MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC))) %>%
  filter(YEAR >= 1972)


## Create manuscript supplementary figure.

EBSplot<-ggplot(SST_EBS,aes(YEAR,JFMAM)) +
  geom_abline(intercept=mean(SST_EBS$JFMAM,na.rm=T),slope=0,color="gray60",linetype=2) + 
  geom_path(color="gray",linewidth=0.6) +
  geom_point(data=SST_EBS,aes(shape=inData),color="steelblue3") +
  geom_smooth(method='lm',linetype=2,color=1,se=FALSE) +
  geom_smooth(method='lm',data=SST_EBS[SST_EBS$inData==T,],aes(YEAR,JFMAM),color="steelblue2",fill="steelblue2",alpha=0.2) +
  geom_smooth(method='lm',data=SST_EBS,aes(YEAR,JFMAM),color="gray60",fill="gray60",alpha=0.2) +
  theme(legend.position="none") +
  scale_shape_manual(values = c(NA,19)) +
  ylab("JFMAM Depth-Avg Temp. (°C)") +
  ggtitle(label="Bering Sea") +
  xlab("Year")

EBSplot


GOAplot<-ggplot(SST_GOA,aes(YEAR,JFMAM)) +
  geom_abline(intercept=mean(SST_GOA$JFMAM,na.rm=T),slope=0,color="gray60",linetype=2) +
  geom_path(color="gray",linewidth=0.6) +
  geom_point(data=SST_GOA,aes(shape=inData),color="steelblue3") +
  geom_smooth(method='lm',linetype=2,color=1,se=FALSE) +
  geom_smooth(method='lm',data=SST_GOA[SST_GOA$inData==T,],aes(YEAR,JFMAM),color="steelblue2",fill="steelblue2",alpha=0.2) +
  geom_smooth(method='lm',data=SST_GOA,aes(YEAR,JFMAM),color="gray60",fill="gray80",alpha=0.2) +
  theme(legend.position="none") +
  scale_shape_manual(values = c(NA,19)) +
  ylab("JFMAM Surface Temp. (°C)") +
  ggtitle(label="Gulf of Alaska") +
  xlab("")

GOAplot


plot2<-grid.arrange(GOAplot,EBSplot,ncol=1)

ggsave(here("Figures",paste0("TemperatureTimeseries_",Sys.Date(),".png")),plot2,width=6,height=6)


### Regress temperature by year to test for long-term linear trends.
## For all years 1970-2022, and for years with data included in this study.

EBS_lm1<-lm(JFMAM ~ YEAR, data=SST_EBS[SST_EBS$inData ==T,])
EBS_lm2<-lm(JFMAM ~ YEAR, data=SST_EBS)

GOA_lm1<-lm(JFMAM ~ YEAR, data=SST_GOA[SST_GOA$inData ==T,])
GOA_lm2<-lm(JFMAM ~ YEAR, data=SST_GOA)


