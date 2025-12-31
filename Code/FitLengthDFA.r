## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Estimate single-trend DFA on mean estimated length-at-date, for 28 spp (GOA)
## and 6 spp (EBS)
## 
## 

source("Code/MyFunctions.r")
library(ggpubr)
today<-Sys.Date()

# #Load year effects from Year Factor models

date<-"2025-11-17"
myregion<-"goa"
myfoldername=paste0(myregion,"_YearFModelFits_k5_noSP_gauss")
mod_YearF_lists_goa<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
#date<-"2025-09-26"
#myfoldername=paste0(myregion,"_YearFModelFits_k5_noSP_genGam")
#mod_YearF_lists_goa<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))

date<-"2025-11-17"
myregion<-"ebs"
myfoldername=paste0(myregion,"_YearFModelFits_k5_noSP_gauss")
mod_YearF_lists_ebs<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))

## Temperature data for comparison with trends
# GOA SST
SST_GOA<-read.csv(here("Data","Environmental","NCEP_reanalysisGOAmonthlySST_thru2022.csv"))
colnames(SST_GOA)[1]<-"YEAR"
SST_GOA <- SST_GOA %>%
  #  mutate(inData = YEAR %in% goayrs) %>%
  rowwise() %>%
  mutate(JFMAM=mean(c(JAN,FEB,MAR,APR,MAY))) %>%
  filter(YEAR >= 1970)
# EBS depth-averaged temp from Bering10K model
SST_EBS<-read.csv(here("Data","Environmental","EBStemps.csv"))
SST_EBS <- SST_EBS %>%
  #  mutate(inData = YEAR %in% goayrs) %>%
  rowwise() %>%
  mutate(JFMAM=mean(c(JAN,FEB,MAR,APR,MAY))) %>%
  filter(YEAR >= 1970)

## Species names and common names table (all)

spp_all<-read.csv(here("Data","SpeciesNamesCodes.csv"))


## Set DFA params 
chains <- 4
iter <- 4000

#################
## GOA DFA
#################

# Merge Year Effects into one dataframe, prep data for DFA
inds<-which(mod_YearF_lists_goa %>% map(1)!="error")
allfits_YearF_g <- mod_YearF_lists_goa[inds] %>% map(6) %>% dplyr::bind_rows()

yrs_g<-as.numeric(as.character(unique(allfits_YearF_g$YEARF)))
spp_g<-unique(allfits_YearF_g$SPECIES)

yrvec_g<-min(yrs_g):max(yrs_g)
yrgrid_g<-expand.grid(yrvec_g,spp_g) %>% 
  rename(YEARF=Var1,SPECIES=Var2) %>%
  mutate(YEARF=as.factor(YEARF))
  
yrcomp_wide_g<-allfits_YearF_g %>%
  right_join(yrgrid_g) %>%
  select(YEARF,visregFit,SPECIES) %>%
  mutate(YEARF=as.numeric(as.character(YEARF))) %>%
  arrange(SPECIES,YEARF) %>%
  pivot_wider(names_from=SPECIES,values_from=visregFit) %>% 
  select(-YEARF) %>%
  t() 

sppnames_g<-rownames(yrcomp_wide_g)
common_names_g<-data.frame("SPECIES_NAME"=sppnames_g) %>% 
  left_join(spp_all) %>%
  select(COMMON_NAME)

yrcomp_long_g<-allfits_YearF_g %>%
  mutate(Year=as.numeric(as.character(YEARF))) %>%
  mutate(se=(visregUpr-visregLwr)/3.92) %>% ## 95%CI is +/- 1.96 SE
  mutate(weights = (1 / se)^2) %>%
  mutate(YearNum=Year-min(Year)+1) %>%
  select(YearNum,Year,visregFit,SPECIES, se,weights) %>%
  mutate(SPECIES=as.factor(SPECIES)) %>%
  rename(ts=SPECIES,obs=visregFit,se=se,time=YearNum) %>%
  as.data.frame()

# Quick check on weights  
allwts<- yrcomp_long_g %>%
  group_by(ts) %>%
  summarise(semean=mean(se))

## Fit DFA
## One Trend Model
## Note that labels for loadings plots are easy to mismatch.
## Use the wide format, and select spp names directly from that dataframe.
## If using weights, need long format 


set.seed(103)  #103 gives positive loading for GOA unweighted. 

f1_wide_g <- fit_dfa(
  y = yrcomp_wide_g, num_trends = 1, scale="zscore",
  iter = iter, chains = chains, thin = 1
) 

saveRDS(f1_wide_g,file=here("Results","DFAResults",paste0(today(),"_GOA_DFA_",iter,"iter.rds")))

set.seed(105)
f1_long_wt_g <- fit_dfa(
  y = yrcomp_long_g, num_trends = 1, scale="zscore",data_shape = "long",
  iter = iter, chains = chains, thin = 1, inv_var_weights = "weights")

saveRDS(f1_long_wt_g,file=here("Results","DFAResults",paste0(today(),"_GOA_DFA_wt",iter,"iter.rds")))


## To save time, once DFA has been fit,
## read in previous DFA, selected weighted or unweighted, then produce plots.

f1_wide_g<-readRDS(file=here("Results","DFAResults",paste0("2025-11-17_GOA_DFA_6000iter.rds")))
is_converged(f1_wide_g, threshold = 1.05)

f1_long_wt_g<-readRDS(file=here("Results","DFAResults",paste0("2025-11-18_GOA_DFA_wt4000iter.rds")))
is_converged(f1_long_wt_g, threshold = 1.05)


rot_g <- rotate_trends(f1_wide_g) #unweighted DFA
#rot_g <- rotate_trends(f1_long_wt_g) #weighted DFA

#test plot
plot_loadings(rot_g,names=common_names_g$COMMON_NAME)

####### extract trends for comparison with temperature
tr<-t(rot_g$trends_mean)
trends_temp<-data.frame("spptrend"=tr,"YEAR"=yrvec_g)
colnames(trends_temp)[1]<-"spptrend"

trends_temp_GOA<-trends_temp %>% 
  left_join(SST_GOA) %>%
  filter(YEAR %in% yrs_g)
  
cor(trends_temp_GOA)

ggplot(trends_temp_GOA %>% select(spptrend,YEAR,JFMAM), aes(YEAR,spptrend)) +
  geom_path() +
  geom_path(aes(YEAR,JFMAM),col=2)

goa_temp<-ggplot(trends_temp_GOA, aes(JFMAM,spptrend)) +
  geom_point() +
  geom_text_repel(aes(label=YEAR),size=3,color="gray60") +
  ylab("DFA Trend") +
  xlab("JFMAM Temperature (°C)")

goa_temp
#ggsave(here("Figures",paste0("GOA_DFA_Trend_vTemp",today(),".png")))

uniqueyrs_g<-yrs_g
allyrs_g<-min(yrs_g):max(yrs_g)

##
#W Plot loadings and trends
goa_l<-plot_loadings_LR(rot_g,names=common_names_g$COMMON_NAME)
goa_l
goa_l2 <- goa_l +  theme(legend.position = "none")
goa_l2
#ggsave(here("Figures",paste0(myregion,"_DFA_Loadings_",today(),".png")),width=8,height=4)

goa_tr<-plot_trends_LR(rot_g,allyears=allyrs_g,datayears=uniqueyrs_g) +
  ggtitle("GOA") +
  theme(plot.title = element_text(hjust = 0.5))
goa_tr
#ggsave(here("Figures",paste0(myregion,"_DFA_Trend_",today(),".png")),width=8,height=6)


#################
##  EBS
#################

# Merge Year Effects into one dataframe, prep data for DFA 
inds<-which(mod_YearF_lists_ebs %>% map(1)!="error")
allfits_YearF_e <- mod_YearF_lists_ebs[inds] %>% map(6) %>% dplyr::bind_rows()

yrs_e<-as.numeric(as.character(unique(allfits_YearF_e$YEARF)))
spp_e<-unique(allfits_YearF_e$SPECIES)

yrvec_e<-min(yrs_e):max(yrs_e)
yrgrid_e<-expand.grid(yrvec_e,spp_e) %>% 
  rename(YEARF=Var1,SPECIES=Var2) %>%
  mutate(YEARF=as.factor(YEARF))

yrcomp_wide_e<-allfits_YearF_e %>%
  right_join(yrgrid_e) %>%
  select(YEARF,visregFit,SPECIES) %>%
  mutate(YEARF=as.numeric(as.character(YEARF))) %>%
  arrange(SPECIES,YEARF) %>%
  pivot_wider(names_from=SPECIES,values_from=visregFit) %>% 
  select(-YEARF) %>%
  t() 

sppnames_e<-rownames(yrcomp_wide_e)
common_names_e<-data.frame("SPECIES_NAME"=sppnames_e) %>% 
  left_join(spp_all) %>%
  select(COMMON_NAME)

yrcomp_long_e<-allfits_YearF_e %>%
  mutate(Year=as.numeric(as.character(YEARF))) %>%
  mutate(se=(visregUpr-visregLwr)/3.92) %>% ## 95%CI is +/- 1.96 SE
  mutate(weights = (1 / se)^2) %>%
  mutate(YearNum=Year-min(Year)+1) %>%
  select(YearNum,Year,visregFit,SPECIES, se,weights) %>%
  mutate(SPECIES=as.factor(SPECIES)) %>%
  rename(ts=SPECIES,obs=visregFit,se=se,time=YearNum) %>%
  as.data.frame()

## Fit DFA
## One Trend Model
## Unweighted and weighted

set.seed(102)  #102 gives positive loading for EBS, neg for GOA. 

f1_wide_e <- fit_dfa(
  y = yrcomp_wide_e, num_trends = 1, scale="zscore",
  iter = iter, chains = chains, thin = 1
) 

saveRDS(f1_wide_e,file=here("Results","DFAResults",paste0(today(),"_EBS_DFA_",iter,"iter.rds")))

set.seed(100)  #102 gives positive loading for EBS, neg for GOA. 

f1_long_wt_e <- fit_dfa(
  y = yrcomp_long_e, num_trends = 1, scale="zscore",data_shape = "long",
  iter = iter, chains = chains, thin = 1, inv_var_weights = "weights")

saveRDS(f1_long_wt_e,file=here("Results","DFAResults",paste0(today(),"_EBS_DFA_wt",iter,"iter.rds")))

## Read in DFA fits, selected weighted or unweighted, then produce plots.

f1_wide_e<-readRDS(file=here("Results","DFAResults",paste0("2025-11-17_EBS_DFA_4000iter.rds")))
is_converged(f1_wide_e, threshold = 1.05)

f1_long_wt_e<-readRDS(file=here("Results","DFAResults",paste0("2025-11-18_EBS_DFA_wt4000iter.rds")))
is_converged(f1_long_wt_e, threshold = 1.05)

rot_e <- rotate_trends(f1_wide_e) #unweighted
#rot_e <- rotate_trends(f1_long_wt_e)  #weighted

#test plot
plot_loadings(rot_e,names=common_names_e$COMMON_NAME)

####### extract trends for comparison with temperature
tr<-t(rot_e$trends_mean)
trends_temp<-data.frame("spptrend"=tr,"YEAR"=yrvec_e)
colnames(trends_temp)[1]<-"spptrend"

### FOR EBS TREND ###

trends_temp_EBS<-trends_temp %>% 
  left_join(SST_EBS) %>%
  arrange(YEAR) %>%
  filter(YEAR %in% yrs_e)

cor(trends_temp_EBS)

ggplot(trends_temp_EBS %>% select(spptrend,YEAR,JFMAM), aes(YEAR,spptrend)) +
  geom_path() +
  geom_path(aes(YEAR,JFMAM),col=2)

ebs_temp<-ggplot(trends_temp_EBS, aes(JFMAM,spptrend)) +
  geom_point() +
  geom_text_repel(aes(label=YEAR),size=3,color="gray60") +
  ylab("DFA Trend") +
  xlab("JFMAM Temperature (°C)")

ebs_temp
#ggsave(here("Figures",paste0("EBS_DFA_Trend_vTemp",today(),".png")),width=6,height=5)

ebs_l<-plot_loadings_LR(rot_e,names=common_names_e$COMMON_NAME)
ebs_l2<-ebs_l +  theme(legend.position = "top",legend.direction = "horizontal")
ebs_l2
#ggsave(here("Figures",paste0(myregion,"_DFA_Loadings_",today(),".png")),width=8,height=4)

uniqueyrs_e<-yrs_e
allyrs_e<-min(yrs_e):max(yrs_e) # but use GOA years if plotting together
ebs_tr<-plot_trends_LR(rot_e,allyears=allyrs_e,datayears=uniqueyrs_e) +
  ggtitle("EBS") +
  theme(plot.title = element_text(hjust = 0.5))
ebs_tr
ebs_tr<-ebs_tr + xlim(min(yrs_g),max(yrs_g))
#ggsave(here("Figures",paste0(myregion,"_DFA_Trend_",today(),".png")),width=8,height=6)

##############
##############
## Plot timeseries and loadings, with temp relationship, for GOA, EBS together 
##############
##############

ggarrange(goa_tr,goa_l2,goa_temp,ebs_tr,ebs_l2,ebs_temp,ncol=3,nrow=2)

ggsave(here("Figures","DFA",paste0("DFA_Compiled_",iter,"iter_",today(),".png")),width=14,height=8)
#ggsave(here("Figures","DFA",paste0("DFA_Compiled_",iter,"iter_weighted_",today(),".png")),width=14,height=8)

ggarrange(goa_tr,ebs_tr,goa_temp,ebs_temp,goa_l2,ebs_l2,ncol=2,nrow=3,
          heights=c(1.5,1.8,2),labels = c("(a)","(b)","(c)","(d)","(e)","(f)"))

ggsave(here("Figures","DFA",paste0("DFA_Compiled_vert_",iter,"iter_",today(),".png")),width=10,height=12)
#ggsave(here("Figures","DFA",paste0("DFA_Compiled_vert_",iter,"iter_weighted_",today(),".png")),width=10,height=12)
ggsave(here("Figures","DFA",paste0("DFA_Compiled_vert_",iter,"iter_",today(),".pdf")),width=10,height=12)


