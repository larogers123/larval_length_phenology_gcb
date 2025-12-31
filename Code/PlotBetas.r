## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Plot estimated slopes by species from models with temperature or linear year covariates
## 

source("Code/MyFunctions.r")
library(gridExtra);library(ggpubr);library(ggh4x)


## Species names and common names table
spp_all<-read.csv(here("Data","SpeciesNamesCodes.csv")) %>%
  rename(SPECIES=SPECIES_NAME) %>% select(SPECIES,COMMON_NAME)

#Specify which model family from which to extract/plot estimates
famname<-"gauss"

########### GOA #############

#load files with saved model run summaries to extract coefficients

if(famname=="gauss"){

  date<-"2025-09-30"
  myfoldername="goa_TempModelFits_k5_noSP_gauss_JFMAM"
  mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
  date<-"2025-09-30"
  myfoldername="goa_LinearYearModelFits_k5_noSP_gauss"
  mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
}

if(famname=="genGam"){
  
  ## GenGamma models
  date<-"2024-09-20"
  myfoldername="goa_TempModelFits_k5_noSP_genGam_JFMAM"
  mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
  myfoldername="goa_LinearYearModelFits_k5_noSP_genGam_"
  mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
}


if(famname=="logN"){
  ## logNormal models
  #load files with saved model run summaries to extract coefficients
  date<-"2025-08-11"
  myfoldername="goa_TempModelFits_k5_noSP_logN_JFMAM"
  mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
  myfoldername="goa_LinearYearModelFits_k5_noSP_logN"
  mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
}

inds<-which(mod_year_lists %>% map(3)!="error")
allcoefs_year <- mod_year_lists[inds] %>% map(1) %>% dplyr::bind_rows() %>%
  left_join(spp_all)

YEAReffects<-allcoefs_year %>%
  filter(term=="YEARSC") %>%
  arrange(estimate) %>%
  mutate(SPECIES=factor(SPECIES, levels=SPECIES)) %>%
  mutate(COMMON_NAME=factor(COMMON_NAME, levels=COMMON_NAME)) 

inds<-which(mod_temp_lists %>% map(3)!="error")
allcoefs_temp <- mod_temp_lists[inds] %>% map(1) %>% dplyr::bind_rows() %>%
  left_join(spp_all)

tempvar="JFMAM"
TEMPeffects<-allcoefs_temp %>% 
  filter(term==tempvar) %>%
  arrange(estimate) %>%
  mutate(SPECIES=factor(SPECIES, levels=SPECIES)) %>%
  mutate(COMMON_NAME=factor(COMMON_NAME, levels=COMMON_NAME)) 

ggplot(YEAReffects, aes(x=estimate,y=COMMON_NAME)) +
  geom_point() +
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(xintercept=0,color="gray",lty=2)


ggplot(TEMPeffects, aes(x=estimate,y=COMMON_NAME)) +
  geom_point() +
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(xintercept=0,color="gray",lty=2)


ALLeffects<-bind_rows(TEMPeffects, YEAReffects) #%>%
  #arrange(term, estimate) %>%
#  mutate(SPECIES=factor(SPECIES, levels=SPECIES)) 
#  mutate(SPECIES = fct_reorder(SPECIES, desc(val))) %>%

ALLeffects<-ALLeffects %>%
  mutate(ispos = estimate>0,
         signif = sign(conf.low)==sign(conf.high)) %>%
  mutate(slope = case_when(ispos == TRUE ~ 'positive',
                           ispos == FALSE ~ 'negative')) %>%
  mutate(significance = case_when(signif ==TRUE ~ 'p < 0.05',
                                  signif ==FALSE ~ 'NS')) %>%
  mutate(term=case_match(term,
                    "JFMAM" ~ 'JFMAM Temp',
                    "YEARSC" ~ 'Year'))



#coloring by signif/slope

# color palette
col1<-'#F87505'
col2<-'#048F96'

p1<-ggplot(ALLeffects, aes(x=estimate,y=COMMON_NAME, color=slope, alpha=significance)) +
  geom_point() +
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(xintercept=0,color="gray",lty=2) +
  facet_wrap(~term,scales="free_x") +
  ylab("Species") +
  xlab("Beta") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_color_manual(values = c("negative" = col1, "positive" = col2))


p1

# get x axis limits from p1 to apply to histogram below
xrange1 <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
xrange2 <- ggplot_build(p1)$layout$panel_params[[2]]$x.range

# include below this the histogram of effects sizes
 

# 
# p2<-ggplot(ALLeffects, aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
#   geom_histogram(show.legend=FALSE,  boundary=0, bins=9)+
#   facet_wrap(~term,scales="free_x") +
#   geom_vline(xintercept=0,color="gray",lty=2) +
#   ylab("Count") +
#   xlab("Beta") +
#   scale_alpha_discrete(range=c(0.1, 0.6)) +
#   scale_color_manual(values = c("negative" = col1, "positive" = col2)) +
#   scale_fill_manual(values = c("negative" = col1, "positive" = col2))
# 
# p2
# p3<- p2 + facetted_pos_scales(x=list(
#   term == "JFMAM Temp" ~ scale_x_continuous(limits=xrange1, expand=c(0,0)),
#   term == "Year" ~ scale_x_continuous(limits=xrange2, expand = c(0,0))
# ))
# p3
#grid.arrange(p1, p2, heights = c(2,1))

if(famname=="gauss"){
  binw1<-0.5
  binw2<-0.05
}
if(famname=="genGam"){
  binw1<-0.05
  binw2<-0.005
}
if(famname=="logN"){
  binw1<-0.05
  binw2<-0.005
}



#plot histograms in different layers to control binwidth (from https://groups.google.com/g/ggplot2/c/aQQ2hTYRQF8)
p2b<-ggplot(ALLeffects, aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
  geom_histogram(data = ALLeffects %>% filter(term=="JFMAM Temp"), show.legend=FALSE,  boundary=0, binwidth=binw1)+
  geom_histogram(data = ALLeffects %>% filter(term=="Year"), show.legend=FALSE,  boundary=0, binwidth=binw2)+
    facet_wrap(~term,scales="free_x") +
  geom_vline(xintercept=0,color="gray",lty=2) +
  ylab("Count") +
  xlab("Beta") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_color_manual(values = c("negative" = col1, "positive" = col2)) +
  scale_fill_manual(values = c("negative" = col1, "positive" = col2))


p2b
#Use x-axis limits from top panels.
p3<- p2b + facetted_pos_scales(x=list(
  term == "JFMAM Temp" ~ scale_x_continuous(limits=xrange1, expand=c(0,0)),
  term == "Year" ~ scale_x_continuous(limits=xrange2, expand = c(0,0))
))
p3

ggarrange(p1,p3,ncol=1,heights = c(3,1),align="v")

#try ggh4x package to add unique scales to histograms


ggsave(here("Figures",paste0("Temp_Year_Slopes_GOA_",famname,"_",today(),".png")),width=9,height=8)
ggsave(here("Figures",paste0("Temp_Year_Slopes_GOA_",famname,"_",today(),".pdf")),width=9,height=8)

############ EBS ########################
############ EBS ########################
############ EBS ########################
############ EBS ########################

famname<-"gauss"
#load files with saved model run summaries to extract coefficients

if(famname=="gauss"){

  date<-"2025-12-03"
  
  myfoldername="ebs_TempModelFits_k5_noSP_gauss_JFMAM"
  mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
  myfoldername="ebs_LinearYearModelFits_k5_noSP_gauss"
  mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
}

if(famname=="genGam"){
  ## GenGamma models
  #load files with saved model run summaries to extract coefficients
  date<-"2025-08-11"
  myfoldername="ebs_TempModelFits_k5_noSP_genGam_JFMAM"
  mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
  myfoldername="ebs_LinearYearModelFits_k5_noSP_genGam"
  mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
}

if(famname=="logN"){
  ## logNormal models
  #load files with saved model run summaries to extract coefficients
  date<-"2025-08-11"
  myfoldername="ebs_TempModelFits_k5_noSP_logN_JFMAM"
  mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
  myfoldername="ebs_LinearYearModelFits_k5_noSP_logN"
  mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
}

inds<-which(mod_year_lists %>% map(3)!="error")
allcoefs_year <- mod_year_lists[inds] %>% map(1) %>% dplyr::bind_rows() %>%
  left_join(spp_all)

YEAReffects<-allcoefs_year %>%
  filter(term=="YEARSC") %>%
  arrange(estimate) %>%
  mutate(SPECIES=factor(SPECIES, levels=SPECIES)) %>%
  mutate(COMMON_NAME=factor(COMMON_NAME, levels=COMMON_NAME)) 

inds<-which(mod_temp_lists %>% map(3)!="error")
allcoefs_temp <- mod_temp_lists[inds] %>% map(1) %>% dplyr::bind_rows() %>%
  left_join(spp_all)

tempvar="JFMAM"
TEMPeffects<-allcoefs_temp %>% 
  filter(term==tempvar) %>%
  arrange(estimate) %>%
  mutate(SPECIES=factor(SPECIES, levels=SPECIES)) %>%
  mutate(COMMON_NAME=factor(COMMON_NAME, levels=COMMON_NAME)) 

ggplot(YEAReffects, aes(x=estimate,y=COMMON_NAME)) +
  geom_point() +
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(xintercept=0,color="gray",lty=2)


ggplot(TEMPeffects, aes(x=estimate,y=COMMON_NAME)) +
  geom_point() +
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(xintercept=0,color="gray",lty=2)


ALLeffects<-bind_rows(TEMPeffects, YEAReffects) #%>%
#arrange(term, estimate) %>%
#  mutate(SPECIES=factor(SPECIES, levels=SPECIES)) 
#  mutate(SPECIES = fct_reorder(SPECIES, desc(val))) %>%

ALLeffects<-ALLeffects %>%
  mutate(ispos = estimate>0,
         signif = sign(conf.low)==sign(conf.high)) %>%
  mutate(slope = case_when(ispos == TRUE ~ 'positive',
                           ispos == FALSE ~ 'negative')) %>%
  mutate(significance = case_when(signif ==TRUE ~ 'p < 0.05',
                                  signif ==FALSE ~ 'NS')) %>%
  mutate(term=case_match(term,
                         "JFMAM" ~ 'JFMAM Temp',
                         "YEARSC" ~ 'Year'))



#coloring by signif/slope

# color palette
col1<-'#F87505'
col2<-'#048F96'

p1<-ggplot(ALLeffects, aes(x=estimate,y=COMMON_NAME, color=slope, alpha=significance)) +
  geom_point() +
  geom_pointrange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(xintercept=0,color="gray",lty=2) +
  facet_wrap(~term,scales="free_x") +
  ylab("Species") +
  xlab("Beta") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_color_manual(values = c("negative" = col1, "positive" = col2))


p1

# get x axis limits from p1 to apply to histogram below
xrange1 <- ggplot_build(p1)$layout$panel_params[[1]]$x.range
xrange2 <- ggplot_build(p1)$layout$panel_params[[2]]$x.range

# include below this the histogram of effects sizes


# 
# p2<-ggplot(ALLeffects, aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
#   geom_histogram(show.legend=FALSE,  boundary=0, bins=9)+
#   facet_wrap(~term,scales="free_x") +
#   geom_vline(xintercept=0,color="gray",lty=2) +
#   ylab("Count") +
#   xlab("Beta") +
#   scale_alpha_discrete(range=c(0.1, 0.6)) +
#   scale_color_manual(values = c("negative" = col1, "positive" = col2)) +
#   scale_fill_manual(values = c("negative" = col1, "positive" = col2))
# 
# p2
# p3<- p2 + facetted_pos_scales(x=list(
#   term == "JFMAM Temp" ~ scale_x_continuous(limits=xrange1, expand=c(0,0)),
#   term == "Year" ~ scale_x_continuous(limits=xrange2, expand = c(0,0))
# ))
# p3
#grid.arrange(p1, p2, heights = c(2,1))

#plot histograms in different layers to control binwidth (from https://groups.google.com/g/ggplot2/c/aQQ2hTYRQF8)

if(famname=="gauss"){
  binw1<-0.5
  binw2<-0.05
}

if(famname=="genGam"){
  binw1<-0.1
  binw2<-0.005
}
if(famname=="logN"){
  binw1<-0.1
  binw2<-0.005
}

p2b<-ggplot(ALLeffects, aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
  geom_histogram(data = ALLeffects %>% filter(term=="JFMAM Temp"), show.legend=FALSE,  boundary=0, binwidth=binw1)+ #0.5 for gauss
  geom_histogram(data = ALLeffects %>% filter(term=="Year"), show.legend=FALSE,  boundary=0, binwidth=binw2)+  #0.05 for gauss
  facet_wrap(~term,scales="free_x") +
  geom_vline(xintercept=0,color="gray",lty=2) +
  ylab("Count") +
  xlab("Beta") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_color_manual(values = c("negative" = col1, "positive" = col2)) +
  scale_fill_manual(values = c("negative" = col1, "positive" = col2))


p2b
#Use x-axis limits from top panels to specify scales for histograms
p3<- p2b + facetted_pos_scales(x=list(
  term == "JFMAM Temp" ~ scale_x_continuous(limits=xrange1, expand=c(0,0)),
  term == "Year" ~ scale_x_continuous(limits=xrange2, expand = c(0,0))
))
p3

ggarrange(p1,p3,ncol=1,heights = c(1.2,1),align="v")


ggsave(here("Figures",paste0("Temp_Year_Slopes_EBS_",famname,"_",today(),".png")),width=9,height=4.5)
ggsave(here("Figures",paste0("Temp_Year_Slopes_EBS_",famname,"_",today(),".pdf")),width=9,height=4.5)

       