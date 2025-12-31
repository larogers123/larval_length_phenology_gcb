## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
##
## Length DFA sample size testing
## Fit DFA to estimates of size-at-date for larvae from the GOA. 
## Subset data to include only the longest timeseries, versus all,
## and compare trends. 
##

source("Code/MyFunctions.r")
library(ggpubr)
today<-Sys.Date()

flip_trend<-function(trend){
  r <- trend %>% 
    mutate(estimate=-1*estimate,
           lower=-1*lower,
           upper=-1*upper)
  r
}

# #Load year effects from Year Factor models

date<-"2025-11-17"
myregion<-"goa"
myfoldername=paste0(myregion,"_YearFModelFits_k5_noSP_gauss")
mod_YearF_lists_goa<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))

## Species names and common names table (all)
spp_all<-read.csv(here("Data","SpeciesNamesCodes.csv"))


## Set DFA params 
chains <- 4
iter <- 4000

## GOA long format (but don't use weights)

# Merge Year Effects into one dataframe, prep data for DFA
inds<-which(mod_YearF_lists_goa %>% map(1)!="error")
allfits_YearF_g <- mod_YearF_lists_goa[inds] %>% map(6) %>% dplyr::bind_rows()

yrs_g<-as.numeric(as.character(unique(allfits_YearF_g$YEARF)))
spp_g<-unique(allfits_YearF_g$SPECIES)

yrvec_g<-min(yrs_g):max(yrs_g)

#long format data
yrcomp_long<-allfits_YearF_g %>%
  mutate(Year=as.numeric(as.character(YEARF))) %>%
  mutate(se=(visregUpr-visregLwr)/3.92) %>% ## 95%CI is +/- 1.96 SE
  mutate(weights = (1 / se)^2) %>%
  mutate(YearNum=Year-min(Year)+1) %>%
  select(YearNum,Year,visregFit,SPECIES, se,weights) %>%
  mutate(SPECIES=as.factor(SPECIES)) %>%
  group_by(SPECIES) %>%
  mutate(nYrs=n()) %>%
  rename(ts=SPECIES,obs=visregFit,se=se,time=YearNum) %>%
  as.data.frame()

sort(table(yrcomp_long$ts))

#time series length cutoffs for comparison
yrcomp_10y<-yrcomp_long
yrcomp_15y<-yrcomp_long %>% filter(nYrs>=15) %>% mutate(ts=droplevels(ts))
yrcomp_20y<-yrcomp_long %>% filter(nYrs>=20) %>% mutate(ts=droplevels(ts))
yrcomp_25y<-yrcomp_long %>% filter(nYrs>=25) %>% mutate(ts=droplevels(ts))

set.seed(103)  #103 gives positive loading for GOA. 

# base DFA includes all species
f1_full <- fit_dfa(
  y = yrcomp_long, num_trends = 1, scale="zscore",data_shape = "long",
  iter = iter, chains = chains, thin = 1)#, inv_var_weights = "weights")
f10 <- f1_full

f15 <- fit_dfa(
  y = yrcomp_15y, num_trends = 1, scale="zscore",data_shape = "long",
  iter = iter, chains = chains, thin = 1)#, inv_var_weights = "weights")

f20 <- fit_dfa(
  y = yrcomp_20y, num_trends = 1, scale="zscore",data_shape = "long",
  iter = iter, chains = chains, thin = 1)#, inv_var_weights = "weights")

f25 <- fit_dfa(
  y = yrcomp_25y, num_trends = 1, scale="zscore",data_shape = "long",
  iter = iter, chains = chains, thin = 1)#, inv_var_weights = "weights")

######

r_10 <- rotate_trends(f10)
r_15 <- rotate_trends(f15)
r_20 <- rotate_trends(f20)
r_25 <- rotate_trends(f25)

trends_10<-dfa_trends(r_10,years=min(yrcomp_10y$Year):max(yrcomp_10y$Year)) %>%
  filter(time %in% unique(yrcomp_10y$Year))
trends_15<-dfa_trends(r_15,years=min(yrcomp_15y$Year):max(yrcomp_15y$Year)) %>%
  filter(time %in% unique(yrcomp_15y$Year))
trends_20<-dfa_trends(r_20,years=min(yrcomp_20y$Year):max(yrcomp_20y$Year)) %>%
  filter(time %in% unique(yrcomp_20y$Year))
trends_25<-dfa_trends(r_25,years=min(yrcomp_25y$Year):max(yrcomp_25y$Year)) %>%
  filter(time %in% unique(yrcomp_25y$Year))



## For MS supplement, show trend as estimated from all TS, and then from
## stricter subsets.
## 
## 


trends_10$type<-">=10 yrs"
trends_15$type<-">=15 yrs"
trends_20$type<-">=20 yrs"
trends_25$type<-">=25 yrs"

#trends_10<-flip_trend(trends_10)
#trends_15<-flip_trend(trends_15)
#trends_25<-flip_trend(trends_25)


alltrends<-bind_rows(trends_10,trends_15) %>% 
  bind_rows(trends_20) %>%
  bind_rows(trends_25)

ggplot(alltrends,aes(time, estimate,fill=type,color=type)) +
  geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), alpha=0.2,color="NA") +
  geom_line() +
  geom_point(alpha=0.6) +
  xlab("Year") +
  ylab("DFA Trend") +
  scale_fill_brewer(palette = "BrBG")+
  scale_color_brewer(palette = "BrBG") +
  labs(color="Timeseries length") +
  labs(fill="Timeseries length")

ggsave(here("Figures","DFA","DFA_TimeseriesLengthComp_noWts.png"),width=8,height=5)

trends_wide<-pivot_wider(alltrends,names_from=type,values_from=c(estimate,lower,upper))

cor(trends_wide[,3:6]) # ALL are >0.98

comptrends<-trends_wide %>% left_join(SST_GOA,join_by(time==YEAR))

cor(comptrends$JFMAM,comptrends$`estimate_>=10 yrs`,use="p")
cor(comptrends$JFMAM,comptrends$`estimate_>=15 yrs`,use="p")
cor(comptrends$JFMAM,comptrends$`estimate_>=20 yrs`,use="p")
cor(comptrends$JFMAM,comptrends$`estimate_>=25 yrs`,use="p")


plot(comptrends$JFMAM,comptrends$`estimate_>=10 yrs`)


