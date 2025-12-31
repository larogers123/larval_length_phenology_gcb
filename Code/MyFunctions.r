## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
### Series of functions to call to prep data, create custom plots, etc
### 
### 
library(sdmTMB)
library(mgcv)
library(gratia)
library(visreg)
library(dplyr)
library(ggplot2)
library(purrr)
library(lubridate)
library(here)
#library(sdmTMBextra)
library(tidyr)
library(bayesdfa)
library(ggrepel)
library(stringr)

#library(INLA)
#
source("Code/PrepMapData.r")

options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
theme_set(theme_bw())

#options(ggplot2.continuous.colour = "scale_color_nmfs")
#options(ggplot2.continuous.fill = "scale_fill_nmfs")
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

prepSpeciesData<-function(dat,spp,utm_crs=32604){
    sppdat<-dat %>%
    filter(SPECIES_NAME==spp) %>%
    rename(DENSITY = LARVALCATCHPER10M2) %>%
    mutate(DATETIME = ymd_hm(GMT_DATE_TIME_TXT)) %>%
    filter(!is.na(DENSITY)) %>%
    mutate(YEARF= factor(YEAR)) %>%
    mutate(YEARSC= YEAR - 2000) %>%
    mutate(HAUL_ID = factor(HAUL_ID)) %>%
    mutate(MESH = factor(MESH)) %>%
    mutate(YDAYSC=YDAY-145) %>%
    mutate(LOG_LENGTH=log(CORRECTED_LENGTH)) %>%
      rowwise() %>%
    mutate(JFMAM=mean(c(JAN,FEB,MAR,APR,MAY)),
           prevNOVDEC=mean(c(prevNOV,prevDEC)),
           JANFEB=mean(c(JAN,FEB)),
           MARAPR=mean(c(MAR,APR)),
           MAYJUN=mean(c(MAY,JUN)),
           JULAUG=mean(c(JUL,AUG))) %>%
    group_by(YEAR) %>%
#    filter(is.na(APR)==FALSE) %>%
    filter(n()>=10) %>%
    ungroup() %>%
    droplevels() %>%
    add_utm_columns(ll_names = c("LON", "LAT"),utm_crs=utm_crs)
  return(as.data.frame(sppdat))
}


plotLengthMaps<-function(dat,lme){
  akproj<-switch(lme,
                 goa=ak_coast_goa_proj,
                 ebs=ak_coast_ebs_proj)
  xlims<-switch(lme,
                goa=c(0, 1364000 - 300000),
                ebs=c(-255000, 940000))
  ylims<-switch(lme,
                goa=c(5940000, 6700000),
                ebs=c(5880000, 6760000))
  ggplot() +
    geom_point(data=dat,aes(x=X * 1000, y= Y*1000, size = DENSITY, colour = CORRECTED_LENGTH),alpha = 0.3) +
    facet_wrap(~YEAR) +
#    coord_fixed() +
    geom_sf() +
    geom_sf(data=akproj) +
    xlim(xlims) +
    ylim(ylims) +
    labs(x = "Longitude", y = "Latitude")  +
    ggtitle(dat$SPECIES_NAME[1], subtitle=dat$COMMON_NAME[1])
}



plotSpatialPreds<-function(preds,toplot,lme){ #one of est,epsilon_st,omega_s
  akproj<-switch(lme,
                 goa=ak_coast_goa_proj,
                 ebs=ak_coast_ebs_proj)
  xlims<-switch(lme,
                goa=c(0, 1364000 - 300000),
                ebs=c(-255000, 940000))
  ylims<-switch(lme,
                goa=c(5940000, 6700000),
                ebs=c(5880000, 6760000))
  ggplot() +
    geom_raster(data=preds, aes(x=X * 1000, y=Y*1000, fill = .data[[toplot]])) + #epsilon_st is spatiotemp field
    facet_wrap(~YEAR) +
 #   scale_fill_gradient2() +
    geom_sf() +
    geom_sf(data=akproj) +
    theme_light() +
    xlim(xlims) +
    ylim(ylims) +
    labs(x = "Longitude", y = "Latitude")  
}

# DFA Trends (modified from bayesdfa package)
plot_trends_LR <- function (rotated_modelfit, allyears = NULL, datayears = NULL) 
{
  rotated <- rotated_modelfit
  df <- dfa_trends(rotated, years = allyears) %>% 
    filter(time %in% datayears) %>%
    full_join(data.frame("time"=min(datayears):max(datayears)))
  p1 <- ggplot(df, aes_string(x = "time", y = "estimate")) + 
    geom_point() +
    geom_linerange(aes_string(ymin = "lower", ymax = "upper"), 
                   alpha = 0.4) + 
    geom_line() + ## change to geom_path to connect across missing points
    #    facet_wrap("trend_number") + 
    xlab("Year") + 
    ylab("DFA trend") 
  p1
}


# DFA Factor Loadings (modified from bayesdfa package)
# modified but may only work with one trend
plot_loadings_LR <- function (rotated_modelfit, names = NULL, facet = FALSE, violin = TRUE, 
                              conf_level = 0.95) 
{
  v <- dfa_loadings(rotated_modelfit, summary = FALSE, names = names, 
                    conf_level = conf_level)
  df <- dfa_loadings(rotated_modelfit, summary = TRUE, names = names, 
                     conf_level = conf_level)
  if (!violin) {
    p1 <- ggplot(df, aes_string(x = "name", y = "median", 
              col = "trend", alpha = "prob_diff0")) + 
      geom_point(size = 3, 
      position = position_dodge(0.3)) + 
      geom_errorbar(aes_string(ymin = "lower", ymax = "upper"), position = position_dodge(0.3), width = 0) +
      geom_hline(yintercept = 0, lty = 2) + 
      coord_flip() + xlab("Time Series") + ylab("Loading")
  }
  if (violin) {
    p1 <- ggplot(v, aes_string(x = 'name', y = 'loading', 
                               fill = 'prob_diff0')) + 
      geom_violin(color=NA) + 
      geom_hline(yintercept = 0, lty = 2) + coord_flip() + 
      scale_fill_distiller("Probability different \nfrom zero",
                           type="seq",palette="Blues",aesthetics="fill",
                           direction = 1,limits=c(0.5,1)) + 
      xlab("Species") + ylab("Loading") +
      scale_x_discrete(limits = rev)#+
    #  theme(legend.position = "top")
  }
  if (facet) {
    p1 <- p1 + facet_wrap(~trend, scales = "free_x")
  }
  p1
}

