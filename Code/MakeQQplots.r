# Accompanying code for:
# Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
# temperature in Alaska's marine fishes. Global Change Biology.
# 
# Plot DHARMa residuals for lists of sdmTMB models
# Generates SI figures S10-S13.

source("Code/MyFunctions.r")

today<-Sys.Date()
spp_all<-read.csv(here("Data","SpeciesNamesCodes.csv"))

dharma_res_fun<-function(mymod, nsims){
  mymod %>% 
    simulate(type='mle-mvn', nsim = nsims) %>% 
    dharma_residuals(mymod,return_DHARMa = FALSE,plot=FALSE)
}

make_QQ_ggplots<-function(foldername, modelname, nsims){
  allfiles<-list.files(here("Results","ModelResults",foldername),full.names=TRUE)
  sppnames<-list.files(here("Results","ModelResults",foldername))
  sppnames<-sppnames[!sppnames == "AllSppLists.rds"]
  allfiles<-allfiles[!grepl("AllSppLists.rds", allfiles)]
  rds_list <- lapply(allfiles, readRDS)
  extracted_elements <- lapply(rds_list, function(x) x$MODEL)
  dharma_res_list<-map(extracted_elements,dharma_res_fun,nsims=nsims) #or use lapply
  names(dharma_res_list)<-sppnames
  dharma_res_df <- bind_rows(dharma_res_list, .id = "SPECIES_NAME") %>%
    mutate(SPECIES_NAME = str_extract(SPECIES_NAME, "^[^\\.]*")) %>%
    left_join(spp_all[,c("SPECIES_NAME","COMMON_NAME")])
  ggplot(dharma_res_df,aes(expected,observed)) +
    geom_point() +
    labs(x="Expected",y="Observed",title=paste0("DHARMa residuals from ",modelname)) +
    facet_wrap(~SPECIES_NAME,ncol=6) +
    geom_abline(slope=1,color="red") 
#  return(dharma_res_df)
}

######### GOA ###########

foldername<-"goa_YearFModelFits_k5_noSP_gauss_2025-11-17"
plot1<-make_QQ_ggplots(foldername, "Year Factor models, Gaussian, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot1, width=11, height=9)

foldername<-"goa_LinearYearModelFits_k5_noSP_gauss_2025-09-30"
plot2<-make_QQ_ggplots(foldername, "Linear Year models, Gaussian, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot2, width=11, height=9)

foldername<-"goa_TempModelFits_k5_noSP_gauss_JFMAM_2025-09-30"
plot3<-make_QQ_ggplots(foldername, "Temperature models, Gaussian, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot3, width=11, height=9)

foldername<-"goa_TempModelFits_k5_noSP_genGam_JFMAM_2025-09-26"
plot4<-make_QQ_ggplots(foldername, "Temperature models, Generalized Gamma, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot4, width=11, height=9)

foldername<-"goa_TempModelFits_k5_noSP_logN_JFMAM_2025-08-11"
plot5<-make_QQ_ggplots(foldername, "Temperature models, Lognormal, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot5, width=11, height=9)

foldername<-"goa_LinearYearModelFits_k5_noSP_genGam_2025-08-11"
plot6<-make_QQ_ggplots(foldername, "Linear Year models, Generalized Gamma, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot6, width=11, height=9)

foldername<-"goa_LinearYearModelFits_k5_noSP_logN_2025-08-11"
plot7<-make_QQ_ggplots(foldername, "Linear Year models, Lognormal, GOA", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot7, width=11, height=9)

##### EBS #########

foldername<-"ebs_YearFModelFits_k5_noSP_gauss_2025-12-03"
plot8<-make_QQ_ggplots(foldername, "Year Factor models, Gaussian, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot8, width=11, height=2.3)

foldername<-"ebs_LinearYearModelFits_k5_noSP_gauss_2025-12-03"
plot9<-make_QQ_ggplots(foldername, "Linear Year models, Gaussian, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot9, width=11, height=2.3)

foldername<-"ebs_TempModelFits_k5_noSP_gauss_JFMAM_2025-12-03"
plot10<-make_QQ_ggplots(foldername, "Temperature models, Gaussian, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot10, width=11, height=2.3)

foldername<-"ebs_TempModelFits_k5_noSP_genGam_JFMAM_2025-12-03"
plot11<-make_QQ_ggplots(foldername, "Temperature models, Generalized Gamma, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot11, width=11, height=2.3)

foldername<-"ebs_TempModelFits_k5_noSP_logN_JFMAM_2025-12-03"
plot12<-make_QQ_ggplots(foldername, "Temperature models, Lognormal, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot12, width=11, height=2.3)

foldername<-"ebs_LinearYearModelFits_k5_noSP_genGam_2025-12-03"
plot13<-make_QQ_ggplots(foldername, "Linear Year models, Generalized Gamma, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot13, width=11, height=2.3)

foldername<-"ebs_LinearYearModelFits_k5_noSP_logN_2025-12-03"
plot14<-make_QQ_ggplots(foldername, "Linear Year models, Lognormal, EBS", nsims=250)
ggsave(here("Figures","QQplots",paste0(foldername,".png")),
       plot=plot14, width=11, height=2.3)





