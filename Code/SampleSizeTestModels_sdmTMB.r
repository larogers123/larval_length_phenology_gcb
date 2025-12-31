## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Randomly subset data-rich species length data to test how reduced
## sample size affects estimates of temperature effects on larval size (JFMAM temp model).
## 
## 
source("Code/MyFunctions.r")

#Prep GOA data
goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")
goa_ebs<-bind_rows(goa,ebs)

myregion<-"goa"
regiondat<-goa 

myspp<-unique(regiondat$SPECIES_NAME)

#Filter out species that are consistently not converging
#Also leave out Capelin as they are previous year's spawn

myspp<-myspp[!myspp %in% c("Anoplopoma fimbria","Lumpenella longirostris","Ronquilus jordani","Mallotus villosus")]
utm_crs <-32604

sp_dat_list<-lapply(X=myspp,FUN=prepSpeciesData,dat=regiondat,utm_crs=utm_crs)
names(sp_dat_list)<-myspp

######
###### Model fitting functions
######


######################
### Year Factor Model 
######################


fit_YearF_mod_func_sim<-function(sp,foldername,myregion,myfamily){
  
  #Check for folder to save model fits, create if not existing  
  today<-Sys.Date()
  mypath<-here("Results","ModelResults",paste0(foldername,"_",today)) 
  if(!dir.exists(mypath)) dir.create(mypath)
  
  sppdat<-simdatlist[[sp]]
  
  mesh_cutoff<-switch(myregion,
                      goa=30,
                      ebs=40)
  
  mesh <- make_mesh(sppdat, xy_cols = c("X", "Y"), cutoff = mesh_cutoff)
  #plot(mesh) 
  
  
  fit2 <- try(sdmTMB(
    CORRECTED_LENGTH ~ 0 + s(YDAYSC, k=5) + YEARF + MESH, # 
    data = sppdat,
    mesh = mesh,
    spatial = "off",
    spatiotemporal = "IID",
    family=myfamily,
    time="YEAR"
  )
  )
  if(inherits(fit2, "try-error")){
    print(paste("model 2 error for species",sp))
    reslist<-list("error","error","error")
  } else { if (sanity(fit2)$all_ok ==F) fit2<-try(run_extra_optimization(fit2))}
  if(inherits(fit2, "try-error")){
    print(paste("model 2 error for species",sp))
    reslist<-list("error","error","error")} else {

      COEFS<-tidy(fit2,conf.int=TRUE)
      COEFS$SPECIES<-sp
      res<-residuals(fit2)
      
      vr<-visreg(fit2,xvar="YEARF",cond=list(YDAYSC=0),plot=FALSE)
      partials<-vr$res
      partials$SPECIES<-sp
      fit<-vr$fit
      fit$SPECIES<-sp    
      
      # Save model object separately for each species
      # Compile summary stats into an object to save
      # Function will output visreg, AIC, summary stats, but not model object.
      modreslist<-list("MODEL"=fit2,"COEFS"=COEFS,"AIC"=AIC(fit2),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"MESH_CUTOFF"=mesh_cutoff)
      saveRDS(modreslist,file=here("Results","ModelResults",paste0(foldername,"_",today),paste0(sp,".rds")))
      reslist<-list("COEFS"=COEFS,"AIC"=AIC(fit2),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"MESH_CUTOFF"=mesh_cutoff)
    }
  return(reslist)
}

######################
#### Temperature effect model
######################

fit_temp_mod_func_sim<-function(sp,foldername,myregion,tempvar,myfamily){
  
  #Check for folder to save model fits, create if not existing  
  today<-Sys.Date()
  mypath<-here("Results","ModelResults",paste0(foldername,"_",today)) 
  if(!dir.exists(mypath)) dir.create(mypath)
  
  mesh_cutoff<-switch(myregion,
                      goa=30,
                      ebs=40)
  
  sppdat<-simdatlist[[sp]]
  mesh <- make_mesh(sppdat, xy_cols = c("X", "Y"), cutoff = mesh_cutoff)
  
  myformula<-as.formula(paste0("CORRECTED_LENGTH ~ 0 + s(YDAYSC, k=5) + MESH +",tempvar)) 
  
  fit4 <- try(
    sdmTMB(
      myformula, 
      data = sppdat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "iid",
      family=myfamily,
      time="YEAR",
      anisotropy = FALSE
    )
  )
  
  if(inherits(fit4, "try-error")){
    print(paste("model 4 error for species",sp))
    reslist<-list("error","error","error") 
    
    #} else if (sanity(fit4)$all_ok ==F) {
    #  print(paste("model 4 failed to converge for species", sp)) ### need to exit, or print error, or go to next spp
  } else { if (sanity(fit4)$all_ok ==F) fit4<-try(run_extra_optimization(fit4))}
  if(inherits(fit4, "try-error")){
    print(paste("model 4 error for species",sp))
    reslist<-list("error","error","error")} else {
 
      COEFS<-tidy(fit4,conf.int=TRUE)
      COEFS$SPECIES<-sp
      res<-residuals(fit4)
      
      vr<-visreg(fit4,xvar=tempvar,cond=list(YDAYSC=0),plot=FALSE)
      partials<-vr$res
      partials$SPECIES<-sp
      fit<-vr$fit
      fit$SPECIES<-sp
      
      # Save model object separately for each species
      # Compile summary stats into an object to save
      # Function will output visreg, AIC, summary stats, but not model object.
      modreslist<-list("MODEL"=fit4,"COEFS"=COEFS,"AIC"=AIC(fit4),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"FORMULA"=myformula,"MESH_CUTOFF"=mesh_cutoff)
      saveRDS(modreslist,file=here("Results","ModelResults",paste0(foldername,"_",today),paste0(sp,".rds")))
      reslist<-list("COEFS"=COEFS,"AIC"=AIC(fit4),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"FORMULA"=myformula,"MESH_CUTOFF"=mesh_cutoff)
    }
  return(reslist)
}


###
### Use GOA Pollock as "data rich" species. Subset to lower sample sizes. 
### Randomly choose 10 samples/year for 10 years as the most strict test.
### 

fulldat<-sp_dat_list$`Gadus chalcogrammus`

## Function to randomly subsample n_size samples within each year

subsetdata<-function(sim,fulldat,n_size){
  subdat<- fulldat %>%
  group_by(YEAR) %>%
  sample_n(size=n_size, replace=FALSE) %>%
  ungroup() %>%
  mutate(SIM=sim)
  return(as.data.frame(subdat))
}

# Make a list of subsampled datasets and run models on list
set.seed(999)
sims<-1:250 #need at least 250 to get 100 converged YearF models
simdatlist<-lapply(X=sims,FUN=subsetdata,fulldat=fulldat,n_size=10)

## Fit model to subset data
myfamily<-gaussian(link = "identity"); famname<-"gauss"

myfoldername=paste0("Sim_",myregion,"_YearF_",famname)
myYearFfoldername<-myfoldername
mod_YearF_lists <- lapply(sims, fit_YearF_mod_func_sim,foldername=myfoldername,myregion=myregion,myfamily=myfamily)
#names(mod_YearF_lists)<-myspp#[sp_tofit]
saveRDS(mod_YearF_lists,file=here("Results","ModelResults",paste0(myfoldername,"_",today()),"AllSppLists.rds"))

tempvar="JFMAM"
myfoldername=paste0("Sim_",myregion,"_TempModelFits_k5_noSP_",famname,"_",tempvar)
myTempfoldername<-myfoldername
mod_temp_lists <- lapply(sims, fit_temp_mod_func_sim,foldername=myfoldername,myregion=myregion,tempvar=tempvar,myfamily=myfamily)
#names(mod_temp_lists)<-myspp#[sp_tofit]
saveRDS(mod_temp_lists,file=here("Results","ModelResults",paste0(myfoldername,"_",today()),"AllSppLists.rds"))

### To compare to model of full dataset, look at distribution of estimated effect sizes, including direction and signif.

# Read in yearF model fit to full dataset.
fullyearF<-readRDS("Results/ModelResults/goa_YearFModelFits_k5_noSP_gauss_2025-09-30/Gadus chalcogrammus.rds")
fullyearF$FIT

# Read in temperature model fit to full dataset.
fulltemp<-readRDS("Results/ModelResults/goa_TempModelFits_k5_noSP_gauss_JFMAM_2025-09-30/Gadus chalcogrammus.rds")
fulltempcoef<-fulltemp$COEFS$estimate[3]





### Plotting

inds<-which(mod_YearF_lists %>% map(1)!="error")
allcoefs_YearF <- mod_YearF_lists[inds] %>% map(1) %>% dplyr::bind_rows()
#allpreds <- mod_YearF_lists[inds] %>% map(3) %>% dplyr::bind_rows()
allpartials_YearF <- mod_YearF_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_YearF <- mod_YearF_lists[inds] %>% map(6) %>% dplyr::bind_rows()

allfits_YearF <-allfits_YearF %>%
  mutate(YEAR=as.numeric(levels(YEARF))[YEARF]) 
fullyearfits <- fullyearF$FIT %>%
  mutate(YEAR=as.numeric(levels(YEARF))[YEARF]) 

# Reshape to wide format for correlations
allfits_wide<-allfits_YearF %>%
  select(YEAR,visregFit,SPECIES) %>%
  pivot_wider(names_from=SPECIES,values_from = visregFit)

yrcors<-cor(allfits_wide[,-1])
hist(yrcors[upper.tri(yrcors)])
mean(yrcors[upper.tri(yrcors)])
median(yrcors[upper.tri(yrcors)])


ggplot(allfits_YearF,aes(YEAR,visregFit)) +
  geom_point() +
  geom_point(data=fullyearfits,aes(YEAR,visregFit),color="orange",shape=17,size=2) +
  geom_errorbar(data=fullyearfits,aes(x=YEAR,ymin=visregLwr,ymax=visregUpr),color="orange") +
#  facet_wrap(~SPECIES) +
  ylab("Larval Length (mm)") +
  xlab("Year")

allcoefs_YearF<-allcoefs_YearF %>%
  filter(! term %in% c("MESH505","sYDAYSC"))

ggplot(allcoefs_YearF,aes(term,estimate)) +
  geom_point() +
#  geom_point(data=fullyearfits,aes(YEAR,visregFit),color="orange",shape=17,size=2) +
#  geom_errorbar(data=fullyearfits,aes(x=YEAR,ymin=visregLwr,ymax=visregUpr),color="orange") +
  #  facet_wrap(~SPECIES) +
  ylab("Larval Length (mm)") +
  xlab("Year")




## Temp model plots

# color palette
col1<-'#F87505'
col2<-'#048F96'

## need to exclude species with errors AND without SEs
inds<-which(mod_temp_lists %>% map(3)!="error")
allcoefs_temp <- mod_temp_lists[inds] %>% map(1) %>% dplyr::bind_rows()
#allpreds_temp <- mod_temp_lists[inds] %>% map(3) %>% dplyr::bind_rows()
allpartials_temp <- mod_temp_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_temp <- mod_temp_lists[inds] %>% map(6) %>% dplyr::bind_rows()

#tempvar="JFMAM"
TEMPeffects<-allcoefs_temp %>% filter(term==tempvar) %>% arrange(SPECIES)

histplot<-TEMPeffects %>%
  mutate(ispos = estimate>0, 
         signif = sign(conf.low)==sign(conf.high)) %>%
  mutate(slope = case_when(ispos == TRUE ~ 'positive',
                           ispos == FALSE ~ 'negative')) %>%
  mutate(significance = case_when(signif ==TRUE ~ 'p < 0.05',
                                  signif ==FALSE ~ 'NS')) %>%
  filter(signif %in% c("TRUE","FALSE"))

histplot<-histplot[1:200,]

toplot<- histplot %>%  select(SPECIES,slope,significance) %>%
  full_join(allfits_temp)


ggplot() +
  geom_point(data=allpartials_temp,aes(x=.data[[tempvar]],y=visregRes),alpha=0.2,size=0.2,color="gray") +
  geom_ribbon(data=toplot %>% filter(YDAYSC==0),
              aes(.data[[tempvar]], visregFit, ymin = visregLwr,
                  ymax = visregUpr, fill=slope, alpha=significance)) +
  geom_line(data=toplot %>% filter(YDAYSC==0), aes(.data[[tempvar]], visregFit,linetype=significance)) +
  xlab(paste0(tempvar," Temperature (C)")) +
  ylab("Larval Length (mm, May 25)") +
  facet_wrap(~SPECIES,scales="free_y") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_linetype_manual(values=c('NS' = "dashed", "p < 0.05" = "solid")) +
  scale_fill_manual(values = c("negative" = col1, "positive" = col2))


## Figure for supplement showing temperature effects histogram

ggplot(histplot,aes(x=estimate,fill=slope,alpha=significance)) +
  geom_histogram(binwidth=0.1,center=0.05,show.legend=TRUE,color=1)+
  xlab("Estimated temp. effect (mm/deg C)") +
  geom_vline(xintercept=0,linetype="dashed") +
  scale_y_continuous(expand=c(0,0))+
  #  scale_x_continuous(limits=c(-2.5,2.5))+
  theme_classic() +
#  scale_colour_manual(values = c("negative" = col1, "positive" = col2))+
  scale_fill_manual(values = c("negative" = col1, "positive" = col2)) +
  geom_vline(xintercept=fulltempcoef,linewidth=1.5,color="dark orange")

ggsave(here("Figures","Sims","TemperatureHists.png"))

table(histplot$signif)
