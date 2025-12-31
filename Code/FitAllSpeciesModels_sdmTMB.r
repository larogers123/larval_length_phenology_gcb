## Accompanying code for:
## Rogers, L.A., K.E. Axler, J.S. Bigman. Widespread phenological shifts with
## temperature in Alaska's marine fishes. Global Change Biology.
## 
## Fit models for all species and save results
## 

# Read in data prep and mapping functions, packages, graphics defaults
source("Code/MyFunctions.r")
today<-Sys.Date()


######################
### PREP DATA
######################

goa<-read.csv("Data/Processed/GOA_spring_d117to160_32spp_Jul2024.csv")
ebs<-read.csv("Data/Processed/EBS_spring_d98to162_11spp_Apr2024.csv")

goa_ebs<-bind_rows(goa,ebs)


######################
#### Select region
######################

myregion<-c("ebs","goa")[1]

if(myregion=="ebs"){ 
  regiondat<-ebs 
}else if(myregion=="goa"){ 
  regiondat<-goa 
}

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



# create a list with each item a species-specific dataframe for model fitting.
if(myregion=="ebs"){
  utm_crs <- 32603 
}else{
  utm_crs <-32604
}

sp_dat_list<-lapply(X=myspp,FUN=prepSpeciesData,dat=regiondat,utm_crs=utm_crs)
names(sp_dat_list)<-myspp

#summary of number of specimens
#
nspecs<-sp_dat_list %>% 
  map(dim) %>%
  sapply(function(x) x[1]) %>%
  sum()

#summary of number of hauls
#
nhauls<-sp_dat_list %>% 
  bind_rows() %>%
  summarize(nhauls=length(unique(HAUL_ID))) #2438 for EBS, 6571 for GOA

#summary median sampling date
meddate<-sp_dat_list %>% 
  bind_rows() %>%
  group_by(HAUL_ID) %>%
  summarize(HAUL_DAY = mean(YDAY)) %>%
  summarize(med=median(HAUL_DAY),mean=mean(HAUL_DAY)) #144 for GOA, 134 EBS



## plot species length data maps for every species
#for(mysp in myspp){
# ggsave(paste0("Figures/SpeciesLengthDataMaps/",myregion,"_",mysp,".png"),
#       plotLengthMaps(sp_dat_list[[mysp]],myregion),
#       width=12, height=8, units="in")
#}


# Fit models: 
# 1) Base model: ~ 1 + s(YDAYSC,k=5) + MESH
# 2) Year factor to look at interannual change: ~ 0 + s(YDAYSC, k=5) + YEARF + MESH 
# 4) Temperature effect: ~ 0 + s(YDAYSC, k=5) + MESH + tempvar
# 6) Year linear effect: ~ 0 + s(YDAYSC, k=5) + YEARSC + MESH

# Use 40km cutoff for EBS mesh, and 30 for GOA 
 

######################
### MODEL FITTING FUNCTIONS
######################

######################
### Base Model - No Temp/Year effects
######################

fit_base_mod_func<-function(sp,foldername,myregion, myfamily){
  
  #Check for folder to save model fits, create if not existing  
  today<-Sys.Date()
  mypath<-here("Results","ModelResults",paste0(foldername,"_",today)) 
  if(!dir.exists(mypath)) dir.create(mypath)
  
  mesh_cutoff<-switch(myregion,
                      goa=30,
                      ebs=40)
  
  sppdat<-sp_dat_list[[sp]]
  mesh <- make_mesh(sppdat, xy_cols = c("X", "Y"), cutoff = mesh_cutoff)
  
  fit1 <- try(
    sdmTMB(
      CORRECTED_LENGTH ~ 1 + s(YDAYSC,k=5) + MESH, #+ (1|YEARF), # 
      data = sppdat,
      mesh = mesh,
      spatial = "off",
      spatiotemporal = "iid",
      time="YEAR",
      anisotropy = FALSE,
      family=myfamily
  )
)

  if(inherits(fit1, "try-error")){
      print(paste("model 1 error for species",sp))
      reslist<-list("error","error","error")
      
  } else { if (sanity(fit1)$all_ok ==F) fit1<-try(run_extra_optimization(fit1))}
  if(inherits(fit1, "try-error")){
    print(paste("model 1 error for species",sp))
    reslist<-list("error","error","error")} else {
      
      COEFS<-tidy(fit1,conf.int=TRUE)
      COEFS$SPECIES<-sp
      res<-residuals(fit1)
      
      vr<-visreg(fit1,xvar="YDAYSC",plot=FALSE)
      partials<-vr$res
      partials$SPECIES<-sp
      fit<-vr$fit
      fit$SPECIES<-sp
   
      # Save model object separately for each species
      # Compile summary stats into an object to save
      # Function will output visreg, AIC, summary stats, but not model object.
      modreslist<-list("MODEL"=fit1,"COEFS"=COEFS,"AIC"=AIC(fit1),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"MESH_CUTOFF"=mesh_cutoff)
      saveRDS(modreslist,file=here("Results","ModelResults",paste0(foldername,"_",today),paste0(sp,".rds")))
      reslist<-list("COEFS"=COEFS,"AIC"=AIC(fit1),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"MESH_CUTOFF"=mesh_cutoff)
    }
  return(reslist)
}

######################
### Year Factor Model 
######################


fit_YearF_mod_func<-function(sp,foldername,myregion,myfamily){

#Check for folder to save model fits, create if not existing  
today<-Sys.Date()
mypath<-here("Results","ModelResults",paste0(foldername,"_",today)) 
if(!dir.exists(mypath)) dir.create(mypath)

sppdat<-sp_dat_list[[sp]]

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

    vr<-visreg(fit2,xvar="YEARF",cond=list(YDAYSC=0),plot=FALSE,scale="response")
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

####
####

fit_temp_mod_func<-function(sp,foldername,myregion,tempvar,myfamily){

  #Check for folder to save model fits, create if not existing  
  today<-Sys.Date()
  mypath<-here("Results","ModelResults",paste0(foldername,"_",today)) 
  if(!dir.exists(mypath)) dir.create(mypath)

  mesh_cutoff<-switch(myregion,
                      goa=30,
                      ebs=40)
    
  sppdat<-sp_dat_list[[sp]]
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
  
} else { if (sanity(fit4)$all_ok ==F) fit4<-try(run_extra_optimization(fit4))}
if(inherits(fit4, "try-error")){
  print(paste("model 4 error for species",sp))
  reslist<-list("error","error","error")} else {

COEFS<-tidy(fit4,conf.int=TRUE)
COEFS$SPECIES<-sp
res<-residuals(fit4)

vr<-visreg(fit4,xvar=tempvar,cond=list(YDAYSC=0),plot=FALSE,scale="response")
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

####
####


fit_year_mod_func<-function(sp,foldername,myregion,myfamily){
  
  #Check for folder to save model fits, create if not existing  
  today<-Sys.Date()
  mypath<-here("Results","ModelResults",paste0(foldername,"_",today)) 
  if(!dir.exists(mypath)) dir.create(mypath)
  
  mesh_cutoff<-switch(myregion,
                      goa=30,
                      ebs=40)
  
  sppdat<-sp_dat_list[[sp]]
  mesh <- make_mesh(sppdat, xy_cols = c("X", "Y"), cutoff = mesh_cutoff)

fit6 <- try(sdmTMB(
  CORRECTED_LENGTH ~  0 + s(YDAYSC, k=5) + YEARSC + MESH, # 
  data = sppdat,
  mesh = mesh,
  spatial = "off",
  spatiotemporal = "iid",
  time="YEAR",
  family=myfamily
  )
)
if(inherits(fit6, "try-error")){
  print(paste("model 6 error for species",sp))
  reslist<-list("error","error","error")
} else {if (sanity(fit6)$all_ok ==F) fit6<-try(run_extra_optimization(fit6))}
if(inherits(fit6, "try-error")){
  print(paste("model 6 error for species",sp))
  reslist<-list("error","error","error")} else {

COEFS<-tidy(fit6,conf.int=TRUE)
COEFS$SPECIES<-sp
res<-residuals(fit6)

vr<-visreg(fit6,xvar="YEARSC",cond=list(YDAYSC=0),plot=FALSE)
partials<-vr$res
partials$SPECIES<-sp
fit<-vr$fit
fit$SPECIES<-sp

# Save model object separately for each species
# Compile summary stats into an object to save
# Function will output visreg, AIC, summary stats, but not model object.
modreslist<-list("MODEL"=fit6,"COEFS"=COEFS,"AIC"=AIC(fit6),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"MESH_CUTOFF"=mesh_cutoff)
saveRDS(modreslist,file=here("Results","ModelResults",paste0(foldername,"_",today),paste0(sp,".rds")))
reslist<-list("COEFS"=COEFS,"AIC"=AIC(fit6),"PREDS"=NA,"RES"=res,"PARTIALS"=partials,"FIT"=fit,"MESH_CUTOFF"=mesh_cutoff)
}
return(reslist)

}
######################
######################
### RUN MODELS
######################
######################

# Select error family
myfamily<-gaussian(link = "identity"); famname<-"gauss"
#myfamily<-gengamma(link="log"); famname<- "genGam"
#myfamily<-gaussian(link="log"); famname<- "logN"

myfoldername<-paste0(myregion,"_BaseModelFits_k5_noSP_",famname)
mybasefoldername<-myfoldername
mod_base_lists <- lapply(myspp, fit_base_mod_func,foldername=myfoldername,myregion=myregion,myfamily=myfamily)
names(mod_base_lists)<-myspp#[sp_tofit]
saveRDS(mod_base_lists,file=here("Results","ModelResults",paste0(myfoldername,"_",today()),"AllSppLists.rds"))

myfoldername=paste0(myregion,"_YearFModelFits_k5_noSP_",famname)
myYearFfoldername<-myfoldername
mod_YearF_lists <- lapply(myspp, fit_YearF_mod_func,foldername=myfoldername,myregion=myregion,myfamily=myfamily)
names(mod_YearF_lists)<-myspp#[sp_tofit]
saveRDS(mod_YearF_lists,file=here("Results","ModelResults",paste0(myfoldername,"_",today()),"AllSppLists.rds"))

tempvar="JFMAM"
myfoldername=paste0(myregion,"_TempModelFits_k5_noSP_",famname,"_",tempvar)
myTempfoldername<-myfoldername
mod_temp_lists <- lapply(myspp, fit_temp_mod_func,foldername=myfoldername,myregion=myregion,tempvar=tempvar,myfamily=myfamily)
names(mod_temp_lists)<-myspp#[sp_tofit]
saveRDS(mod_temp_lists,file=here("Results","ModelResults",paste0(myfoldername,"_",today()),"AllSppLists.rds"))

myfoldername=paste0(myregion,"_LinearYearModelFits_k5_noSP_",famname)
myYearfoldername<-myfoldername
mod_year_lists <- lapply(myspp, fit_year_mod_func,foldername=myfoldername,myregion=myregion,myfamily=myfamily)
names(mod_year_lists)<-myspp#[sp_tofit]
saveRDS(mod_year_lists,file=here("Results","ModelResults",paste0(myfoldername,"_",today()),"AllSppLists.rds"))


######################
######################
### PROCESS / PLOT OUTPUT (mostly done in other files, code kept below for reference)
######################
######################
date<-"2024-12-12"
#myfoldername="goa_BaseModelFits_noSpatial"
mod_base_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
#myfoldername="ebs_YearFModelFits_11spp_to2021"
mod_YearF_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
#myfoldername="goa_TempModelFits_27spp_JFMAM"
mod_temp_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))
#myfoldername="goa_LinearYearModelFits_27spp"
mod_year_lists<-readRDS(here("Results","ModelResults",paste0(myfoldername,"_",date),"AllSppLists.rds"))


################################
## PLOTTING YEAR FACTOR MODELS
## 
################################
################################
#need to exclude lists with errors
#then merge into one dataframe for ggplot
inds<-which(mod_YearF_lists %>% map(1)!="error")
allcoefs_YearF <- mod_YearF_lists[inds] %>% map(1) %>% dplyr::bind_rows()
#allpreds <- mod_YearF_lists[inds] %>% map(3) %>% dplyr::bind_rows()
allpartials_YearF <- mod_YearF_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_YearF <- mod_YearF_lists[inds] %>% map(6) %>% dplyr::bind_rows()


#ggplot(allpreds %>% filter(YEAR==2021, YDAY==145),aes(YEARF,est)) +
#  geom_point() +
  #geom_linerange()
#  facet_wrap(~SPECIES)

#ggplot(allpreds %>% filter(YEAR==2021, YEARF==2021), aes(YDAY, (est),
#            ymin = (est - 1.96 * est_se), ymax = (est + 1.96 * est_se))) +
#  geom_line() + geom_ribbon(alpha = 0.4) +
#  facet_wrap(~SPECIES)

#create species/common names d.f., and determine y plot location for printing
spdf<-data.frame("SPECIES"=unique(allfits_YearF$SPECIES))
spdf<-left_join(spdf,goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],by=join_by("SPECIES"=="SPECIES_NAME"),multiple="any")
#yloc<-allpartials_YearF %>% filter(!is.infinite(visregRes)) %>% group_by(SPECIES)  %>% summarize(yloc=max(visregRes)) %>% select(SPECIES,yloc)
#spdf<-full_join(spdf,yloc)
altmax<-allfits_YearF %>% filter(YDAYSC==0) %>% group_by(SPECIES) %>% summarize(maxest= max(visregUpr))
spdf<-full_join(spdf,altmax) %>% group_by(SPECIES) %>% mutate(maxy=maxest+1)

#convert year to numeric for plotting
allfits_YearF <-allfits_YearF %>%
  mutate(YEAR=as.numeric(levels(YEARF))[YEARF]) 

minY<-min(allfits_YearF$YEAR)

ggplot(allfits_YearF,aes(YEAR,visregFit)) +
  geom_point() +
  geom_errorbar(aes(x=YEAR,ymin=visregLwr,ymax=visregUpr)) +
  facet_wrap(~SPECIES,scales="free_y") +
  geom_text(data=spdf,aes(label=COMMON_NAME,y=maxy),x=minY+1,hjust=0,size=3,color="gray30") +
  theme(strip.text = element_text(face = "italic")) +
  ylab("Larval Length (mm)") +
  xlab("Year")


if(myregion=="ebs"){
  ggsave(file=here("Figures",paste0(myYearFfoldername,"_",today,".png")),width=9,height=6)
}else{ggsave(file=here("Figures",paste0(myYearFfoldername,"_",today,".png")),width=14,height=10)}
# 
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
## 
################################
################################

# color palette
col1<-'#F87505'
col2<-'#048F96'

## need to exclude species with errors
inds<-which(mod_temp_lists %>% map(3)!="error")
allcoefs_temp <- mod_temp_lists[inds] %>% map(1) %>% dplyr::bind_rows()
allpreds_temp <- mod_temp_lists[inds] %>% map(3) %>% dplyr::bind_rows()
allpartials_temp <- mod_temp_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_temp <- mod_temp_lists[inds] %>% map(6) %>% dplyr::bind_rows()
allmesh_temp <- mod_temp_lists[inds] %>% map(8) %>% dplyr::bind_rows()

#tempvar="JFMAM"
TEMPeffects<-allcoefs_temp %>% filter(term==tempvar) %>% arrange(SPECIES)

histplot<-TEMPeffects %>%
  mutate(ispos = estimate>0, 
         signif = sign(conf.low)==sign(conf.high)) %>%
  mutate(slope = case_when(ispos == TRUE ~ 'positive',
                           ispos == FALSE ~ 'negative')) %>%
  mutate(significance = case_when(signif ==TRUE ~ 'p < 0.05',
                                  signif ==FALSE ~ 'NS'))

toplot<- histplot %>%  select(SPECIES,slope,significance) %>%
full_join(allfits_temp)

#create species/common names d.f.
spdf<-data.frame("SPECIES"=unique(toplot$SPECIES))
spdf<-left_join(spdf,goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],by=join_by("SPECIES"=="SPECIES_NAME"),multiple="any")
yloc<-allpartials_temp %>% filter(!is.infinite(visregRes)) %>% group_by(SPECIES)  %>% summarize(yloc=max(visregRes)) %>% select(SPECIES,yloc)
spdf<-full_join(spdf,yloc)

#Find max y for each species for labels
altmax<-toplot %>% filter(YDAYSC==0) %>% group_by(SPECIES) %>% summarize(maxest= max(visregUpr))
spdf<-full_join(spdf,altmax) %>% group_by(SPECIES) %>% mutate(maxy=max(yloc,maxest))

#Find x location for label (minimum temp on x scale)
xloc<-min(allfits_temp[,tempvar])

ggplot() +
  geom_point(data=allpartials_temp,aes(x=.data[[tempvar]],y=visregRes),alpha=0.2,size=0.2,color="gray") +
  geom_ribbon(data=toplot %>% filter(YDAYSC==0),
              aes(.data[[tempvar]], visregFit, ymin = visregLwr,
                  ymax = visregUpr, fill=slope, alpha=significance)) +
  geom_line(data=toplot %>% filter(YDAYSC==0), aes(.data[[tempvar]], visregFit,linetype=significance)) +
  facet_wrap(~SPECIES,scales="free_y") +
  geom_text(data=spdf,aes(label=COMMON_NAME,y=maxy),x=xloc,hjust=0,size=3,color="gray30") +
  theme(strip.text = element_text(face = "italic")) +
  xlab(paste0(tempvar," Temperature (C)")) +
  ylab("Larval Length (mm, May 25)") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
  scale_linetype_manual(values=c('NS' = "dashed", "p < 0.05" = "solid")) +
  scale_fill_manual(values = c("negative" = col1, "positive" = col2))

if(myregion=="ebs"){
ggsave(file=here("Figures",paste0(myTempfoldername,"_",today,".png")),width=9,height=6)
}else{ggsave(file=here("Figures",paste0(myTempfoldername,"_",today,".png")),width=14,height=10)}

# color this histogram to match effects plots
ggplot(histplot,aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
  geom_histogram(binwidth=0.5,center=0.25,show.legend=FALSE)+
  xlab("Estimated temp. effect (mm/deg C)") +
  geom_vline(xintercept=0,linetype="dashed") +
#  scale_y_continuous(breaks=c(0,1,2,3),expand=c(0,0))+
#  scale_x_continuous(limits=c(-2.5,2.5))+
  theme_classic() +
  scale_colour_manual(values = c("negative" = col1, "positive" = col2))+
  scale_fill_manual(values = c("negative" = col1, "positive" = col2))

ggsave(file=here("Figures",paste0(myTempfoldername,"_hist_",today,".png")),width=3,height=2.5)


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

####################################
## PLOTTING LINEAR YEAR MODELS
####################################
####################################
####################################
## need to exclude species with errors 
inds<-which(mod_year_lists %>% map(3)!="error")
allcoefs_year <- mod_year_lists[inds] %>% map(1) %>% dplyr::bind_rows()
#allpreds_year <- mod_year_lists[inds] %>% map(3) %>% dplyr::bind_rows()
allpartials_year <- mod_year_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_year <- mod_year_lists[inds] %>% map(6) %>% dplyr::bind_rows() 

YEAReffects<-allcoefs_year %>%
  filter(term=="YEARSC") %>%  
  arrange(SPECIES)  

histplot<-YEAReffects %>%
  mutate(ispos = estimate>0, 
         signif = sign(conf.low)==sign(conf.high)) %>%
  mutate(slope = case_when(ispos == TRUE ~ 'positive',
                           ispos == FALSE ~ 'negative')) %>%
  mutate(significance = case_when(signif ==TRUE ~ 'p < 0.05',
                                  signif ==FALSE ~ 'NS'))
  
toplot<- histplot %>%  select(SPECIES,slope,significance) %>% full_join(allfits_year)

#create species/common names d.f.
spdf<-data.frame("SPECIES"=unique(toplot$SPECIES))
spdf<-left_join(spdf,goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],by=join_by("SPECIES"=="SPECIES_NAME"),multiple="any")
yloc<-allpartials_year %>% filter(!is.infinite(visregRes)) %>% group_by(SPECIES)  %>% summarize(yloc=max(visregRes)) %>% select(SPECIES,yloc)
spdf<-full_join(spdf,yloc)

altmax<-toplot %>% filter(YDAYSC==0) %>% group_by(SPECIES) %>% summarize(maxest= max(visregUpr))
spdf<-full_join(spdf,altmax) %>% group_by(SPECIES) %>% mutate(maxy=max(yloc,maxest))

minY<-min(allfits_year$YEARSC)+2000

ggplot() +
  geom_point(data=allpartials_year,aes(x=YEARSC+2000,y=visregRes),
             alpha=0.2,size=0.2,color="gray") +
  geom_ribbon(data=toplot %>% filter(YDAYSC==0), 
              aes(YEARSC+2000, (visregFit), ymin = visregLwr, 
                  ymax = visregUpr, fill=slope, alpha=significance)) +
  geom_line(data=toplot %>% filter(YDAYSC==0), 
             aes(YEARSC+2000, visregFit,color=slope, linetype=significance)) +
  facet_wrap(~SPECIES,scales="free_y") +
  geom_text(data=spdf,aes(label=COMMON_NAME,y=maxy),x=minY+1,hjust=0,size=3,color="gray30") +
    xlab("Year") +
  theme(strip.text = element_text(face = "italic")) +
  ylab("Larval Length (mm, May 25)") +
  scale_alpha_discrete(range=c(0.1, 0.6)) +
#  scale_fill_discrete(type=scale_fill_brewer) +
  scale_linetype_manual(values=c('NS' = "dashed", "p < 0.05" = "solid"))

#ggsave(file=here("Figures",paste0("GOAYear_30spp_",today,".png")),width=14,height=10)
#ggsave(file=here("Figures",paste0("EBSYear_11spp_",today,".png")),width=9,height=6)
if(myregion=="ebs"){
  ggsave(file=here("Figures",paste0(myYearfoldername,"_",today,".png")),width=9,height=6)
}else{ggsave(file=here("Figures",paste0(myYearfoldername,"_",today,".png")),width=14,height=10)}


#color this histogram to match effects plots
ggplot(histplot,aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
  geom_histogram(binwidth=0.02,center=0.01, show.legend=FALSE)+
 xlab("Estimated year effect (mm/year)") +
 geom_vline(xintercept=0,linetype="dashed") +
#  scale_y_continuous(breaks=c(3,6,9),expand=c(0,0))+
#  scale_x_continuous(breaks=c(-0.06,0,0.06))+
  theme_classic() 

ggsave(file=here("Figures",paste0(myYearfoldername,"_hist_",today,".png")),width=3,height=2.5)

## scaled for EBS
## 
ggplot(histplot,aes(x=estimate,fill=slope,color=slope,alpha=significance)) +
  geom_histogram(binwidth=0.1,center=0.05, show.legend=FALSE)+
  xlab("Estimated year effect (mm/year)") +
  geom_vline(xintercept=0,linetype="dashed") +
  scale_y_continuous(breaks=c(2,4,6),expand=c(0,0))+
  scale_x_continuous(breaks=c(-0.3,0,0.3),limits = c(-0.4,0.4))+
  theme_classic() 
  
ggsave(file=here("Figures",paste0(myYearfoldername,"_hist_",today,".png")),width=3,height=2.5)


######################
######################
##  Plotting YDAY effect from Base model
######################
######################

## need to exclude species with errors 
inds<-which(mod_base_lists %>% map(3)!="error")
allpartials_base <- mod_base_lists[inds] %>% map(5) %>% dplyr::bind_rows() 
allfits_base <- mod_base_lists[inds] %>% map(6) %>% dplyr::bind_rows() 
allfits_base$YDAY <- allfits_base$YDAYSC +145
allpartials_base$YDAY <- allpartials_base$YDAYSC +145

toplot<- allfits_base

#create species/common names d.f.
spdf<-data.frame("SPECIES"=unique(toplot$SPECIES))
spdf<-left_join(spdf,goa_ebs[,c("SPECIES_NAME","COMMON_NAME")],by=join_by("SPECIES"=="SPECIES_NAME"),multiple="any")
yloc<-allpartials_base %>% filter(!is.infinite(visregRes)) %>% group_by(SPECIES)  %>% summarize(yloc=max(visregRes)) %>% select(SPECIES,yloc)
spdf<-full_join(spdf,yloc)

altmax<-toplot %>% group_by(SPECIES) %>% summarize(maxest= max(visregUpr))
spdf<-full_join(spdf,altmax) %>% group_by(SPECIES) %>% mutate(maxy=max(yloc,maxest))

ggplot() +
  geom_point(data=allpartials_base,aes(x=YDAY,y=visregRes),
             alpha=0.2,size=0.2,color="gray") +
  geom_ribbon(data=toplot, 
              aes(YDAY, (visregFit), ymin = visregLwr, 
                  ymax = visregUpr),fill="steelblue2") +
  geom_line(data=toplot,  
            aes(YDAY, visregFit)) +
  facet_wrap(~SPECIES,scales="free_y") +
  geom_text(data=spdf,aes(label=COMMON_NAME,y=maxy),x=1980,hjust=0,size=3,color="gray30") +
  xlab("Day of Year") +
  theme(strip.text = element_text(face = "italic")) +
  ylab("Larval Length (mm, May 25)") 

#ggsave(file=here("Figures",paste0("GOAYear_30spp_",today,".png")),width=14,height=10)
#ggsave(file=here("Figures",paste0("EBSYear_11spp_",today,".png")),width=9,height=6)
if(myregion=="ebs"){
  ggsave(file=here("Figures",paste0(mybasefoldername,"_",today,".png")),width=9,height=6)
}else{ggsave(file=here("Figures",paste0(mybasefoldername,"_",today,".png")),width=14,height=10)}




