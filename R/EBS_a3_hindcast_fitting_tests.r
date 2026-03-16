# Run the ACLIM 2.0 hindcast (1991–2021) for preliminary use with model fitting
library(Rpath)

#source("fits/clim_sim_hcr_obs_err_2_fits.r")
#source("fits/bioen_projections_fits.r")
#setwd("C:/Users/kerim.aydin/Work/src/sean_rpath/aclim2_rpath")
# Helper functions preceded with TMP for eventual package inclusion

source("tmp_ecofitting.R")

# need to have repo permissions to use these files directly
mdir <- "https://raw.githubusercontent.com/aw2113/aclim3_rpath/refs/heads/main/"

# Setup Base Ecopath 
  Ebase <- "model/ebs_aclim3_74bio_base.csv"  # Base biomass, production, fishing, etc.
  Ediet <- "model/ebs_aclim3_74bio_diet.csv"  # Diet matrix
  Eped  <- "model/ebs_aclim3_74bio_pedigree.csv"  # Data pedigree = quality of input data
  Estz  <- "model/ebs_aclim3_74bio_stanzas.csv"  # Stanzas
  Estg  <- "model/ebs_aclim3_74bio_stanza_groups.csv" # Stanza groups

# FITTING DATA FILES
  datfiles <- list(
    catchfile           = "data/ebs_aclim3_74bio_catch_long.csv",
    surveyfile_shelf    = "data/ebs_aclim3_74bio_fitting_survey_index_Q1.csv",
    surveyfile_slope    = "data/EBS_slope_jan2025_fitting_index_aclim3_Q1.csv",
    surveyfile_mammals  = "data/ebs_aclim3_74bio_mammal_index.csv",
    surveyfile_birds    = "data/seabird_colony_data.csv")  

# Download if needed - set to false by default anticipating call from App
UPDATE_FILES <- FALSE
if (UPDATE_FILES){
  download.file(paste(mdir,Ebase,sep=''), Ebase)
  download.file(paste(mdir,Ediet,sep=''), Ediet)
  download.file(paste(mdir,Eped,sep=''), Eped)  
  download.file(paste(mdir,Estz,sep=''), Estz)
  download.file(paste(mdir,Estg,sep=''), Estg)
  for(i in datfiles){download.file(paste(mdir,i,sep=''),i)} 
}


  bio_dat <- rbind(read.csv(datfiles$surveyfile_shelf),read.csv(datfiles$surveyfile_slope),
                   read.csv(datfiles$surveyfile_mammals),read.csv(datfiles$surveyfile_birds))
  
  unbal <- rpath.stanzas(read.rpath.params(Ebase, Ediet, Eped, Estg, Estz)) # unbalanced
  bal   <- rpath(unbal) # balanced

# Setup years and base ecosim scenario
  hind_years <- 1991:2021
  scene0 <- rsim.scenario(bal, unbal, years = hind_years)
  # check if in equilibrium
  # run.hind <- rsim.run(scene0, method='AB', years=hind_years)
  # rsim.plot(run.hind)

  
# Read in biomass to fit 
  scene1 <- TMP.read.fitting.biomass(scene0, bio_dat)
  #scene1 <- read.fitting.biomass(scene0, "fits/EBS_2023_fitting_survey_index_2_Q1.csv")

# Read in catch to fit.  This does not on its own apply read-in catch to catch forcing.  
  scene2 <- read.fitting.catch(scene1, datfiles$catchfile)

# This function copies the fitting catch table into the forced_catch matrix 
# This does not zero out previous Effort-applied (equilibrium) catch
  scene3 <- fitcatch.to.forcecatch(scene2)

# Zeroing out effort matrix manually (this will turn off detritus flows that are by gear)
  scene3$fishing$ForcedEffort[] <- 0  
  
  
scene_fit <- scene3


#scene_fit$fitting$Biomass$Type <- "index"
run0 <- rsim.run(scene_fit, method='AB', years=hind_years)
rsim.plot(run0)
  
rsim.fit.table(scene_fit,run0)
rsim.fit.obj(scene_fit,run0,verbose=T)

scene3$params$
obj_function <- function(mzero, scene, years){
  scene$params$MzeroMort[names(mzero)] <- mzero
  frun <- rsim.run(scene, method='AB', years=years)
  return(rsim.fit.obj(scene,frun,F))
}


fit_species <- c("arrowtooth_adu","shallow_demersals") 
mzero_vector <- scene3$params$MzeroMort[fit_species]

mzero_vector["arrowtooth_adu"] <- 0.005
obj_function(mzero_vector, scene3, hind_years)
obj_function(NULL, scene3, hind_years)

fit_species <- c("shallow_demersals") 
mzero_vector <- scene3$params$MzeroMort[fit_species]

result <- optim(mzero_vector, obj_function, scene=scene3, years=hind_years) 
# First Base run ends here

obj_function(NULL, scene3, hind_years)
obj_function(result$par, scene3, hind_years)

scene4 <- adjust.scenario(scene3,"MzeroMort",group=names(result$par),value=result$par)
obj_function(NULL, scene4, hind_years)

f1 <- rsim.run(scene3, method='AB', years=hind_years)
f2 <- rsim.run(scene4, method='AB', years=hind_years)

delta <- rsim.fit.table(scene4,f2) - rsim.fit.table(scene3,f1)

par(mfrow=c(1,2))
TMP.rsim.plot.full(scene3,f1,"shallow_demersals")
TMP.rsim.plot.full(scene4,f2,"shallow_demersals")

TMP.rsim.plot.full(scene3,f1,"skates")
TMP.rsim.plot.full(scene4,f2,"skates")

scene3$fitting$Biomass$wt[!(scene3$fitting$Biomass$Group %in% fit_species)] <- 0
scene3$fitting$Biomass$Type[!(scene3$fitting$Biomass$Group %in% fit_species)] <- "absolute"

obj_function(NULL, scene3, hind_years)
result <- optim(mzero_vector, obj_function, scene=scene3, years=hind_years,method="Brent",lower=-10,upper=10) 

scene4 <- adjust.scenario(scene3,"MzeroMort",group=fit_species,value=result$par)

obj_function(NULL, scene4, hind_years)
TMP.rsim.plot.full(scene3,f1,"shallow_demersals")
TMP.rsim.plot.full(scene4,f2,"shallow_demersals")
delta <- rsim.fit.table(scene4,f2) - rsim.fit.table(scene3,f1)

frun <- rsim.run(scene3, method='AB', years=hind_years)
rsim.fit.obj(scene3,frun,T)
TMP.rsim.plot.full(scene4,f2,"shallow_demersals")




#rsim.fit.obj(scene_fit, run0,FALSE)

rsim.fit.obj(scene_fit, run0,TRUE)$Biomass

par(mfrow=c(1,2))
rsim.plot.biomass(scene_fit, run0, "pollock_adu")
rsim.plot.catch(scene_fit, run0, "pollock_adu")

slist <- scene_fit$params$spname[scene_fit$params$FishFrom+1]
sum((scene_fit$params$FishQ*
scene_fit$params$B_BaseRef[slist])[slist=="motile_epifauna"])



data.frame(scene_fit$params$spname[scene_fit$params$FishFrom+1], scene_fit$params$FishQ)


## Sense tests
# set-up sense

scene <- rsim.scenario(bal, unbal, years = 1900:1999)

NUM_RUNS <- 500
parlist <- as.list(rep(NA,NUM_RUNS)) # lists to store generated parameter sets
kept <- rep(NA,NUM_RUNS) # vector to track which parameter sets survived burn-in
set.seed(666) # so results can be replicated

ptm <- proc.time()
for (i in 1:NUM_RUNS){
  EBSsense <- scene
  # INSERT SENSE ROUTINE BELOW
  parlist[[i]]<- scene$params # Base ecosim params
  parlist[[i]]<- rsim.sense(scene, unbal, Vvary = c(-4.5,4.5), Dvary = c(0,0)) # Replace the base params with Ecosense params
  # EBSsense$start_state$Biomass <- parlist[[i]]$B_BaseRef
  parlist[[i]]$BURN_YEARS <- 100 # Set Burn Years in parlist
  EBSsense$params <- parlist[[i]]
  EBStest <- rsim.run(EBSsense, method="AB", years=1900:1999)
  failList <- which(is.na(EBStest$end_state$Biomass))
  {if (length(failList)>0)
  {cat(i,": fail in year ",EBStest$crash_year,": ",EBStest$params$spname[failList],
       "\n"); kept[i]<-F; parlist[[i]] <- NULL; flush.console()}
    else
    {cat(i,": success!\n"); kept[i]<-T; parlist[[i]]$BURN_YEARS <- 1; flush.console()}}
  # parlist[[i]]$BURN_YEARS <- 1
}
proc.time() - ptm

# Setup random fishing
par(mfrow=c(1,2))
fish_test  <- runif(100,.01,.04)

# below assumes sense loop was run earlier, so scene, parlist, kept exist
# run with base parameters
scene_test <- scene
scene_test$fishing$ForcedEffort[] <- 0
scene_test$fishing$ForcedCatch[,"skates"] <- fish_test  #"pollock_adu" 
run1<- rsim.run(scene_test, method="AB", 1900:1999)
plot(1900:1999, fish_test, ylim=c(0,.04))
lines(1900:1999, run1$annual_Catch[,"skates"])

# run with random parameters
scene_test <- scene
scene_test$fishing$ForcedEffort[] <- 0
scene_test$fishing$ForcedCatch[,"skates"] <- fish_test #"pollock_adu"
scene_test$params <- parlist[[which(kept==T)[2]]]
run2<- rsim.run(scene_test, method="AB", 1900:1999)
plot(1900:1999, fish_test, ylim=c(0,.04))
lines(1900:1999, run2$annual_Catch[,"skates"])











# write.csv(run.hind$annual_Biomass,"Kerim_biomass_testDec2023.csv",row.names=F)
