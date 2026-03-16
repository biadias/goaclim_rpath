# Run the ACLIM 2.0 hindcast (1991–2021) for preliminary use with model fitting
library(Rpath)

# TMP functions that will move to the Rpath package when finished
# are in the following files
fup <- function(){
  source("R/merge_ecofitting.R")
  source("R/ecofitting_plots.R")
}

fup()

# Base Ecopath Files 
  Ebase <- "ebs/model/ebs_aclim3_76bio_base.csv"  # Base biomass, production, fishing, etc.
  Ediet <- "ebs/model/ebs_aclim3_76bio_diet.csv"  # Diet matrix
  Eped  <- "ebs/model/ebs_aclim3_76bio_pedigree.csv"  # Data pedigree = quality of input data
  Estz  <- "ebs/model/ebs_aclim3_76bio_stanzas.csv"  # Stanzas
  Estg  <- "ebs/model/ebs_aclim3_76bio_stanza_groups.csv" # Stanza groups

# FITTING DATA FILES
  datfiles <- list(
    catchfile           = "ebs/data/ebs_aclim3_76bio_catch_long.csv",
    surveyfile_shelf    = "ebs/data/ebs_shelf_fitting_index_aclim3_Q1.csv",
    surveyfile_slope    = "ebs/data/ebs_slope_fitting_index_aclim3_Q1.csv",
    surveyfile_mammals  = "ebs/data/ebs_aclim3_mammal_index.csv",
    surveyfile_birds    = "ebs/data/seabird_colony_data.csv",
    hindcast_ocean      = "ebs/climate/a2_hind.csv"
    )  

################################################################################  
# Download files if needed - set to false by default (doesn't update from repo)
# need to have repo permissions to use these files directly.  Assumes repo has
# the model/ and data/ subdirectories with files above.

#mdir <- "https://raw.githubusercontent.com/aw2113/aclim3_rpath/refs/heads/main/"  
mdir <- "C:/Users/kerim.aydin/Work/src/aclim3_rpath/"
UPDATE_FILES <- FALSE
if (UPDATE_FILES){ #Use download.file instead of file.copy for remote
  file.copy(paste(mdir,Ebase,sep=''), Ebase, overwrite=T) 
  file.copy(paste(mdir,Ediet,sep=''), Ediet, overwrite=T)
  file.copy(paste(mdir,Eped,sep=''), Eped, overwrite=T)  
  file.copy(paste(mdir,Estz,sep=''), Estz, overwrite=T)
  file.copy(paste(mdir,Estg,sep=''), Estg, overwrite=T)
  for(i in datfiles){file.copy(paste(mdir,i,sep=''),i,overwrite=T)}
}
################################################################################

# Assemble biomass files into a single data table.  All biomass files must have
# the same columns (even colums unused by fitting must match)
  bio_dat <- rbind(read.csv(datfiles$surveyfile_shelf),read.csv(datfiles$surveyfile_slope),
                   read.csv(datfiles$surveyfile_mammals),read.csv(datfiles$surveyfile_birds))

  catch_dat <- read.csv(datfiles$catchfile)
  
  hind_dat  <- read.csv(datfiles$hindcast_ocean)

# Basic Ecopath balance
  unbal <- rpath.stanzas(read.rpath.params(Ebase, Ediet, Eped, Estg, Estz)) # unbalanced
  bal   <- rpath(unbal) # balanced

# Set hindcast years  
  hind_years <- 1991:2021

# Setup base scenario
  scene0 <- rsim.scenario(bal, unbal, years = hind_years)  
  
# Read in biomass to fit 
  scene1 <- read.fitting.biomass(scene0, bio_dat)

# Read in catch to fit.  This does not on its own apply read-in catch to catch forcing.  
  scene2 <- read.fitting.catch(scene1, catch_dat)

# fitcatch.to.forcecatch (1) turns off gear Effort, (2) moves Ecopath baseline 
# fishing to Frate forcing for all years and (3) adds forcecatch for years where
# catch is provided, 0ing the baseline Frate in those years.  
  scene3 <- fitcatch.to.forcecatch(scene2, bal)
# Currently forced catch and forced F do not go to detritus, so need to keep
# that group constant.
  scene3$forcing$ForcedBio[,"discards_offal"] <- bal$Biomass["discards_offal"]  

# Climate forcing of lower trophics
  start_clim <- which(hind_dat$year==min(hind_years) & hind_dat$mo==1)
  end_clim   <- which(hind_dat$year==max(hind_years) & hind_dat$mo==12)
  last_clim  <- (1+max(hind_years)-min(hind_years)) * 12 
  hind_clim  <- hind_dat[start_clim:end_clim,]
  scene4 <- scene3  
  BaseB  <- scene4$params$B_BaseRef
  scene4$forcing$ForcedBio[1:last_clim,"euphausiids"]      <- BaseB["euphausiids"]      * hind_clim[,"eup"]
  scene4$forcing$ForcedBio[1:last_clim,"copepods"]         <- BaseB["copepods"]         * hind_clim[,"cop"]
  scene4$forcing$ForcedBio[1:last_clim,"pelagic_microbes"] <- BaseB["pelagic_microbes"] * hind_clim[,"mzl"]
  scene4$forcing$ForcedBio[1:last_clim,"benthic_microbes"] <- BaseB["benthic_microbes"] * hind_clim[,"mzl"]
  scene4$forcing$ForcedBio[1:last_clim,"lg_phytoplankton"] <- BaseB["lg_phytoplankton"] * hind_clim[,"phl"]
  scene4$forcing$ForcedBio[1:last_clim,"sm_phytoplankton"] <- BaseB["sm_phytoplankton"] * hind_clim[,"phs"]      
  # CHECK WITH ANDY WHY OFFSET 1:371 versus 2:372
  # scene$forcing$ForcedBio[1:371,"euphausiids"]      <- BaseB["euphausiids"]      * hind_clim[2:372,"eup"]
  # scene$forcing$ForcedBio[1:371,"copepods"]         <- BaseB["copepods"]         * hind_clim[2:372,"cop"]
  # scene$forcing$ForcedBio[1:371,"pelagic_microbes"] <- BaseB["pelagic_microbes"] * hind_clim[2:372,"mzl"]
  # scene$forcing$ForcedBio[1:371,"benthic_microbes"] <- BaseB["benthic_microbes"] * hind_clim[2:372,"mzl"]
  # scene$forcing$ForcedBio[1:371,"lg_phytoplankton"] <- BaseB["lg_phytoplankton"] * hind_clim[2:372,"phl"]
  # scene$forcing$ForcedBio[1:371,"sm_phytoplankton"] <- BaseB["sm_phytoplankton"] * hind_clim[2:372,"phs"]     
  
###############  
    
# Finalize scenario into scene_fit
  scene_fit <- scene4

################################################################################
# Run model
# Load Or Update fitting functions
  fup()
  
  run_fit   <- rsim.run(scene_fit, method='AB', years=hind_years)
  
  #render.show.fit(bal, scene_fit, run_fit, "test1")
################################################################################
# M0 Fitting (species one by one)
#  
# This function returns a list of biomass timeseries, species with timeseries, and combo
  bio_series   <- rsim.fit.list.bio.series(scene_fit)
  catch_series <- rsim.fit.list.catch.series(scene_fit)
  
# Set fitting weight of all biomass and catch timeseries to 0 
  scene_zero <- scene_fit
  scene_zero <- rsim.fit.set.bio.wt(scene_zero, bio_series$groups, bio_series$sources, 0)
  scene_zero <- rsim.fit.set.catch.wt(scene_zero, catch_series$groups, 0)
  

this_survey <- "ebs_race_shelf"  
survey_spp  <-rsim.fit.list.bio.bysurvey(scene_fit)[[this_survey]]

final_vartype <- NULL
final_species <- NULL
final_values  <- NULL

for (sp in survey_spp){

  scene_m <- rsim.fit.set.bio.wt(scene_zero, sp, this_survey, 1.0)
  scene_m <- rsim.fit.set.q(scene_m, sp, this_survey, years=1991)

# Vectors of parameters to fit.  Type of variable, species, and starting value  
  fit_vartype <- c("mzero") 
  fit_species <- c(sp)
  fit_start   <- as.numeric(scene_m$params$MzeroMort[sp])

  optim_out <- optim(fit_start, rsim.fit.run, method="Brent", lower=-5, upper=5,
                    species=fit_species, vartype=fit_vartype,
                    scene=scene_m, run_method="AB", run_years=hind_years)
  
  cat(sp," ",optim_out$par," ",optim_out$value,"\n"); flush.console()
  run_out <- rsim.fit.run(optim_out$par, fit_species, fit_vartype, scene_m, "AB", hind_years, verbose=T)
  X11(width=8,height=4)
  rsim.plot.full(scene_m, run_out, sp)    
  final_vartype <- c(final_vartype, fit_vartype)
  final_species <- c(final_species, fit_species)
  final_values  <- c(final_values, optim_out$par)
  
}  

###############################################################################
##
## BELOW THIS LINE ARE OUT-OF-ORDER TESTS NOT "FINAL" fitting sequence.
##
###############################################################################

  # rsim.fit.run produces a run that applies the above fit_vectors to the
  # scenario before running (scenario remains unchanged)
  run_0 <- rsim.fit.run(fit_start, fit_species, fit_vartype, scene_f3, "AB", hind_years, verbose=T) 
  X11(width=8,height=4)
  rsim.plot.full(scene_f3,run_0,"yf_sole")  
  #rsim.plot(run_0)
  
# By default, that biomass time series is an index (q internally calculated)
# Trying that fit.
  optim_f3 <- optim(fit_start, rsim.fit.run, method="Brent", lower=-3, upper=3,
                species=fit_species, vartype=fit_vartype,
                scene=scene_f3, run_method="AB", run_years=hind_years)
  optim_f3$par
  optim_f3$value
  
  run_f3 <- rsim.fit.run(optim_f3$par, fit_species, fit_vartype, scene_f3, "AB", hind_years, verbose=T) 
  X11(width=8,height=4)
  rsim.plot.full(scene_f3,run_f3,"yf_sole") 
  
  
# Trying fit with that same time series with a fixed q of 1 (very silly)
  scene_f4 <- rsim.set.fit.q(scene_f3, "yf_sole", "ebs_race_shelf", q=1)
  optim_f4 <- optim(fit_start, rsim.fit.run, method="Brent", lower=-3, upper=3,
                    species=fit_species, vartype=fit_vartype,
                    scene=scene_f4, run_method="AB", run_years=hind_years)
  optim_f4$par
  optim_f4$value

  run_f4 <- rsim.fit.run(optim_f4$par, fit_species, fit_vartype, scene_f4, "AB", hind_years, verbose=T) 
  X11(width=8,height=4)
  rsim.plot.full(scene_f4,run_f4,"yf_sole") 
  
# Trying fit with that same time series with a fixed q using 1991 and Ecopath starting value
  scene_f5 <- rsim.set.fit.q(scene_f3, "yf_sole", "ebs_race_shelf", years=1991)
  optim_f5 <- optim(fit_start, rsim.fit.run, method="Brent", lower=-3, upper=3,
                    species=fit_species, vartype=fit_vartype,
                    scene=scene_f5, run_method="AB", run_years=hind_years)
  optim_f5$par
  optim_f5$value
  run_f5 <- rsim.fit.run(optim_f5$par, fit_species, fit_vartype, scene_f5, "AB", hind_years, verbose=T) 
  X11(width=8,height=4)
  rsim.plot.full(scene_f5,run_f5,"yf_sole")   
  


  
    
  scene_test <- rsim.set.fit.q(scene3, "sandlance", "ebs_race_shelf", years=1991)
  
  
  
  run_f2   <- rsim.run(scene_f2, method='AB', years=hind_years)

  
  rsim.fit.obj(scene_f2,run_f2,verbose=F)
  rsim.fit.table(scene_f2,run_f2)
  
  X11(width=8,height=4)
  rsim.plot.full(scene_f2,run_f2,"sandlance")

    

  
  
  

  
  fit.res <- rsim.fit.run(test$par,
                          species=fit_species, vartype=fit_vartype,
                          scene=scene_f2, run_method="AB", run_years=hind_years,
                          verbose=T)

  X11(width=8,height=4)
  rsim.plot.full(scene_f2,fit.res,"sandlance")
  rsim.plot(fit.res)
  
  

  scene_test <- rsim.set.fit.q(scene3, "sandlance", "ebs_race_shelf", years=1985)
  scene_test <- rsim.set.fit.q(scene3, "sandlance", "ebs_race_shelf", years=1991:1993)
  
  #
  # #From a2 Repo
  #
  # # HISTORICAL CLIMATE
  # # hind_clim <- read.csv("climate/A2hind.csv", row.names=NULL)
  # hind_clim <- read.csv("climate/a2_hind.csv", row.names=NULL)
  # hind_clim <- hind_clim[253:624,]
  # BaseB        <- scene$params$B_BaseRef
  # names(BaseB) <- scene$params$spname
  # # Set groups to force   
  # scene$forcing$ForcedBio[1:371,"euphausiids"]      <- BaseB["euphausiids"]      * hind_clim[2:372,"eup"]
  # scene$forcing$ForcedBio[1:371,"copepods"]         <- BaseB["copepods"]         * hind_clim[2:372,"cop"]
  # scene$forcing$ForcedBio[1:371,"pelagic_microbes"] <- BaseB["pelagic_microbes"] * hind_clim[2:372,"mzl"]
  # scene$forcing$ForcedBio[1:371,"benthic_microbes"] <- BaseB["benthic_microbes"] * hind_clim[2:372,"mzl"]
  # scene$forcing$ForcedBio[1:371,"lg_phytoplankton"] <- BaseB["lg_phytoplankton"] * hind_clim[2:372,"phl"]
  # scene$forcing$ForcedBio[1:371,"sm_phytoplankton"] <- BaseB["sm_phytoplankton"] * hind_clim[2:372,"phs"]     
  
  
  # # Set-up bioenergetic forcing
  # bioen_sp_noceph <- bioen_sp[! bioen_sp %in% c("octopus","squids")]
  # # print(bioen_sp_noceph)
  # # consumption multiplier
  # if (cons == TRUE) {
  #   if(esm == "none" & ssp == "persist") {
  #     # hindcast
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedSearch[1:372,i] <- tdc_hind_bt[253:624,i]
  #     }
  #     # projection
  #     # do nothing? or default?
  #     # for(i in bioen_sp_noceph){
  #     #   scene$forcing$ForcedSearch[373:1308,i] <- 1
  #     # }
  #   }
  #   
  #   if(esm == "c" & ssp == 45){
  #     # hindcast
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedSearch[1:372,i] <- tdc_hind_bt[253:624,i]
  #     }
  #     # projection
  #     bioen_proj <- bioen_params(cmip, esm, ssp)
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedSearch[373:1068,i] <- bioen_proj[[1]][,i]
  #     }
  #     
  #   }
  #   else {
  #     # hindcast
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedSearch[1:372,i] <- tdc_hind_bt[253:624,i]
  #     }
  #     # projection
  #     bioen_proj <- bioen_params(cmip, esm, ssp)
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedSearch[373:1308,i] <- bioen_proj[[1]][,i]
  #     }
  #     
  #   }
  # }
  # 
  # # respiration multiplier ---------- ---------- #
  # if (resp == TRUE) {
  #   if(esm == "none" & ssp == "persist") {
  #     # hindcast
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedActresp[1:372,i] <- tdr_hind_bt[253:624,i]
  #     }
  #     # projection
  #     # do nothing? or default?
  #     # for(i in bioen_sp_noceph){
  #     #   scene$forcing$ForcedActresp[373:1308,i] <- 1
  #     # }
  #   }
  #   
  #   if(esm == "c" & ssp == 45){
  #     # hindcast
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedActresp[1:372,i] <- tdr_hind_bt[253:624,i]
  #     }
  #     # projection
  #     bioen_proj <- bioen_params(cmip, esm, ssp)
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedActresp[373:1068,i] <- bioen_proj[[2]][,i]
  #     }
  #   }
  #   else {
  #     # hindcast
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedActresp[1:372,i] <- tdr_hind_bt[253:624,i]
  #     }
  #     # projection
  #     bioen_proj <- bioen_params(cmip, esm, ssp)
  #     for(i in bioen_sp_noceph){
  #       scene$forcing$ForcedActresp[373:1308,i] <- bioen_proj[[2]][,i]
  #     }
  #     # print(scene$forcing$ForcedActresp[373:1308,"pcod_adu"])
  #   }
  # }
  
  