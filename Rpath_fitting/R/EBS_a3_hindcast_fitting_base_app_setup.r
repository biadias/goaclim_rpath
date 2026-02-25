# Run the ACLIM 2.0 hindcast (1991–2021) for preliminary use with model fitting
library(Rpath)

# TMP functions that will move to the Rpath package when finished
fup <- function(){
  source("R/merge_ecofitting.R")
  source("R/ecofitting_plots.R")
}

fup()

# Base Ecopath Files 
  Ebase <- "model/ebs_aclim3_76bio_base.csv"  # Base biomass, production, fishing, etc.
  Ediet <- "model/ebs_aclim3_76bio_diet.csv"  # Diet matrix
  Eped  <- "model/ebs_aclim3_76bio_pedigree.csv"  # Data pedigree = quality of input data
  Estz  <- "model/ebs_aclim3_76bio_stanzas.csv"  # Stanzas
  Estg  <- "model/ebs_aclim3_76bio_stanza_groups.csv" # Stanza groups

# FITTING DATA FILES
  datfiles <- list(
    catchfile           = "data/ebs_aclim3_76bio_catch_long.csv",
    surveyfile_shelf    = "data/ebs_shelf_6mar2025_fitting_index_aclim3_Q1.csv",
    surveyfile_slope    = "data/ebs_slope_6mar2025_fitting_index_aclim3_Q1.csv",
    surveyfile_mammals  = "data/ebs_aclim3_mammal_index.csv",
    surveyfile_birds    = "data/seabird_colony_data.csv")  

################################################################################  
# Download files if needed - set to false by default (doesn't update from repo)
# need to have repo permissions to use these files directly.  Assumes repo has
# the model/ and data/ subdirectories with files above.

mdir <- "https://raw.githubusercontent.com/aw2113/aclim3_rpath/refs/heads/main/"  
UPDATE_FILES <- FALSE
if (UPDATE_FILES){
  download.file(paste(mdir,Ebase,sep=''), Ebase)
  download.file(paste(mdir,Ediet,sep=''), Ediet)
  download.file(paste(mdir,Eped,sep=''), Eped)  
  download.file(paste(mdir,Estz,sep=''), Estz)
  download.file(paste(mdir,Estg,sep=''), Estg)
  for(i in datfiles){download.file(paste(mdir,i,sep=''),i)} 
}
################################################################################

# Assemble biomass files into a single data table.  All biomass files must have
# the same columns (even colums unused by fitting must match)
  bio_dat <- rbind(read.csv(datfiles$surveyfile_shelf),read.csv(datfiles$surveyfile_slope),
                   read.csv(datfiles$surveyfile_mammals),read.csv(datfiles$surveyfile_birds))
  
# Basic Ecopath balance
  unbal <- rpath.stanzas(read.rpath.params(Ebase, Ediet, Eped, Estg, Estz)) # unbalanced
  bal   <- rpath(unbal) # balanced

# Setup years and base ecosim scenario
  hind_years <- 1991:2023
  scene0 <- rsim.scenario(bal, unbal, years = hind_years)

# Read in biomass to fit 
  scene1 <- read.fitting.biomass(scene0, bio_dat)

# Read in catch to fit.  This does not on its own apply read-in catch to catch forcing.  
  scene2 <- read.fitting.catch(scene1, datfiles$catchfile)

# fitcatch.to.forcecatch (1) turns off gear Effort, (2) moves Ecopath baseline 
# fishing to Frate forcing for all years and (3) adds forcecatch for years where
# catch is provided, 0ing the baseline Frate in those years.  
  scene3 <- fitcatch.to.forcecatch(scene2, bal)
# Currently forced catch and forced F do not go to detritus, so need to keep
# that group constant.
  scene3$forcing$ForcedBio[,"discards_offal"] <- bal$Biomass["discards_offal"]  
  
# Finalize scenario into scene_fit
  scene_fit <- scene3

################################################################################
  
# Run model

  run_fit   <- rsim.run(scene_fit, method='AB', years=hind_years)
  
  #render.show.fit(bal, scene_fit, run_fit, "test1")
  
# This function returns a list of biomass timeseries, species with timeseries, and combo
  bio_series   <- rsim.fit.list.bio.series(scene_fit)
  catch_series <- rsim.fit.list.catch.series(scene_fit)
  
  # Set fitting weight of all biomass and catch timeseries to 0 
  scene_f0 <- scene_fit
  scene_f1 <- rsim.fit.set.bio.wt(scene_f0, bio_series$groups, bio_series$sources, 0)
  scene_f2 <- rsim.fit.set.catch.wt(scene_f1, catch_series$groups, 0)
  
  # Now we're going to fit sandlance, just to the sandlance survey data
  # First turn on the weight for that timeseries
  scene_f3 <- rsim.fit.set.bio.wt(scene_f2, "yf_sole", "ebs_race_shelf", 1.0)

# Vectors of parameters to fit.  Type of variable, species, and starting value  
  fit_vartype <- c("mzero") 
  fit_species <- c("yf_sole")
  fit_start   <- as.numeric(scene_fit$params$MzeroMort["yf_sole"])

  fit_start
  
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
  
  # #From a2 Repo
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
  
  