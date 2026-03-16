# Run the ACLIM 3.0 hindcast (1991–2021) for preliminary use with model fitting
library(Rpath)
library(dplyr)

# Fitting functions that will move to Rpath package are in the following files
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
  surveyfile_shelf    = "ebs/data/ebs_fitting_shelf_index.csv",
  surveyfile_slope    = "ebs/data/ebs_slope_fitting_index_aclim3_Q1.csv",
  surveyfile_ai_sebs  = "ebs/data/ai_sebs_fitting_index.csv",
  surveyfile_mammals  = "ebs/data/ebs_aclim3_mammal_index.csv",
  surveyfile_birds    = "ebs/data/seabird_colony_data.csv",
  surveyfile_zoop     = "ebs/data/ebs_zooplankton_pedigree.csv",
  # hindcast_ocean      = "ebs/climate/a2_hind.csv",
  hindcast_a3         = "ebs/climate/a3_hindcast.csv",
  hindcast_consump    = "ebs/data/cons_force_hind_bt_all.csv",
  hindcast_repsir     = "ebs/data/resp_force_hind_bt_all.csv",
  bioenergetics       = "ebs/data/ebs_aclim3_bioen.csv"
)  

# OLD FILE COPY DELETED - TOO COMPLICATED, UPDATE above datfiles by hand
#mdir <- "https://raw.githubusercontent.com/aw2113/aclim3_rpath/refs/heads/main/"  
#mdir <- "C:/Users/kerim.aydin/Work/src/aclim3_rpath/"

################################################################################

# Assemble biomass files into a single data table.  All biomass files must have
# the same columns (even colums unused by fitting must match)
bio_dat <- rbind(read.csv(datfiles$surveyfile_shelf),
                 read.csv(datfiles$surveyfile_slope),
                 read.csv(datfiles$surveyfile_ai_sebs),
                 read.csv(datfiles$surveyfile_mammals),
                 read.csv(datfiles$surveyfile_birds),
                 read.csv(datfiles$surveyfile_zoop))

catch_dat <- read.csv(datfiles$catchfile)

# hind_dat  <- read.csv(datfiles$hindcast_ocean)
hind_dat  <- read.csv(datfiles$hindcast_a3)

hind_cons <- read.csv(datfiles$hindcast_consump)
hind_resp <- read.csv(datfiles$hindcast_repsir)
bioen_pars<- read.csv(datfiles$bioenergetics)

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

# Climate forcing of lower trophics --------------------------------------------
scene4a <- scene3
# # 4.Add PP forcing (WGOA code from Bia modified for EBS) 

hind_clim <- hind_dat[253:624,]
# as of September 4,2025 hindcast prod only goes through 2019 (1:348)
scene4a$forcing$ForcedSearch[1:348,"lg_phytoplankton"] <- hind_clim$phl_anom_7605[1:348]
scene4a$forcing$ForcedSearch[1:348,"sm_phytoplankton"] <- hind_clim$phs_anom_7605[1:348]

# zero-trap: if biomass is <= 0 then make it = epsilon, otherwise leave it alone.
epsilon <- 1 * 10 ^ -15 # a really small number that is > 0
scene4a$forcing$ForcedSearch[, "lg_phytoplankton"]      <-
  ifelse(
    scene4a$forcing$ForcedSearch[, "lg_phytoplankton"] < epsilon, epsilon,
    scene4a$forcing$ForcedSearch[, "lg_phytoplankton"]
  )
scene4a$forcing$ForcedSearch[, "sm_phytoplankton"]      <-
  ifelse(
    scene4a$forcing$ForcedSearch[, "sm_phytoplankton"] < epsilon, epsilon,
    scene4a$forcing$ForcedSearch[, "sm_phytoplankton"]
  )


# Adding hindcast bioenergetics forcing
scene4b <- scene4a
bioen_sp <- bioen_pars$Species
# removing cephalopods for now
# bioen_sp_noceph <- bioen_sp[! bioen_sp %in% c("octopus","squids")]
# consumption multiplier
for(i in bioen_sp){
  scene4b$forcing$ForcedSearch[1:372,i] <- hind_cons[253:624,i]
}
# respiration multiplier
for(i in bioen_sp){
  scene4b$forcing$ForcedActresp[1:372,i] <- hind_resp[253:624,i]
}

# Manual adjustments
scene5 <- scene4b
# kamchatka exploitation rate: 0.02 from 2018 Kamchatcka SA  
scene5$fishing$ForcedFRate[as.character(1991:2010),"kamchatka"] <- 0.02

# Skate ID before 1999 is not good enough to use - removing from fit data
scene5$fitting$Biomass <- scene5$fitting$Biomass[
  !(scene5$fitting$Biomass$Group == "oth_skate" & scene5$fitting$Biomass$Year %in% 1991:1998),]
scene5$fitting$Biomass <- scene5$fitting$Biomass[
  !(scene5$fitting$Biomass$Group == "ak_skate" & scene5$fitting$Biomass$Year %in% 1991:1998),]  

# Remove biomass data series we aren't using/don't trust
removed_series <- c(
  "atka:ebs_race_shelf",      # atka seldom caught on ebs shelf, AI better survey
  "ben_zooplankton:ebs_race_shelf", # taxa in this group not reliably caught in any survey
  "capelin:ebs_race_shelf",         # capelin not reliably sampled by BTS's
  "gr_turbot_adu:ebs_race_slope",   # ebs shelf more comprehensive time series
  "infauna:ebs_race_shelf",         # infauna not well caught by any BTS
  "infauna:ebs_race_slope",         # infauna not well caught by any BTS
  "oth_rockfish:ebs_race_shelf", # rarely encountered in EBS surveys, AI better survey
  "ak_skate:ebs_race_slope",        # bulk of biomass on EBS shelf
  "pac_ocean_perch:ebs_race_shelf", # poorly sampled by EBS surveys, AI better survey
  "salmon_returning:ebs_race_shelf",# not well sampled by BTS's
  "squids:ebs_race_shelf",          # not well sampled by BTS 
  "deep_demersals:ebs_race_shelf",  # not well sampled by shelf BTS
  "north_rockfish:ebs_race_shelf", # poorly sampled by EBS surveys, AI better survey
  "rougheye_rock:ebs_race_shelf",  # poorly sampled by EBS surveys, AI better survey
  "mycto_bathy:ebs_race_shelf",     # not sampled by EBS shelf survey
  "pel_zooplankton:ebs_race_shelf", # not well sampled by BTS's
  "atka:ebs_race_slope",           # atka seldom caught in ebs, AI better survey
  "ben_zooplankton:ebs_race_slope", # taxa in this group not reliably caught in any survey
  "eelpouts:ebs_race_slope",        # different slope vs shelf spp
  "herring:ebs_race_slope",         # too spotty
  "jellyfish:ebs_race_slope",       # catching pelagic strays
  "king_crabs:ebs_race_slope",      # shelf suvey a better index
  # "opilio:ebs_race_slope",      # shelf suvey a better index
  "misc_crabs:ebs_race_slope",      # shelf suvey a better index
  "nr_sole:ebs_race_slope",         # primarily a shelf species
  "oth_flatfish:ebs_race_slope",    # primarily a shelf stock
  "oth_forage:ebs_race_slope",      # catching pelagic strays
  "pandalid:ebs_race_slope",        # poorly sampled by slope survey
  "pollock_adu:ebs_race_slope",     # poorly sampled by slope survey
  "salmon_returning:ebs_race_slope",# poorly sampled by BTS's
  "shallow_demersals:ebs_race_slope", # primarily a shelf group
  "sharks:ebs_race_slope",          # not well sampled by slope BTS
  "squids:ebs_race_slope",          # not well sampled by BTS
  "structural_epifauna:ebs_race_slope", # different species on slope with patchy distribution, poorly sampled
  "euphausiids:foci_summer_abundance", # poorly sampled by FOCI gear
  "euphausiids:foci_spring_abundance" # poorly sampled by FOCI gear
)

scene6 <- rsim.fit.remove.bio.timeseries(scene5, removed_series)


################################################################################
# Run model with all Q's internally calculated
#  fup()
#  #scene_fit <- scene5  
#  run_fit1   <- rsim.runandplot(scene5, hind_years, c(rpath.living(bal),rpath.detrital(bal)))

################################################################################
# Q adjustments
#
# Set Qs based on 1991 biomass
# scene6 <- scene5
# all_series <- strsplit(rsim.fit.list.bio.series(scene6)$all,":")
# for (i in 1:length(all_series)){
#   scene6 <- rsim.fit.set.q(scene6, all_series[[i]][1], all_series[[i]][2], years=1991, type="1991")  
# }
# run_fit2   <- rsim.runandplot(scene6, hind_years)

# set Qs based on first existing year of biomass data (1991 for many), but
# 2002 for slope series, some other series differ
scene7 <- scene6

all_series <- strsplit(rsim.fit.list.bio.series(scene7)$all,":")  
for (i in 1:length(all_series)){
  # loop through each time series
  species <- all_series[[i]][1]
  survey <- all_series[[i]][2]
  # find first year of time series
  first_year <- min(as.numeric(scene7$fitting$Biomass$Year[
    (scene7$fitting$Biomass$Group == species & scene7$fitting$Biomass$Source==survey)]))
  # set survey q so that Ecopath biomass equals first year of survey data
  scene7 <- rsim.fit.set.q(scene7, species, survey, years=first_year, type=first_year)
}

# END MODEL AND DATA SETUP
################################################################################  
#
scene_fit <- scene7


# At end of this script, following variables should be defined:
# bal
# unbal
# scene_fit
# hind_years
  
  
  
  