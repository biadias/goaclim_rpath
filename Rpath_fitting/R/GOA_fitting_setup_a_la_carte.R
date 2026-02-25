#------------------------------------------------------------------------------#
#AUTHORS: Kerim Aydin, Bridget Ferrisss, Bia Dias, Andy Whitehouse
#AFFILIATIONS: Alaska Fisheries Science Center / CICOES University of Washington
#E-MAIL OF CORRESPONDENCE AUTHOR: kerim.aydin@noaa.gov bridget.ferriss@noaa.gov,
# bia.dias@noaa.gov
#
# #### Fit WGOA model with separate options (scenes) to turn on and off bioenergetics and primary production forcing
#code modified from "GOA_fitting_groups.R retired code on google drive folder"
#------------------------------------------------------------------------------#



#library(devtools)
#install_github('NOAA-EDAB/Rpath', ref="dev")
library(Rpath)
library(dplyr)
library(tidyverse)
library(janitor)

# Support function files that may be moved to Rpath when finalized
fup <- function() {
  #source("R/xml_convert.r") # comment out this line
  source("R/merge_ecofitting.R")
  source("R/ecofitting_plots.R")
}
fup()

#### A. LOAD WGOA and EGOA models #############################################################

##### 1.Load and test WGOA model ####
WGOA_EwE_file <- "GOA/WGOA_23Jan2026_simpleDet.eiixml" 
# changed following two lines
w.unbal  <-  create.rpath.from.eiixml(WGOA_EwE_file)  #xml_unbal(WGOA_EwE_file)
w.unbal  <-  rpath.stanzas(w.unbal)
w.bal    <-  rpath(w.unbal) # balanced
# Test equilibrium if desired
#w.scene0 <- rsim.scenario(w.bal, w.unbal, years=1990:2089)
#w.run0   <- rsim.run(w.scene0, method="AB", years = 1990:2089)
#rsim.plot(w.run0) # should be flatline

##### 2.Load and test EGOA model - commenting out for now ####
#    EGOA_EwE_file <- "GOA/EGOA_20250514_simpleDet.eiixml"
#    e.unbal  <- xml_unbal(EGOA_EwE_file)
#    e.bal    <- rpath(e.unbal) # balanced
#    # Test equilibrium if desired
#      # e.scene0 <- rsim.scenario(e.bal, e.unbal, years=1990:2089)
#      # e.run0   <- rsim.run(e.scene0, method="AB", years = 1990:2089)
#      # rsim.plot(e.run0) # should be flatline

#### B. WGOA MODEL AND DATA SETUP #############################################################
datfiles <- list(
  catchfile             = "GOA/wgoa_data_rpath_fitting/wgoa_catches_ft_cas_long_haladded.csv",
  # Halibut from IPHC added
  surveyfile_shelf      = "GOA/wgoa_data_rpath_fitting/wgoa_race_biomass_ts_fitting_index_v2_tons_ka.csv",
  surveyfile_nonrace    = "GOA/wgoa_data_rpath_fitting/wgoa_nonrace_biomass_ts_fitting_index.csv",
  surveyfile_gak        = "GOA/wgoa_data_rpath_fitting/gak_zooplankton_b_ts_v2.csv",
  surveyfile_ecofoci    = "GOA/wgoa_data_rpath_fitting/goaecofoci_zooplankton_b_ts_v2.csv",
  surveyfile_wintshelik = "GOA/wgoa_data_rpath_fitting/wgoa_pollock_shelikof_v2_biomass_ts_fitting_index_v2_tons_ka.csv"
  #surveyfile_slope    = "data/ebs_slope_fitting_index_aclim3_Q1.csv",
  #surveyfile_mammals  = "data/ebs_aclim3_mammal_index.csv",
  #surveyfile_birds    = "data/seabird_colony_data.csv",
  #hindcast_ocean      = "climate/a2_hind.csv"
)

ecofoci <- read.csv(datfiles$surveyfile_ecofoci) %>%
  mutate(
    Group = case_when(
      Group == "Cop" ~ "small_copepods",
      Group == "NCa" ~ "large_copepods",
      Group == "Eup" ~ "euphausiids"
    )
  )

wintshelik <- read.csv(datfiles$surveyfile_wintshelik) %>%
  mutate(Value = as.numeric(Value))


##### 1.Read biomass and catch timeseries into data frames ####
#do any needed data cleanup  #
bio_dat <- rbind(
  read.csv(datfiles$surveyfile_shelf),
  read.csv(datfiles$surveyfile_gak),
  wintshelik,
  ecofoci
)
bio_dat$Group <- make_clean_names(bio_dat$Group, allow_dupes = TRUE)


nr <- read.csv(datfiles$surveyfile_nonrace)
nr$Group <- make_clean_names(nr$Group, allow_dupes = TRUE)
nr$Stdev = nr$Value * nr$CV
bio_dat <- rbind(bio_dat, nr)

catch_dat <- read.csv(datfiles$catchfile)
catch_dat$Group <- make_clean_names(catch_dat$Group, allow_dupes = TRUE)

##### 2.Balance model for fitting  ####
unbal <- w.unbal
bal   <- rpath(unbal)

# Setup base scenario
hind_years <- 1990:2020
scene0 <- rsim.scenario(bal, unbal, years = hind_years)
# Read in biomass to fit from data frame
scene1 <- read.fitting.biomass.cv(scene0, bio_dat)
# Read in catch to fit from data frame. This does not apply catch to catch forcing.
scene2 <- read.fitting.catch(scene1, catch_dat)

# The fitcatch.to.forcecatch function has three effects:
# (1) turns off gear Effort,
# (2) moves Ecopath baseline fishing to Frate forcing for all years
# (3) adds forcecatch for years where catch is provided, 0ing the baseline Frate in those years.
scene3 <- fitcatch.to.forcecatch(scene2, bal)

# Currently forced catch and forced F do not go to detrital groups (discard/offal),
# so need to keep discards/offal locked at equilibrium for simulations.
scene3$forcing$ForcedBio[, "discards"] <- bal$Biomass["discards"]
scene3$forcing$ForcedBio[, "offal"]    <- bal$Biomass["offal"]

###revised WGOA set survey biomass Q's to fixed first year of time series or index for all indices (same as EBS)
# Originally Set all survey biomass q's to fixed 90-93 avg only for surveys (BTS) and not indices (Q is estimated)
##  scene4 <- scene3
##   this_survey <- "race_wgoa"
##   survey_spp  <-rsim.fit.list.bio.bysurvey(scene4)[[this_survey]]
##   for (sp in survey_spp){
##     scene4 <- rsim.fit.set.q(scene4, sp, this_survey, q=NULL, years=1990:1993, type="90-93")
##   }
scene4 <- scene3
all_series <- strsplit(rsim.fit.list.bio.series(scene4)$all, ":")
for (i in 1:length(all_series)) {
  # loop through each time series
  species <- all_series[[i]][1]
  survey <- all_series[[i]][2]
  # find first year of time series
  first_year <- min(as.numeric(scene4$fitting$Biomass$Year[(scene4$fitting$Biomass$Group == species &
                                                              scene4$fitting$Biomass$Source == survey)]))
  # set survey q so that Ecopath biomass equals first year of survey data
  scene4 <- rsim.fit.set.q(scene4, species, survey, years = first_year, type =
                             first_year)
}

##### 3. Remove biomass data series we aren't using/don't trust ####
removed_series <- c(
  "atka_mackerel:race_wgoa",
  "lingcod:race_wgoa",
  "pacific_hake:race_wgoa",
  "pacific_herring_adult:craig_index",
  "salmon_shark:race_wgoa",
  "walleye_pollock_adult:race_wgoa",
  "miscellaneous_deep_sea_fish:race_wgoa", #removing this time series_replacing by shelikoff index
  "squid:race_wgoa",
  "nonpandalid_shrimp:race_wgoa",
  "other_gelatinous_zooplankton:race_wgoa",
  "squid:race_wgoa",
  "infauna:race_wgoa",
  "mysids:race_wgoa",
  "euphausiids:race_wgoa",
  "euphausiids:summer_ecofoci",
  "large_copepods:summer_ecofoci",
  "small_copepods:summer_ecofoci"
)

scene5 <- rsim.fit.remove.bio.timeseries(scene4, removed_series)

##### 4.Rsim & plot ####
plot.species <- c(rpath.living(bal),rpath.detrital(bal))  
# Run and calculate without changing any fitting values  
run.base <- rsim.fit.run(NA, NA, NA, scene=scene5, run_method="AB", 
                         run_years=hind_years, verbose=T)
# test to show output negative log likelihood
run.base.nll <- rsim.fit.run(NA, NA, NA, scene=scene5, run_method="AB", 
                             run_years=hind_years, verbose=F); run.base.nll # default verbosity



##Base scene (scene_base) of fit model (no bioenergetics, no primary production forcing)
scene_base <- scene5

# ============================================================================ #

##### 5. Only add bioenergetics forcing ####
# Start with base scene and add bioenergetics forcing

scene_bioen <- scene_base

source("GOA/wgoa_bioenergetics_code/wgoa_add_bioen_to_scene.r")
scene_bd <- wgoa_add_bioenergetics("ssp126")
# Since scene_bd was created in bioenergetics functions without having other
# forcing added, copy the affected search parameters from scene_bd to scene_bioen
hindmonths <- dim(scene_bioen$forcing$ForcedSearch)[1]
scene_bioen$forcing$ForcedSearch[1:hindmonths, ]  <- scene_bd$forcing$ForcedSearch[1:hindmonths, ]
scene_bioen$forcing$ForcedActresp[1:hindmonths, ] <- scene_bd$forcing$ForcedActresp[1:hindmonths, ]

#=================================================================#
##### 6. Only add primary production forcing ####
# Start with base scene and add primary production forcing

scene_primprod <- scene_base

hindcast_pp <- read.csv("GOA/wgoa_data_rpath_fitting/Long_WGOA_NPZ_PP_monthly_B_added.csv")
ppl_totP <- as.numeric(bal$Biomass["large_phytoplankton"] * bal$PB["large_phytoplankton"])
ppl_force <- hindcast_pp %>%
  filter(simulation == "ssp126" &
           varname == "prod_PhL" & year %in% hind_years) %>%
  group_by(month) %>%
  mutate(Pmean = mean(P_tkm2_mon)) %>%
  ungroup() %>%
  mutate(P_anom = P_tkm2_mon / Pmean,
         P_add  = 1 + (P_tkm2_mon - Pmean) / ppl_totP)

pps_totP <- as.numeric(bal$Biomass["small_phytoplankton"] * bal$PB["small_phytoplankton"])
pps_force <- hindcast_pp %>%
  filter(simulation == "ssp126" &
           varname == "prod_PhS" & year %in% hind_years) %>%
  group_by(month) %>%
  mutate(Pmean = mean(P_tkm2_mon)) %>%
  ungroup() %>%
  mutate(P_anom = P_tkm2_mon / Pmean,
         P_add  = 1 + (P_tkm2_mon - Pmean) / pps_totP)

scene_primprod$forcing$ForcedSearch[, "large_phytoplankton"] <- ppl_force$P_anom #ppl_force$P_add
scene_primprod$forcing$ForcedSearch[, "small_phytoplankton"] <- pps_force$P_anom #pps_force$P_add
#scene_fit$forcing$ForcedBio[,"large_phytoplankton"] <- scene_fit$params$B_BaseRef["large_phytoplankton"] * ppl_force$anomaly_ratio
#scene_fit$forcing$ForcedBio[,"small_phytoplankton"] <- scene_fit$params$B_BaseRef["small_phytoplankton"] * pps_force$anomaly_ratio


# ============================================================================ #
##### 7.Add both primary production and bioenergetics ####
# Start with the scene_bioen (adds to all groups) then add primary production forcing (only to the phyto groups)
#if start with scene_primprod and then add scene_bioen it overwrites the primprod groups
scene_full <- scene_bioen

#see scene_primprod above for all the calcs leading up to this
scene_full$forcing$ForcedSearch[, "large_phytoplankton"] <- ppl_force$P_anom #ppl_force$P_add
scene_full$forcing$ForcedSearch[, "small_phytoplankton"] <- pps_force$P_anom #pps_force$P_add


#################################################################################
# END MODEL AND DATA SETUP
################################################################################
#
