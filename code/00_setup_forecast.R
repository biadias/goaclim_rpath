# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 25 February 2026
#
# code/00_setup_forecast.R
# Purpose: Script for setting up the ssp scenarios for forecasting for the best model
# 1) We need to recreate the GOA_fitting_setup_a_la_carte.R and expand the years 
# beyond the hindcast to match the bioenergetic projections script (1990:2099).
# 2) Then we need to import the vul and M0 from the best fitted model into the scenarios for the forecast.
#
# Rpath_fitting/ has all the fitting repo as a subtree
# ---------------------------------------------------------------------------- #


#library(devtools)
#install_github('NOAA-EDAB/Rpath', ref="dev")
library(Rpath)
library(dplyr)
library(tidyverse)
library(janitor)
library(here)

# Support function files that may be moved to Rpath when finalized
fup <- function() {
  #source("R/xml_convert.r") # comment out this line
  source("Rpath_fitting/R/merge_ecofitting.R")
  source("Rpath_fitting/R/ecofitting_plots.R")
}
fup()

# ---------------------------------------------------------------------------- #
# 1.Load and test WGOA model ####
# ---------------------------------------------------------------------------- #
WGOA_EwE_file <- "Rpath_fitting/GOA/WGOA_19March2026_simpleDet.eiixml" 
# changed following two lines
w.unbal  <-  create.rpath.from.eiixml(WGOA_EwE_file)  
w.unbal  <-  rpath.stanzas(w.unbal)
w.bal    <-  rpath(w.unbal) # balanced
# Test equilibrium if desired
#w.scene0 <- rsim.scenario(w.bal, w.unbal, years=1990:2099)
#w.run0   <- rsim.run(w.scene0, method="AB", years = 1990:2099)
#rsim.plot(w.run0) # should be flatline

## a. WGOA MODEL AND DATA SETUP ####
datfiles <- list(
  catchfile             = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/wgoa_catches_ft_cas_long_haladded.csv",# Halibut from IPHC added
  surveyfile_shelf      = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/wgoa_race_biomass_ts_fitting_index_v2_tons_ka.csv",
  surveyfile_nonrace    = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/wgoa_nonrace_biomass_ts_fitting_index.csv",
  surveyfile_gak        = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/gak_zooplankton_b_ts_v2.csv",
  surveyfile_ecofoci    = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/goaecofoci_zooplankton_b_ts_v2.csv",
  surveyfile_wintshelik = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/wgoa_pollock_shelikof_v2_biomass_ts_fitting_index_v2_tons_ka.csv",
  surveyfile_sable_ll   = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/wgoa_sablefish_ll_v2_biomass_ts_fitting_index_v3.csv"
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

sable_ll <- read.csv(datfiles$surveyfile_sable_ll) %>%
  mutate(Value = as.numeric(Value))


## b. Read biomass and catch timeseries into data frames ####
bio_dat <- rbind(
  read.csv(datfiles$surveyfile_shelf),
  read.csv(datfiles$surveyfile_gak),
  wintshelik,
  sable_ll,
  ecofoci
)
bio_dat$Group <- make_clean_names(bio_dat$Group, allow_dupes = TRUE)


nr <- read.csv(datfiles$surveyfile_nonrace)
nr$Group <- make_clean_names(nr$Group, allow_dupes = TRUE)
nr$Stdev = nr$Value * nr$CV
bio_dat <- rbind(bio_dat, nr)

catch_dat <- read.csv(datfiles$catchfile)
catch_dat$Group <- make_clean_names(catch_dat$Group, allow_dupes = TRUE)

# ---------------------------------------------------------------------------- #
# 2.Balance model for fitting  ####
# ---------------------------------------------------------------------------- #
unbal <- w.unbal
bal   <- rpath(unbal)

# Useful groupings 
all_living   <- bal$Group[1:bal$NUM_LIVING]
all_detritus <- bal$Group[(bal$NUM_LIVING+1):(bal$NUM_LIVING+bal$NUM_DEAD)]
all_gears    <- bal$Group[(bal$NUM_LIVING+bal$NUM_DEAD+1):(bal$NUM_GROUPS)]
# Equilibrium F and Biomass
F_equil <- (rowSums(bal$Landings)+rowSums(bal$Discards))/(bal$Biomass)
names(F_equil) <- bal$Group
F_equil  <- F_equil[c(all_living,all_detritus)]
# F_zero scenario
F_zero <- rep(0,(bal$NUM_LIVING + bal$NUM_DEAD + 1))
names(F_zero) <- c("Outside",bal$Group[c(all_living,all_detritus)])
# Ecopath equilibrium biomass
B_equil  <- bal$Biomass; names(B_equil) <- bal$Group; 
B_equil  <- B_equil[c(all_living,all_detritus)] 


# Setup base scenario years
hind_years <- 1991:2020 #1:360 here and 133:492 months tstep from ssp126_wide_WGOA_temp_1000.csv
fore_years <- 2021:2099 #361:1308 here and 493:1440 months tstep from ssp126_wide_WGOA_temp_1000.csv
all_years  <- c(hind_years,fore_years) #1:1308 here and 133:1440 tstep

# Setup base scenario
scene0 <- rsim.scenario(bal, unbal, years = all_years)
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

# ---------------------------------------------------------------------------- #
# 3. Remove biomass data series we aren't using/don't trust ####
# ---------------------------------------------------------------------------- #
removed_series <- c(
  "atka_mackerel:race_wgoa",
  "lingcod:race_wgoa",
  "pacific_hake:race_wgoa",
  "pacific_herring_adult:craig_index",
  "salmon_shark:race_wgoa",
  "walleye_pollock_adult:race_wgoa",#removing this time series_replacing by shelikoff index
  "sablefish_adult:race_wgoa",
  "miscellaneous_deep_sea_fish:race_wgoa", 
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
# ---------------------------------------------------------------------------- #
# 4.Rsim & plot ####
# ---------------------------------------------------------------------------- #
plot.species <- c(rpath.living(bal),rpath.detrital(bal))  
# Run and calculate without changing any fitting values  
run.base <- rsim.fit.run(NA, NA, NA, scene=scene5, run_method="AB", 
                         run_years=all_years, verbose=T)
# test to show output negative log likelihood
run.base.nll <- rsim.fit.run(NA, NA, NA, scene=scene5, run_method="AB", 
                             run_years=all_years, verbose=F); run.base.nll # default verbosity



##Base scene (scene_base) of fit model (no bioenergetics, no primary production forcing)
scene_base <- scene5

# ---------------------------------------------------------------------------- #
# 5. Bioenergetics forcing ####
# ---------------------------------------------------------------------------- #
# Start with base scene and add bioenergetics forcing

scene_bioen <- scene_base
source("Rpath_fitting/GOA/wgoa_bioenergetics_code/wgoa_add_bioen_to_scene.R")
scene_bd <- wgoa_add_bioenergetics("ssp126")
# Since scene_bd was created in bioenergetics functions without having other
# forcing added, copy the affected search parameters from scene_bd to scene_bioen
hind_years
start_hind_year <- min(hind_years)
end_hind_year <- max(hind_years)
start_hind_date <- ymd(paste0(start_hind_year,"-01-01"))
end_hind_date <- ymd(paste0(end_hind_year,"-12-31"))
end_hind <- interval(start_hind_date, end_hind_date) %/% months(1) +1

fore_years
start_fore_year <- min(fore_years)
end_fore_year <- max(fore_years)
start_fore_date <- ymd(paste0(start_fore_year,"-01-01"))
end_fore_date <- ymd(paste0(end_fore_year,"-12-31"))
end_fore <- interval(start_fore_date, end_fore_date) %/% months(1) +1

scene_bioen$forcing$ForcedSearch[1:end_hind, ]  <- scene_bd$forcing$ForcedSearch[1:end_hind, ]
scene_bioen$forcing$ForcedActresp[1:end_hind, ] <- scene_bd$forcing$ForcedActresp[1:end_hind, ]

# ---------------------------------------------------------------------------- #
# 6. Only add primary production forcing ####
# ---------------------------------------------------------------------------- #
# Start with base scene and add primary production forcing

scene_primprod <- scene_base
hindcast_pp <- read.csv("Rpath_fitting/GOA/wgoa_data_rpath_fitting/Long_WGOA_NPZ_PP_monthly_B_added.csv")
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

scene_primprod$forcing$ForcedSearch[1:end_hind, "large_phytoplankton"] <- ppl_force$P_anom #ppl_force$P_add
scene_primprod$forcing$ForcedSearch[1:end_hind, "small_phytoplankton"] <- pps_force$P_anom #pps_force$P_add
#scene_fit$forcing$ForcedBio[,"large_phytoplankton"] <- scene_fit$params$B_BaseRef["large_phytoplankton"] * ppl_force$anomaly_ratio
#scene_fit$forcing$ForcedBio[,"small_phytoplankton"] <- scene_fit$params$B_BaseRef["small_phytoplankton"] * pps_force$anomaly_ratio


# ---------------------------------------------------------------------------- #
# 7.Add both primary production and bioenergetics ####
# ---------------------------------------------------------------------------- #

# Start with the scene_bioen (adds to all groups) then add primary production forcing (only to the phyto groups)
#if start with scene_primprod and then add scene_bioen it overwrites the primprod groups
scene_full <- scene_bioen
#see scene_primprod above for all the calcs leading up to this
scene_full$forcing$ForcedSearch[1:end_hind, "large_phytoplankton"] <- ppl_force$P_anom #ppl_force$P_add
scene_full$forcing$ForcedSearch[1:end_hind, "small_phytoplankton"] <- pps_force$P_anom #pps_force$P_add


# ---------------------------------------------------------------------------- #
# END MODEL AND DATA SETUP
# ---------------------------------------------------------------------------- #
