# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 03 March 2026
#
# code/01_load_best_model.R
# Purpose: Loading the best scenario to the set up scene_bioen
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# 1. Load the best model scenario ####
# ---------------------------------------------------------------------------- #
result_59_M04 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par.rds")
source("code/00_setup_forecast.R")

mzero_groups_4     <- c("walleye_pollock_adult", "pacific_cod_adult",
                        "arrowtooth_flounder_adult", "deep_water_flatfish" )
par59M04_groups        <- c(
  "other_skates",
  "walleye_pollock_adult",
  "pandalid_shrimp",
  "pacific_sandlance",
  "tanner_crab" ,
  "benthic_zooplankton",
  "pacific_halibut_adult",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult" ,
  "salmon_returning",
  "shelf_demersal_fish",
  "shelf_forage_fish",
  "offal",#only preyvul for this group
  "infauna",
  "pacific_herring_juvenile",
  "sablefish_adult",
  "pacific_herring_adult",
  "arrowtooth_flounder_juvenile",
  "mysids",
  "deep_water_flatfish",
  "pacific_capelin",
  "large_microzooplankton" ,
  "pacific_ocean_perch_adult",
  "squid", 
  "shallow_water_flatfish",
  "motile_epifauna",
  "octopus",
  "nonpandalid_shrimp",
  "walleye_pollock_juvenile",
  "euphausiids"
  #"rex_sole_adult",
  #"pacific_cod_juvenile"
)

par59M04_vartype <- c(rep("mzero", length(mzero_groups_4)),
                      rep("preyvul", length(par59M04_groups)),
                      rep("predvul",(length(par59M04_groups)-1))) 
par59M04_species <- c(mzero_groups_4,
                      par59M04_groups, 
                      par59M04_groups[1:(length(par59M04_groups)-1)])

par59M04_start   <- c(scene_bioen$params$MzeroMort[mzero_groups_4],
                      rep(0,(length(par59M04_species)-length(mzero_groups_4))))


# create scenario object to apply the fit vulnerabilities
combined_pars <- c(
  result_59_M04$Bioen$opt_object_s2$par, 
  result_59_M04$Bioen$opt_object_s1$par
)

scene_bioen_best <- scene_bioen


scene_bioen_best$params <- rsim.fit.apply(
  values = combined_pars, # Uses the fitted vector from your RDS
  species = par59M04_species,  # Your combined species vector (M0 + prey + pred)
  vartype = par59M04_vartype,  # Your combined vartype vector
  scene.params = scene_bioen$params
)
