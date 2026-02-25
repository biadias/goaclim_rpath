# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 26 January 2026
#
# R/GOA_02_fitting_setup.R
# Purpose: Script for setting up the fitting scenarios for the GOA model.
# 
# Rpath_fitting/R/GOA_fitting_setup_a_la_carte.r to set up base scenario objects.
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
source("R/GOA_fitting_setup_a_la_carte.r")
source("R/fitting_aic.R")
source("R/ecofitting_plot_multi.R")

# ---------------------------------------------------------------------------- #
# 63 Parameters Setup ####
# ---------------------------------------------------------------------------- #
sv_groups_63par        <- c(
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
  "euphausiids",
  "rex_sole_adult",
  "pacific_cod_juvenile"
)

#Vectors

par63_vartype <- c(rep("preyvul", length(sv_groups_63par)),
                   rep("predvul",(length(sv_groups_63par)-1)))       #(bal$NUM_DEAD+2)))) bal$NUM_DEAD is the whole detritus pool
par63_species <- c(sv_groups_63par, 
                   sv_groups_63par[1:(length(sv_groups_63par)-1)])   #(bal$NUM_DEAD+2))]) bal$NUM_DEAD is the whole detritus pool
par63_start   <- rep(0,length(par63_species))

# ---------------------------------------------------------------------------- #
# 59 Parameters + 2 M0 Setup (Pollock and Cod) ####
# ---------------------------------------------------------------------------- #

mzero_groups_2     <- c("walleye_pollock_adult", "pacific_cod_adult" )
par61M02_groups        <- c(
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
  "euphausiids",
  "rex_sole_adult"
  #"pacific_cod_juvenile"
)
par61M02_vartype <- c(rep("mzero", length(mzero_groups_2)),
                      rep("preyvul", length(par61M02_groups)),
                      rep("predvul",(length(par61M02_groups)-1))) 
par61M02_species <- c(mzero_groups_2,
                      par61M02_groups, 
                      par61M02_groups[1:(length(par61M02_groups)-1)])

par61M02_start   <- c(scene_base$params$MzeroMort[mzero_groups_2],
                      rep(0,(length(par61M02_species)-length(mzero_groups_2))))

# ---------------------------------------------------------------------------- #
# 59 Parameters + 3 M0 Setup (Pollock, Cod, Arrowtooth) ####
# ---------------------------------------------------------------------------- #

mzero_groups_3     <- c("walleye_pollock_adult", "pacific_cod_adult","arrowtooth_flounder_adult" )
par59M03_groups        <- c(
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


par59M03_vartype <- c(rep("mzero", length(mzero_groups_3)),
                      rep("preyvul", length(par59M03_groups)),
                      rep("predvul",(length(par59M03_groups)-1))) 
par59M03_species <- c(mzero_groups_3,
                      par59M03_groups, 
                      par59M03_groups[1:(length(par59M03_groups)-1)])
par59M03_start   <- c(scene_base$params$MzeroMort[mzero_groups_3],
                      rep(0,(length(par59M03_species)-length(mzero_groups_3))))

