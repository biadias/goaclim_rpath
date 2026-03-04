# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 25 February 2026
#
# code/00_setup_forecast.R
# Purpose: Script for setting up the ssp scenarios for forecasting for the best model
# 
# Rpath_fitting/ has all the fitting repo as a subtree
# ---------------------------------------------------------------------------- #

library(Rpath)

# ---------------------------------------------------------------------------- #
# 1) Load the best model ####
# ---------------------------------------------------------------------------- #

goa_fit <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par.rds")
scene_bioen <- goa_fit$Bioen

# NOTE: I manually set the years in the bioenergetic projections script to 1990:2089, so we need to match that here.
# all_years <- 1:1308
# hind_years <- 1:372
# managed_sp <- c("pollock", "pcod", "arrowtooth", ...) # GOA specific managed species
# F_equil <- 0 # Or your calculated F values

# Source your bioenergetics calculation script (adapted for GOA species/temps)
# source("Rpath_fitting/GOA/bioenergetic_projections.r)") 
# source("code/F_clim_sim_scene_prim_prod.R")

# ---------------------------------------------------------------------------- #
# 2) Build GFDL-Specific Climate Scenarios ####
# ---------------------------------------------------------------------------- #
# Note: You may need to replace 'scene_refbase' inside your F_clim_sim_scene 
# function with 'scene_base' so it references the unfished/base GOA state properly.

# A) Persistence (Base expectation holding recent climate constant)
scene_gfdl_persist <- F_clim_sim_scene(scene = scene_base,
                                       cmip = "CMIP6", 
                                       esm = "none", 
                                       ssp = "persist",
                                       cons = TRUE, resp = TRUE, buf = TRUE)

# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_gfdl_126 <- F_clim_sim_scene(scene = scene_base,
                                   cmip = "CMIP6", 
                                   esm = "GFDL", 
                                   ssp = 126,
                                   cons = TRUE, resp = TRUE, buf = TRUE)

# C) GFDL SSP 585 (High emission fossil-fueled development scenario)
scene_gfdl_585 <- F_clim_sim_scene(scene = scene_base,
                                   cmip = "CMIP6", 
                                   esm = "GFDL", 
                                   ssp = 585,
                                   cons = TRUE, resp = TRUE, buf = TRUE)

# ---------------------------------------------------------------------------- #
# 3) Run Forecast Simulations
# ---------------------------------------------------------------------------- #
# Run the dynamic simulations using Adams-Bashforth (AB) method

forecast_gfdl_persist <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_126     <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_585     <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)

# ---------------------------------------------------------------------------- #
# 4) Calculate B0 and Reference Points (Optional)
# ---------------------------------------------------------------------------- #
# source("R/goa_bref.R")

# Calculate Unfished Biomass (B0)
bzero_persist <- bzero_func(scene_gfdl_persist, managed_sp, F_equil, hind_years, all_years)
bzero_g126    <- bzero_func(scene_gfdl_126, managed_sp, F_equil, hind_years, all_years)
bzero_g585    <- bzero_func(scene_gfdl_585, managed_sp, F_equil, hind_years, all_years)

# Calculate Biomass Reference Points
persist_bref <- goa_bref(bzero_persist)
g126_bref    <- goa_bref(bzero_g126)
g585_bref    <- goa_bref(bzero_g585)

# Save Reference Points
# save(bzero_persist, bzero_g126, bzero_g585, file = "data/brp/goa_gfdl_brps.RData")

