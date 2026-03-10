# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 03 March 2026
#
# code/02_run_sim.R
# Purpose: Running the F_clim_sim_scene function to set up the scenarios and then 
# running the simulations for the best model scenario (GOA_fit_results_59M04par.rds) 
# for the GFDL SSP 126 and SSP 585 scenarios, as well as a persistence scenario.
# ---------------------------------------------------------------------------- #
source("code/00_setup_forecast.R")
source("code/01_load_best_model.R") # Loads the best model scenario and sets up scene_bioen
source("Rpath_fitting/GOA/wgoa_bioenergetics_code/bioenergetic_projections.r") # Make sure this is adapted for GOA and matches the years used in the projections
source("code/Function_F_clim_sim_scene.R")




# ---------------------------------------------------------------------------- #
# 1. Build GFDL-Specific Climate Scenarios ####
# ---------------------------------------------------------------------------- #

# Persistence (Base expectation holding recent climate constant)

scene_gfdl_persist <- F_clim_sim_scene(scene = scene_bioen_best,
                                            ssp = "persist",
                                            cons = TRUE, resp = TRUE, buf = FALSE,
                                            bioen_sp = bioen_sp,
                                            tdc_hind = tdc_hind_bt, 
                                            tdr_hind = tdr_hind_bt,
                                            climate_dir= "data/climate/",
                                            hind_yrs = hind_years,
                                            proj_yrs = 2021:2099,
                                            hind_data_start_yr = 1991,
                                            climate_data_start_yr = 1980,
                                            verbose = TRUE) 


# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_gfdl_126 <- F_clim_sim_scene(scene = scene_bioen_best,
                                        ssp = "126",
                                        cons = TRUE, resp = TRUE, buf = FALSE,
                                        bioen_sp = bioen_sp,
                                        tdc_hind = tdc_hind_bt, 
                                        tdr_hind = tdr_hind_bt,
                                        climate_dir= "data/climate/",
                                        hind_yrs = hind_years,
                                        proj_yrs = 2021:2099,
                                        hind_data_start_yr = 1991,
                                        climate_data_start_yr = 1980,
                                        verbose = TRUE)

# C) GFDL SSP 245 (middle of the road fossil-fueled development scenario)
scene_gfdl_245 <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "245",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# D) GFDL SSP 585 (High emission fossil-fueled development scenario)
scene_gfdl_585 <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "585",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# ---------------------------------------------------------------------------- #
# 2. Run Forecast Simulations
# ---------------------------------------------------------------------------- #
# Run the dynamic simulations using Adams-Bashforth (AB) method

forecast_gfdl_persist <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_126     <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_245     <- rsim.run(scene_gfdl_245, method = "AB", years = all_years)
forecast_gfdl_585     <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)

# ---------------------------------------------------------------------------- #
# 3. Calculate B0
# ---------------------------------------------------------------------------- #
# source("R/goa_bref.R")

# Calculate Unfished Biomass (B0)
bzero_persist <- bzero_func(scene_gfdl_persist, managed_sp, F_equil, hind_years, all_years)
bzero_g126    <- bzero_func(scene_gfdl_126, managed_sp, F_equil, hind_years, all_years)
bzero_g585    <- bzero_func(scene_gfdl_585, managed_sp, F_equil, hind_years, all_years)


# ---------------------------------------------------------------------------- #
# 4. Calculate Reference Points (Optional)
# ---------------------------------------------------------------------------- #
# Calculate Biomass Reference Points #WORK IN PROGRESS
#persist_bref <- goa_bref(bzero_persist)
#g126_bref    <- goa_bref(bzero_g126)
#g585_bref    <- goa_bref(bzero_g585)

# Save Reference Points
# save(bzero_persist, bzero_g126, bzero_g585, file = "data/brp/goa_gfdl_brps.RData")

