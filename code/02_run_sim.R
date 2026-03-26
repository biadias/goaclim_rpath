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
#source("code/00_setup_forecast.R")
source("code/01_load_best_model.R") # Loads the best model scenario and sets up scene_bioen
source("Rpath_fitting/GOA/wgoa_bioenergetics_code/bioenergetic_projections.r") # Make sure this is adapted for GOA and matches the years used in the projections
source("code/Function_F_clim_sim_scene.R")




# ---------------------------------------------------------------------------- #
# 1. Build GFDL-Specific Climate Scenarios ####
# ---------------------------------------------------------------------------- #

# Managed scpecies
managed_sp_list <- c("walleye_pollock_adult", "pacific_cod_adult", "sablefish_adult")


# Persistence (Base expectation holding recent climate constant)

scene_gfdl_persist <- F_clim_sim_scene(scene = scene_bioen_best,
                                       ssp = "persist",
                                       cons = TRUE,resp = TRUE,buf = FALSE,
                                       bioen_sp = bioen_sp,
                                       tdc_hind = tdc_hind_bt,
                                       tdr_hind = tdr_hind_bt,
                                       managed_sp = managed_sp_list,
                                       f_equil = F_equil,
                                       f_zero = F_zero,
                                       f_ref_yrs = 2016:2020,
                                       climate_dir = "data/climate/",
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
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_ref_yrs = 2016:2020,
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
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_ref_yrs = 2016:2020,
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
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# ---------------------------------------------------------------------------- #
# 2. Run Forecast Simulations ####
# ---------------------------------------------------------------------------- #
# Run the dynamic simulations using Adams-Bashforth (AB) method

forecast_gfdl_persist <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_126     <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_245     <- rsim.run(scene_gfdl_245, method = "AB", years = all_years)
forecast_gfdl_585     <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)

# ---------------------------------------------------------------------------- #
# 3. Calculate B0 ####
# ---------------------------------------------------------------------------- #
# source("R/goa_bref.R")

# Calculate Unfished Biomass (B0)
bzero_persist <- bzero_func(scene_gfdl_persist, managed_sp, F_equil, hind_years, all_years)
bzero_g126    <- bzero_func(scene_gfdl_126, managed_sp, F_equil, hind_years, all_years)
bzero_g585    <- bzero_func(scene_gfdl_585, managed_sp, F_equil, hind_years, all_years)


# ---------------------------------------------------------------------------- #
# 4. Calculate Reference Points (Optional) ####
# ---------------------------------------------------------------------------- #
# Calculate Biomass Reference Points #WORK IN PROGRESS
#persist_bref <- goa_bref(bzero_persist)
#g126_bref    <- goa_bref(bzero_g126)
#g585_bref    <- goa_bref(bzero_g585)

# Save Reference Points
# save(bzero_persist, bzero_g126, bzero_g585, file = "data/brp/goa_gfdl_brps.RData")


# ---------------------------------------------------------------------------- #
# 5. Test fitted model cons res for sll ssps ####
# ---------------------------------------------------------------------------- #

ssp_scenarios <- c("persist", "126", "245", "585")
bio_modes <- list(
  cons = list(cons = TRUE,  resp = FALSE),
  resp = list(cons = FALSE, resp = TRUE),
  none = list(cons = FALSE, resp = FALSE)
)

all_scenes <- list()
for(s in ssp_scenarios){
  for(m in names(bio_modes)){
    run_name <- paste0(s,"_", m)
    
    all_scenes[[run_name]] <- F_clim_sim_scene(
      scene=scene_bioen_best,
      ssp=s,
      cons=bio_modes[[m]]$cons,
      resp=bio_modes[[m]]$resp,
      buf= FALSE,
      bioen_sp = bioen_sp,
      tdc_hind = tdc_hind_bt, 
      tdr_hind = tdr_hind_bt,
      managed_sp = managed_sp_list,
      f_equil = F_equil,
      f_zero = F_zero,
      f_ref_yrs = 2016:2020,
      climate_dir = "data/climate/",
      hind_yrs = hind_years,
      proj_yrs = 2021:2099,
      hind_data_start_yr = 1991,
      climate_data_start_yr = 1980,
      verbose = FALSE
      
    )
  }
}
#saveRDS(all_scenes, file="results/diagnostics/all_scenes_cons_resp_ssps.rds")

# ---------------------------------------------------------------------------- #
# 6. Topic 2  ####
# ---------------------------------------------------------------------------- #

# Define the species we want to calculate B0 for
managed_sp_list <- c("walleye_pollock_adult", "pacific_cod_adult", "sablefish_adult")

# Create a container for results
b0_results <- data.frame(Species = managed_sp_list, B0_2099 = NA)

for (sp in managed_sp_list) {
  message(paste("Processing B0 for:", sp))
  
  # Create scenario where ONLY 'sp' has F=0
  this_scene <- F_clim_sim_scene(
    scene           = scene_bioen_best, 
    cons = FALSE, resp = FALSE, buf = FALSE,
    ssp             = "126", 
    managed_sp      = managed_sp_list, # All other managed species stay at mean F
    f_equil         = F_equil, 
    f_zero          = F_zero,
    zero_fishing_sp = sp,         # Toggle this specific species to 0
    bioen_sp        = bioen_sp_noceph,
    tdc_hind_bt     = tdc_hind_bt,
    tdr_hind_bt     = tdr_hind_bt,
    hind_yrs        = 1991:2020,         # Good practice to pass these explicitly
    proj_yrs        = 2021:2099,
    verbose         = FALSE              # Set to FALSE to keep console clean during loops
  )
  
  # Run simulation
  this_run <- rsim.run(this_scene, method = "AB", years = all_years)
  
  # Extract mean biomass of the last 10 years
  final_b0 <- mean(tail(this_run$annual_Biomass[, sp], 10))
  b0_results[b0_results$Species == sp, "B0_2099"] <- final_b0
}

print(b0_results)
