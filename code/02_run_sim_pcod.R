# ---------------------------------------------------------------------------- #
# 1. Build GFDL-Specific Climate Scenarios ####
# ---------------------------------------------------------------------------- #

# Managed scpecies
managed_sp_list <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")


# Persistence (Base expectation holding recent climate constant)

scene_gfdl_persist_pcod <- F_clim_sim_scene(scene = scene_bioen_best,
                                       ssp = "persist",
                                       cons = TRUE,resp = TRUE,buf = FALSE,
                                       pcod_rec= TRUE,
                                       pcod_rec_method= "cauchy",
                                       bioen_sp = bioen_sp,
                                       tdc_hind = tdc_hind_bt,
                                       tdr_hind = tdr_hind_bt,
                                       managed_sp = managed_sp_list,
                                       f_equil = F_equil,
                                       f_zero = F_zero,
                                       f_scenario ="zero",
                                       #zero_fishing_sp = managed_sp_list,
                                       f_ref_yrs = 2016:2020,
                                       climate_dir = "data/climate/",
                                       hind_yrs = hind_years,
                                       proj_yrs = 2021:2099,
                                       hind_data_start_yr = 1991,
                                       climate_data_start_yr = 1980,
                                       verbose = TRUE) 


# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_gfdl_126_pcod <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "126",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= TRUE,
                                   pcod_rec_method= "cauchy",
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
scene_gfdl_245_pcod <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "245",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= TRUE,
                                   pcod_rec_method= "cauchy",
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
scene_gfdl_585_pcod <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "585",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= TRUE,
                                   pcod_rec_method= "cauchy",
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

forecast_gfdl_persist_pcod <- rsim.run(scene_gfdl_persist_pcod, method = "AB", years = all_years)
forecast_gfdl_126_pcod     <- rsim.run(scene_gfdl_126_pcod, method = "AB", years = all_years)
forecast_gfdl_245_pcod     <- rsim.run(scene_gfdl_245_pcod, method = "AB", years = all_years)
forecast_gfdl_585_pcod     <- rsim.run(scene_gfdl_585_pcod, method = "AB", years = all_years)

my_scenarios <- list(
  "Persist_pcod" = forecast_gfdl_persist_pcod,
  #"Persist" = forecast_gfdl_persist,
  "SSP 126" = forecast_gfdl_126_pcod,
  "SSP 245" = forecast_gfdl_245_pcod,
  "SSP 585" = forecast_gfdl_585_pcod
)

plot_ground <- c(
  "arrowtooth_flounder_adult",
  "pacific_cod_adult",
  "walleye_pollock_adult",
  "arrowtooth_flounder_juvenile",
  "pacific_cod_juvenile",
  "walleye_pollock_juvenile")
plot_ground_juv <- c(
  "arrowtooth_flounder_juvenile",
  "pacific_cod_juvenile",
  "walleye_pollock_juvenile")


plot_scenario_comparison(
  sim_list = my_scenarios, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)
