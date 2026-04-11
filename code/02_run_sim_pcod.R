# ---------------------------------------------------------------------------- #
# 1. Build GFDL-Specific Climate Scenarios ####
# ---------------------------------------------------------------------------- #

# Managed scpecies
managed_sp_list <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")

## Pcod rec  F0####
# Persistence (Base expectation holding recent climate constant)

scene_persist_pcod_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                       f_scenario ="zero_all",
                                       #zero_fishing_sp = managed_sp_list,
                                       f_ref_yrs = 2016:2020,
                                       climate_dir = "data/climate/",
                                       hind_yrs = hind_years,
                                       proj_yrs = 2021:2099,
                                       hind_data_start_yr = 1991,
                                       climate_data_start_yr = 1980,
                                       verbose = TRUE) 


# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_126_pcod_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                   f_scenario ="zero_all",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# C) GFDL SSP 245 (middle of the road fossil-fueled development scenario)
scene_245_pcod_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                   f_scenario ="zero_all",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# D) GFDL SSP 585 (High emission fossil-fueled development scenario)
scene_585_pcod_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                   f_scenario ="zero_all",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

## Scenarios w/o cod rec F0####
scene_persist_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
                                       ssp = "persist",
                                       cons = TRUE,resp = TRUE,buf = FALSE,
                                       pcod_rec= FALSE,
                                       #pcod_rec_method= "cauchy",
                                       bioen_sp = bioen_sp,
                                       tdc_hind = tdc_hind_bt,
                                       tdr_hind = tdr_hind_bt,
                                       managed_sp = managed_sp_list,
                                       f_equil = F_equil,
                                       f_zero = F_zero,
                                       f_scenario ="zero_all",
                                       #zero_fishing_sp = managed_sp_list,
                                       f_ref_yrs = 2016:2020,
                                       climate_dir = "data/climate/",
                                       hind_yrs = hind_years,
                                       proj_yrs = 2021:2099,
                                       hind_data_start_yr = 1991,
                                       climate_data_start_yr = 1980,
                                       verbose = TRUE) 


# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_126_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "126",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= FALSE,
                                   #pcod_rec_method= "cauchy",
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_scenario ="zero_all",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# C) GFDL SSP 245 (middle of the road fossil-fueled development scenario)
scene_245_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "245",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= FALSE,
                                   #pcod_rec_method= "cauchy",
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_scenario ="zero_all",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# D) GFDL SSP 585 (High emission fossil-fueled development scenario)
scene_585_f0 <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "585",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= FALSE,
                                   #pcod_rec_method= "cauchy",
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_scenario ="zero_all",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

## Pcod rec  Fmean 2016-2020####
scene_persist_pcod_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                       f_scenario ="mean",
                                       #zero_fishing_sp = managed_sp_list,
                                       f_ref_yrs = 2016:2020,
                                       climate_dir = "data/climate/",
                                       hind_yrs = hind_years,
                                       proj_yrs = 2021:2099,
                                       hind_data_start_yr = 1991,
                                       climate_data_start_yr = 1980,
                                       verbose = TRUE) 


# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_126_pcod_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                   f_scenario ="mean",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# C) GFDL SSP 245 (middle of the road fossil-fueled development scenario)
scene_245_pcod_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                   f_scenario ="mean",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# D) GFDL SSP 585 (High emission fossil-fueled development scenario)
scene_585_pcod_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
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
                                   f_scenario ="mean",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)





## Scenarios w/o cod rec Fmean####

scene_persist_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
                                  ssp = "persist",
                                  cons = TRUE,resp = TRUE,buf = FALSE,
                                  pcod_rec= FALSE,
                                  #pcod_rec_method= "cauchy",
                                  bioen_sp = bioen_sp,
                                  tdc_hind = tdc_hind_bt,
                                  tdr_hind = tdr_hind_bt,
                                  managed_sp = managed_sp_list,
                                  f_equil = F_equil,
                                  f_zero = F_zero,
                                  f_scenario ="mean",
                                  #zero_fishing_sp = managed_sp_list,
                                  f_ref_yrs = 2016:2020,
                                  climate_dir = "data/climate/",
                                  hind_yrs = hind_years,
                                  proj_yrs = 2021:2099,
                                  hind_data_start_yr = 1991,
                                  climate_data_start_yr = 1980,
                                  verbose = TRUE) 


# B) GFDL SSP 126 (Low emission mitigation scenario)
scene_126_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "126",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= FALSE,
                                   #pcod_rec_method= "cauchy",
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_scenario ="mean",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# C) GFDL SSP 245 (middle of the road fossil-fueled development scenario)
scene_245_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "245",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= FALSE,
                                   #pcod_rec_method= "cauchy",
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_scenario ="mean",
                                   f_ref_yrs = 2016:2020,
                                   climate_dir= "data/climate/",
                                   hind_yrs = hind_years,
                                   proj_yrs = 2021:2099,
                                   hind_data_start_yr = 1991,
                                   climate_data_start_yr = 1980,
                                   verbose = TRUE)

# D) GFDL SSP 585 (High emission fossil-fueled development scenario)
scene_585_fmean <- F_clim_sim_scene(scene = scene_bioen_best,
                                   ssp = "585",
                                   cons = TRUE, resp = TRUE, buf = FALSE,
                                   pcod_rec= FALSE,
                                   #pcod_rec_method= "cauchy",
                                   bioen_sp = bioen_sp,
                                   tdc_hind = tdc_hind_bt, 
                                   tdr_hind = tdr_hind_bt,
                                   managed_sp = managed_sp_list,
                                   f_equil = F_equil,
                                   f_zero = F_zero,
                                   f_scenario ="mean",
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

# Persistence Runs
forecast_persist_pcod_f0    <- rsim.run(scene_persist_pcod_f0, method = "AB", years = all_years)
forecast_persist_f0         <- rsim.run(scene_persist_f0, method = "AB", years = all_years)
forecast_persist_pcod_fmean <- rsim.run(scene_persist_pcod_fmean, method = "AB", years = all_years)
forecast_persist_fmean      <- rsim.run(scene_persist_fmean, method = "AB", years = all_years)

# SSP 126 Runs
forecast_126_pcod_f0    <- rsim.run(scene_126_pcod_f0, method = "AB", years = all_years)
forecast_126_f0         <- rsim.run(scene_126_f0, method = "AB", years = all_years)
forecast_126_pcod_fmean <- rsim.run(scene_126_pcod_fmean, method = "AB", years = all_years)
forecast_126_fmean      <- rsim.run(scene_126_fmean, method = "AB", years = all_years)

# SSP 245 Runs
forecast_245_pcod_f0    <- rsim.run(scene_245_pcod_f0, method = "AB", years = all_years)
forecast_245_f0         <- rsim.run(scene_245_f0, method = "AB", years = all_years)
forecast_245_pcod_fmean <- rsim.run(scene_245_pcod_fmean, method = "AB", years = all_years)
forecast_245_fmean      <- rsim.run(scene_245_fmean, method = "AB", years = all_years)

# SSP 585 Runs
forecast_585_pcod_f0    <- rsim.run(scene_585_pcod_f0, method = "AB", years = all_years)
forecast_585_f0         <- rsim.run(scene_585_f0, method = "AB", years = all_years)
forecast_585_pcod_fmean <- rsim.run(scene_585_pcod_fmean, method = "AB", years = all_years)
forecast_585_fmean      <- rsim.run(scene_585_fmean, method = "AB", years = all_years)

# ---------------------------------------------------------------------------- #
# 3. Compile and Plot by SSP Scenario ####
# ---------------------------------------------------------------------------- #
source("code/Function_plot_scenario_comparison.R")

plot_ground <- c(
  "arrowtooth_flounder_adult",
  "pacific_cod_adult",
  "walleye_pollock_adult",
  "arrowtooth_flounder_juvenile",
  "pacific_cod_juvenile",
  "walleye_pollock_juvenile")

## a) Plot Persistence Scenarios ####
list_persist <- list(
  "F0 w/ Pcod Rec pers"    = forecast_persist_pcod_f0,
  "F0 w/o Pcod Rec"   = forecast_persist_f0,
  "Fmean w/ Pcod Rec" = forecast_persist_pcod_fmean,
  "Fmean w/o Pcod Rec"= forecast_persist_fmean
)
plot_scenario_comparison(sim_list = list_persist, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")

## b) Plot SSP 126 Scenarios ####
list_126 <- list(
  "F0 w/ Pcod Rec 126"    = forecast_126_pcod_f0,
  "F0 w/o Pcod Rec"   = forecast_126_f0,
  "Fmean w/ Pcod Rec" = forecast_126_pcod_fmean,
  "Fmean w/o Pcod Rec"= forecast_126_fmean
)

plot_scenario_comparison(sim_list = list_126, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")

## c) Plot SSP 245 Scenarios ####
list_245 <- list(
  "F0 w/ Pcod Rec 245"    = forecast_245_pcod_f0,
  "F0 w/o Pcod Rec"   = forecast_245_f0,
  "Fmean w/ Pcod Rec" = forecast_245_pcod_fmean,
  "Fmean w/o Pcod Rec"= forecast_245_fmean
)
plot_scenario_comparison(sim_list = list_245, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")

# 4. Plot SSP 585 Scenarios
list_585 <- list(
  "F0 w/ Pcod Rec 585"    = forecast_585_pcod_f0,
  "F0 w/o Pcod Rec"   = forecast_585_f0,
  "Fmean w/ Pcod Rec" = forecast_585_pcod_fmean,
  "Fmean w/o Pcod Rec"= forecast_585_fmean
)
plot_scenario_comparison(sim_list = list_585, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")


# 5. Plot SSPs Scenarios
list_ssps <- list(
  "F0 w/ Pcod Rec 126"    = forecast_126_pcod_f0,
  "F0 w/ Pcod Rec 245"    = forecast_245_pcod_f0,
  "F0 w/ Pcod Rec 585"    = forecast_585_pcod_f0
  #"Fmean w/ Pcod Rec 126" = forecast_126_pcod_fmean,
  #"Fmean w/ Pcod Rec 245" = forecast_245_pcod_fmean,
  #"Fmean w/ Pcod Rec 585" = forecast_585_pcod_fmean
)

plot_scenario_comparison(sim_list = list_ssps, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")


# 5. Plot SSPs Scenarios
list_ssps_nocod <- list(
  "F0 126"    = forecast_126_f0,
  "F0 245"    = forecast_245_f0,
  "F0 585"    = forecast_585_f0,
  "Fmean 126" = forecast_126_fmean,
  "Fmean 245" = forecast_245_fmean,
  "Fmean 585" = forecast_585_fmean
)

plot_scenario_comparison(sim_list = list_ssps_nocod, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")

