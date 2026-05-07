# ---------------------------------------------------------------------------- #
# 1. Build GFDL-Specific Climate Scenarios ####
# ---------------------------------------------------------------------------- #
source("code/Function_F_clim_sim_scene.R")
# Managed scpecies
managed_sp_list <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult",
                     "flathead_sole_adult","octopus","deep_water_flatfish","sablefish_adult",
                     "shallow_water_flatfish", "rex_sole_adult","pacific_ocean_perch_adult",
                     "slope_rockfish", "demersal_shelf_rockfish", "pelagic_shelf_rockfish",
                     "atka_mackerel", "big_skate", "longnose_skate", "other_skates", 
                     "pacific_halibut_adult", "salmon_shark","pacific_sleeper_shark")
                     
                     
## Pcod rec  F0####
# Persistence (Base expectation holding recent climate constant)


scene_persist_pcod_f0c <- F_clim_sim_scene(scene = scene_bioen_best,
                                          ssp = "persist",
                                          cons = TRUE,resp = TRUE,buf = FALSE,
                                          pcod_rec= FALSE,
                                          pcod_rec_method= "cauchy",
                                          bioen_sp = bioen_sp,
                                          tdc_hind = tdc_hind,
                                          tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                       tdc_hind = tdc_hind,
                                       tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                       tdc_hind = tdc_hind,
                                       tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                  tdc_hind = tdc_hind,
                                  tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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
                                   tdc_hind = tdc_hind, 
                                   tdr_hind = tdr_hind,
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

scenario_grid <- expand.grid(
  ssp        = c("persist", "126", "245", "585"),
  pcod_rec   = c(TRUE, FALSE),
  buf        = c(TRUE, FALSE),
  f_scenario = c("zero_all", "mean"),
  stringsAsFactors = FALSE
)

scenes <- purrr::pmap(scenario_grid, \(ssp, pcod_rec, buf, f_scenario) {
  F_clim_sim_scene(
    scene      = scene_bioen_best,
    ssp        = ssp,
    pcod_rec   = pcod_rec,
    cons = TRUE, resp = TRUE, buf = buf,
    bioen_sp = bioen_sp,
    tdc_hind = tdc_hind, 
    tdr_hind = tdr_hind,
    managed_sp = managed_sp_list,
    pcod_rec_method = if (pcod_rec) "cauchy" else NULL,
    f_scenario = f_scenario,
    f_equil = F_equil,
    f_zero = F_zero,
    f_ref_yrs = 2016:2020,
    climate_dir= "data/climate/",
    hind_yrs = hind_years,
    proj_yrs = 2021:2099,
    hind_data_start_yr = 1991,
    climate_data_start_yr = 1980,
    verbose = TRUE
  )
})
names(scenes) <- paste(scenario_grid$ssp, 
                       ifelse(scenario_grid$pcod_rec, "pcod", "npcod"),
                       ifelse(scenario_grid$buf, "buf", "nbuf"), 
                       scenario_grid$f_scenario, sep = "_")
  


# ---------------------------------------------------------------------------- #
# 2. Run Forecast Simulations ####
# ---------------------------------------------------------------------------- #
# Run the dynamic simulations using Adams-Bashforth (AB) method


# Persistence Runs
#forecast_persist_pcod_f0c    <- rsim.run(scene_persist_pcod_f0c, method = "AB", years = all_years)
fore_pers_buf_pcod_f0    <- rsim.run(scenes$persist_pcod_buf_zero_all , method = "AB", years = all_years)
fore_pers_buf_f0         <- rsim.run(scenes$persist_npcod_buf_zero_all, method = "AB", years = all_years)
fore_pers_buf_pcod_fmean <- rsim.run(scenes$persist_pcod_buf_mean,      method = "AB", years = all_years)
fore_pers_buf_fmean      <- rsim.run(scenes$persist_npcod_buf_mean,     method = "AB", years = all_years)
fore_pers_pcod_f0    <- rsim.run(scenes$persist_pcod_nbuf_zero_all , method = "AB", years = all_years)
fore_pers_f0         <- rsim.run(scenes$persist_npcod_nbuf_zero_all, method = "AB", years = all_years)
fore_pers_pcod_fmean <- rsim.run(scenes$persist_pcod_nbuf_mean,      method = "AB", years = all_years)
fore_pers_fmean      <- rsim.run(scenes$persist_npcod_nbuf_mean,     method = "AB", years = all_years)


# SSP 126 Runs
forecast_126_buf_pcod_f0    <- rsim.run(scenes$`126_pcod_buf_zero_all`,  method = "AB", years = all_years)
forecast_126_buf_f0         <- rsim.run(scenes$`126_npcod_buf_zero_all`, method = "AB", years = all_years)
forecast_126_buf_pcod_fmean <- rsim.run(scenes$`126_pcod_buf_mean`,      method = "AB", years = all_years)
forecast_126_buf_fmean      <- rsim.run(scenes$`126_npcod_buf_mean`,     method = "AB", years = all_years)
forecast_126_pcod_f0    <- rsim.run(scenes$`126_pcod_nbuf_zero_all`,  method = "AB", years = all_years)
forecast_126_f0         <- rsim.run(scenes$`126_npcod_nbuf_zero_all`, method = "AB", years = all_years)
forecast_126_pcod_fmean <- rsim.run(scenes$`126_pcod_nbuf_mean`,      method = "AB", years = all_years)
forecast_126_fmean      <- rsim.run(scenes$`126_npcod_nbuf_mean`,     method = "AB", years = all_years)

# SSP 245 Runs
forecast_245_buf_pcod_f0    <- rsim.run(scenes$`245_pcod_buf_zero_all`,  method = "AB", years = all_years)
forecast_245_buf_f0         <- rsim.run(scenes$`245_npcod_buf_zero_all`, method = "AB", years = all_years)
forecast_245_buf_pcod_fmean <- rsim.run(scenes$`245_pcod_buf_mean`,      method = "AB", years = all_years)
forecast_245_buf_fmean      <- rsim.run(scenes$`245_npcod_buf_mean`,     method = "AB", years = all_years)
forecast_245_pcod_f0    <- rsim.run(scenes$`245_pcod_nbuf_zero_all`,  method = "AB", years = all_years)
forecast_245_f0         <- rsim.run(scenes$`245_npcod_nbuf_zero_all`, method = "AB", years = all_years)
forecast_245_pcod_fmean <- rsim.run(scenes$`245_pcod_nbuf_mean`,      method = "AB", years = all_years)
forecast_245_fmean      <- rsim.run(scenes$`245_npcod_nbuf_mean`,     method = "AB", years = all_years)

# SSP 585 Runs
forecast_585_buf_pcod_f0    <- rsim.run(scenes$`585_pcod_buf_zero_all`,  method = "AB", years = all_years)
forecast_585_buf_f0         <- rsim.run(scenes$`585_npcod_buf_zero_all`, method = "AB", years = all_years)
forecast_585_buf_pcod_fmean <- rsim.run(scenes$`585_pcod_buf_mean`,      method = "AB", years = all_years)
forecast_585_buf_fmean      <- rsim.run(scenes$`585_npcod_buf_mean`,     method = "AB", years = all_years)
forecast_585_pcod_f0    <- rsim.run(scenes$`585_pcod_nbuf_zero_all`,  method = "AB", years = all_years)
forecast_585_f0         <- rsim.run(scenes$`585_npcod_nbuf_zero_all`, method = "AB", years = all_years)
forecast_585_pcod_fmean <- rsim.run(scenes$`585_pcod_nbuf_mean`,      method = "AB", years = all_years)
forecast_585_fmean      <- rsim.run(scenes$`585_npcod_nbuf_mean`,     method = "AB", years = all_years)

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
  #"F0 w/ Pcod Rec persc"    = forecast_persist_pcod_f0c,
  "F0 w/ Pcod Rec persist"    = forecast_persist_pcod_f0,
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
  "F0 w/ Pcod Rec 585"    = forecast_585_pcod_f0,
  "Fmean w/ Pcod Rec 126" = forecast_126_pcod_fmean,
  "Fmean w/ Pcod Rec 245" = forecast_245_pcod_fmean,
  "Fmean w/ Pcod Rec 585" = forecast_585_pcod_fmean
)

plot_scenario_comparison(sim_list = list_ssps, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")


# 6. Plot SSPs NOCOD_rec Scenarios
list_ssps_fmean <- list(
  #"F0 126"    = forecast_126_f0,
  #"F0 245"    = forecast_245_f0,
  #"F0 585"    = forecast_585_f0,
  "Fmean 126" = forecast_126_fmean,
  "Fmean 245" = forecast_245_fmean,
  "Fmean 585" = forecast_585_fmean,
  "Fmean w/ Pcod Rec 126" = forecast_126_pcod_fmean,
  "Fmean w/ Pcod Rec 245" = forecast_245_pcod_fmean,
  "Fmean w/ Pcod Rec 585" = forecast_585_pcod_fmean
)

plot_scenario_comparison(sim_list = list_ssps_fmean, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")

# 7. Plot SSPs F0 Scenarios
list_ssps_f0 <- list(
  "F0 126"    = forecast_126_f0,
  "F0 245"    = forecast_245_f0,
  "F0 585"    = forecast_585_f0,
  "F0 w/ Pcod Rec 126"    = forecast_126_pcod_f0,
  "F0 w/ Pcod Rec 245"    = forecast_245_pcod_f0,
  "F0 w/ Pcod Rec 585"    = forecast_585_pcod_f0
)

plot_scenario_comparison(sim_list = list_ssps_f0, species_to_plot = plot_ground, 
                         start_year = 1991, variable = "Biomass")

#BIA YOU ARE HERE ####
#Question: How can i cap pollock biomass? is there a density dependence function I can use? ####
