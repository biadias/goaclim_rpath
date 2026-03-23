# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 03 March 2026
#
# code/02_run_sim.R
# Purpose: Running the F_clim_sim_scene function to set up the scenarios and then 
# running the simulations for the best model scenario (GOA_fit_results_59M04par.rds) 
# for the GFDL SSP 126, SSP 425, SSP 585 scenarios, as well as a persistence scenario.
# ---------------------------------------------------------------------------- #
source("code/02_run_sim.R") 
source("code/Function_plot_scenario_comparison.R") 
# ---------------------------------------------------------------------------- #
# 1. Run the simulations (this generates the actual output matrices)
# ---------------------------------------------------------------------------- #
#This is saved in the 02_run_sim.R
forecast_gfdl_persist <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_persist_cons <- rsim.run(scene_gfdl_persist_cons, method = "AB", years = all_years)
forecast_gfdl_persist_res <- rsim.run(scene_gfdl_persist_res, method = "AB", years = all_years)
forecast_gfdl_persist_none <- rsim.run(scene_gfdl_persist_none, method = "AB", years = all_years)
forecast_gfdl_126     <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_245     <- rsim.run(scene_gfdl_245, method = "AB", years = all_years)
forecast_gfdl_585     <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)
# ---------------------------------------------------------------------------- #
# 2. Plot 
# ---------------------------------------------------------------------------- #
# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios <- list(
  "Persist" = forecast_gfdl_persist,
  #"Persist (consumption only)" = forecast_gfdl_persist_cons,
  #"Persist (respiration only)" = forecast_gfdl_persist_res,
  #"Persist (no bioenergetic modifiers)" = forecast_gfdl_persist_none
  "SSP 126" = forecast_gfdl_126,
  "SSP 245" = forecast_gfdl_245,
  "SSP 585" = forecast_gfdl_585
)

#plot_spp <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
plot_arrow_prey <- c("euphausiids", "pandalid_shrimp", "pacific_capelin","pacific_sandlance",
                     "walleye_pollock_adult", "walleye_pollock_juvenile", "pacific_cod_adult",
                     "arrowtooth_flounder_adult", "arrowtooth_flounder_juvenile")
plot_arrow_pred <- c("steller_sea_lion",
                     "longnose_skate",  
                      "pacific_halibut_adult", "arrowtooth_flounder_adult", "shallow_water_flatfish", 
                     "pacific_cod_adult"
                     )

plot_ground <- c("arrowtooth_flounder_adult", "pacific_cod_adult","walleye_pollock_adult" )
my_plot <- plot_scenario_comparison(
  sim_list = my_scenarios, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(my_plot)
