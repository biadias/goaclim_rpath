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
# 1. Run the simulations (this generates the actual output matrices)####
# ---------------------------------------------------------------------------- #
#This is saved in the 02_run_sim.R
all_scenes <- readRDS("results/diagnostics/all_scenes_cons_resp_ssps.rds")

forecast_gfdl_persist      <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_persist_cons <- rsim.run(all_scenes$persist_cons, method = "AB", years = all_years)
forecast_gfdl_persist_resp <- rsim.run(all_scenes$persist_resp, method = "AB", years = all_years)
forecast_gfdl_persist_none <- rsim.run(all_scenes$persist_none, method = "AB", years = all_years)
forecast_gfdl_126          <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_126_cons     <- rsim.run(all_scenes$`126_cons`, method = "AB", years = all_years)
forecast_gfdl_126_resp     <- rsim.run(all_scenes$`126_resp`, method = "AB", years = all_years)
forecast_gfdl_126_none     <- rsim.run(all_scenes$`126_none`, method = "AB", years = all_years)
forecast_gfdl_245          <- rsim.run(scene_gfdl_245, method = "AB", years = all_years)
forecast_gfdl_245_cons     <- rsim.run(all_scenes$`245_cons`, method = "AB", years = all_years)
forecast_gfdl_245_resp     <- rsim.run(all_scenes$`245_resp`, method = "AB", years = all_years)
forecast_gfdl_245_none     <- rsim.run(all_scenes$`245_none`, method = "AB", years = all_years)
forecast_gfdl_585          <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)
forecast_gfdl_585_cons     <- rsim.run(all_scenes$`585_cons`, method = "AB", years = all_years)
forecast_gfdl_585_resp     <- rsim.run(all_scenes$`585_resp`, method = "AB", years = all_years)
forecast_gfdl_585_none     <- rsim.run(all_scenes$`585_none`, method = "AB", years = all_years)


# ---------------------------------------------------------------------------- #
# 2. Plot ssps####
# ---------------------------------------------------------------------------- #
# The names here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios <- list(
  "Persist" = forecast_gfdl_persist,
  "SSP 126" = forecast_gfdl_126,
  "SSP 245" = forecast_gfdl_245,
  "SSP 585" = forecast_gfdl_585
)

#plot_spp <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
plot_arrow_prey <- c(
  "euphausiids",
  "pandalid_shrimp",
  "pacific_capelin",
  "pacific_sandlance",
  "walleye_pollock_adult",
  "walleye_pollock_juvenile",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult",
  "arrowtooth_flounder_juvenile"
)
plot_arrow_pred <- c(
  "steller_sea_lion",
  "longnose_skate",
  "pacific_halibut_adult",
  "arrowtooth_flounder_adult",
  "shallow_water_flatfish",
  "pacific_cod_adult"
)

plot_ground <- c(
  "arrowtooth_flounder_adult",
  "pacific_cod_adult",
  "walleye_pollock_adult")

plot_ground_juv <- c(
  "arrowtooth_flounder_juvenile",
  "pacific_cod_juvenile",
  "walleye_pollock_juvenile")

fg1 <- plot.species[1:12]
fg2 <- plot.species[13:24]
fg3 <- plot.species[25:36]
fg4 <- plot.species[37:48]
fg5 <- plot.species[49:60]
fg6 <- plot.species[61:72]
fg6 <- plot.species[73:84]

vector_list <- list(fg1, fg2, fg3, fg4, fg5, fg6)


for(i in 1:length(vector_list)){
  
  ssp_comps[[i]] <- plot_scenario_comparison(
    sim_list = my_scenarios, 
    species_to_plot = vector_list[[i]], 
    start_year = 1991, 
    variable = "Biomass"  # You can also change this to "Catch"
  )
  print(ssp_comps[[i]])
}

print(ssp_comps)

plot_scenario_comparison(
  sim_list = my_scenarios, 
  species_to_plot = plot_ground_juv, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)
# ---------------------------------------------------------------------------- #
# 2. Plot ssp 126 res con ####
# ---------------------------------------------------------------------------- #

# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios_126 <- list(
  "SSP 126" = forecast_gfdl_126,
  "SSP 126 (consumption only)" = forecast_gfdl_126_cons,
  "SSP 126 (respiration only)" = forecast_gfdl_126_resp,
  "SSP 126 (no bioenergetic modifiers)" = forecast_gfdl_126_none

)

#plot_arrow_prey <- c("euphausiids", "pandalid_shrimp", "pacific_capelin","pacific_sandlance",
#                     "walleye_pollock_adult", "walleye_pollock_juvenile", "pacific_cod_adult",
#                     "arrowtooth_flounder_adult", "arrowtooth_flounder_juvenile")
#plot_arrow_pred <- c("steller_sea_lion",
#                     "longnose_skate",  
#                     "pacific_halibut_adult", "arrowtooth_flounder_adult", "shallow_water_flatfish", 
#                     "pacific_cod_adult"
#)

plot_ground <- c("arrowtooth_flounder_adult", "pacific_cod_adult","walleye_pollock_adult" )
ssp126 <- plot_scenario_comparison(
  sim_list = my_scenarios_126, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(ssp126)

# ---------------------------------------------------------------------------- #
# 2. Plot ssp 245 res con####
# ---------------------------------------------------------------------------- #

# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios_245 <- list(
  "SSP 245" = forecast_gfdl_245,
  "SSP 245 (consumption only)" = forecast_gfdl_245_cons,
  "SSP 245 (respiration only)" = forecast_gfdl_245_resp,
  "SSP 245 (no bioenergetic modifiers)" = forecast_gfdl_245_none
  
)

ssp245 <- plot_scenario_comparison(
  sim_list = my_scenarios_245, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(ssp245)

# ---------------------------------------------------------------------------- #
# 2. Plot ssp 585 res con####
# ---------------------------------------------------------------------------- #

# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios_585 <- list(
  "SSP 585" = forecast_gfdl_585,
  "SSP 585 (consumption only)" = forecast_gfdl_585_cons,
  "SSP 585 (respiration only)" = forecast_gfdl_585_resp,
  "SSP 585 (no bioenergetic modifiers)" = forecast_gfdl_585_none
  
)

ssp585 <- plot_scenario_comparison(
  sim_list = my_scenarios_585, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(ssp585)
