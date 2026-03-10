# 1. Set up your scenarios using the function we just refined
source("code/02_run_sim.R") # Make sure this script runs without errors and that the scene objects are created successfullyn")
source("code/Function_plot_scenario_comparison.R")  # Make sure to source the function if it's not already in your environment

# 2. Run the simulations (this generates the actual output matrices)
forecast_gfdl_persist <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_126     <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_245     <- rsim.run(scene_gfdl_245, method = "AB", years = all_years)
forecast_gfdl_585     <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)

# 3. Create a NAMED LIST of your results. 
# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios <- list(
  "Persist" = forecast_gfdl_persist,
  "SSP 126" = forecast_gfdl_126,
  "SSP 245" = forecast_gfdl_245,
  "SSP 585" = forecast_gfdl_585
)

# 4. Define the species you want to compare
plot_spp <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")

# 5. Generate the plot
my_plot <- plot_scenario_comparison(
  sim_list = my_scenarios, 
  species_to_plot = plot_spp, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

# Show the plot
print(my_plot)
