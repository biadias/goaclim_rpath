# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 18 March 2026
#
# code/04_run_diagnostics.R
# Purpose: Running diagnostics (hindcast and forecast) on the model's output
# ---------------------------------------------------------------------------- #

source("code/03_run_sim_plot.R")
source("code/Function_diagnostics_TEST.R")

# ---------------------------------------------------------------------------- #
# Vulnerabilities ####
# ---------------------------------------------------------------------------- #

my_vv_table <- get_vv_comparison(original_scene = scene_bioen,
                                 fit_results = result_59_M04$Bioen, 
                                 plot = TRUE, 
                                 scene_name = "bioen_base vs Fitted_bioen",
                                 limit_type="log")

#write.csv(my_vv_table, "results/diagnostics/result_59_M04_vv_comparison_bioen.csv", row.names = FALSE)


# ---------------------------------------------------------------------------- #
# Set up global variables ####
# ---------------------------------------------------------------------------- #

years <- all_years
target_species <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
predator <- "arrowtooth_flounder_adult"

# ---------------------------------------------------------------------------- #
# Scenarios list ####
# ---------------------------------------------------------------------------- #

scenarios <- list(
  gfdl_persist      = scene_gfdl_persist,
  #bioen             = scene_bioen,
  gfdl_126          = scene_gfdl_126,
  gfdl_245          = scene_gfdl_245,
  gfdl_585          = scene_gfdl_585
  #gfdl_persist_cons = scene_gfdl_persist_cons,
  #gfdl_persist_res  = scene_gfdl_persist_res
)

# ---------------------------------------------------------------------------- #
# Plot function loop derivs and diet####
# ---------------------------------------------------------------------------- #


all_results <- list()

if(!dir.exists("results/diagnostics")) dir.create("results/diagnostics")
if(!dir.exists("figures/diagnostics")) dir.create("figures/diagnostics")

for(scenario_name in names(scenarios)) {
  cat("Running diagnostics for scenario:", scenario_name, "\n")

  current_scene <- scenarios[[scenario_name]]
  
  cat("extracting diet info...\n")
  diet_info <- get_predator_diet(scene=current_scene, 
                                 years=years, 
                                 predator=predator)
  
  cat("extracting derivates...\n")
  deriv_info <- get_species_derivatives(scene = current_scene,
                                        years = years,
                                        target_species = target_species)
  
  cat("Plotting derivates to PDF...\n")
  pdf_path <- file.path("figures/diagnostics", paste0("derivatives_", scenario_name, ".pdf"))
  pdf(file = pdf_path, width = 11, height = 8.5)
  
  plot_sp_derivates(output_list=deriv_info,
                         target_species=target_species,
                         scenario_name=scenario_name)
  dev.off()
  
  all_results[[scenario_name]] <- list(
    diets = diet_info,
    derivatives = deriv_info
  )
                        
}
#saveRDS(all_results,"results/diagnostics/derivates_59vM04_ssps.rds" )


# ---------------------------------------------------------------------------- #
# Plot Heatmap ####
# ---------------------------------------------------------------------------- #

# ATF
my_sp <- "arrowtooth_flounder_adult"

for (scenario_name in names(all_results)) {
  
  # 1. Extract the specific data dataframe
  species_data <- all_results[[scenario_name]]$derivatives[[my_sp]]
  
  # 2. Generate the plot
  p <- plot_deriv_change_heatmap(
    species_data = species_data, 
    species_name = paste(scenario_name, my_sp)
  )
  
  # 3. Create a dynamic file name
  # Example: "figures/diagnostics/ATF_heatmap_deriv_gfdl_126.pdf"
  file_name <- paste0("figures/diagnostics/ATF_heatmap_deriv_", scenario_name, ".pdf")
  
  # 4. Save the plot
  ggsave(
    filename = file_name, 
    plot = p, 
    width = 11, 
    height = 8.5
  )
  
  cat("Saved:", file_name, "\n")
}


# POL
my_sp <- "walleye_pollock_adult"

for (scenario_name in names(all_results)) {
  
  # 1. Extract the specific data dataframe
  species_data <- all_results[[scenario_name]]$derivatives[[my_sp]]
  
  # 2. Generate the plot
  p <- plot_deriv_change_heatmap(
    species_data = species_data, 
    species_name = paste(scenario_name, my_sp)
  )
  
  # 3. Create a dynamic file name
  # Example: "figures/diagnostics/ATF_heatmap_deriv_gfdl_126.pdf"
  file_name <- paste0("figures/diagnostics/POL_heatmap_deriv_", scenario_name, ".pdf")
  
  # 4. Save the plot
  ggsave(
    filename = file_name, 
    plot = p, 
    width = 11, 
    height = 8.5
  )
  
  cat("Saved:", file_name, "\n")
}

# COD
my_sp <- "pacific_cod_adult"

for (scenario_name in names(all_results)) {
  
  # 1. Extract the specific data dataframe
  species_data <- all_results[[scenario_name]]$derivatives[[my_sp]]
  
  # 2. Generate the plot
  p <- plot_deriv_change_heatmap(
    species_data = species_data, 
    species_name = paste(scenario_name, my_sp)
  )
  
  # 3. Create a dynamic file name
  # Example: "figures/diagnostics/ATF_heatmap_deriv_gfdl_126.pdf"
  file_name <- paste0("figures/diagnostics/COD_heatmap_deriv_", scenario_name, ".pdf")
  
  # 4. Save the plot
  ggsave(
    filename = file_name, 
    plot = p, 
    width = 11, 
    height = 8.5
  )
  
  cat("Saved:", file_name, "\n")
}
