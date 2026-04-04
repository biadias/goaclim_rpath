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

all_years
target_species <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
target_sp_juv <- c("walleye_pollock_juvenile", "pacific_cod_juvenile",
                   "arrowtooth_flounder_juvenile")
predator <- "arrowtooth_flounder_adult"

# ---------------------------------------------------------------------------- #
# Scenarios list ####
# ---------------------------------------------------------------------------- #

scenarios <- list(
  null_F5         = scene_bioen_best,
  persist_F5      = scene_gfdl_persist,
  #bioen             = scene_bioen,
  ssp126_F5          = scene_gfdl_126,
  ssp245_F5          = scene_gfdl_245,
  ssp585_F5          = scene_gfdl_585
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
                                 years=all_years, 
                                 predator=predator)
  
  cat("extracting derivates...\n")
  deriv_info <- get_species_derivatives(scene = current_scene,
                                        years = all_years,
                                        target_species = target_species)
  
  all_results[[scenario_name]] <- list(
    diets = diet_info,
    derivatives = deriv_info
  )
                        
}
#saveRDS(all_results,"results/diagnostics/derivates_59vM04_ssps.rds" )

for(scenario_name in names(scenarios)) {
  cat("Plotting derivates to PDF...\n")
  pdf_path <- file.path("figures/diagnostics", paste0("derivatives_", scenario_name, ".pdf"))
  pdf(file = pdf_path, width = 11, height = 8.5)
  
  current_deriv_info <- all_results[[scenario_name]]$derivatives
  
  
  plot_sp_derivates(output_list=current_deriv_info,
                  target_species=target_species,
                  scenario_name=scenario_name)
  dev.off()
}
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


# ATF juv
my_sp <- "arrowtooth_flounder_juvenile"

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
  file_name <- paste0("figures/diagnostics/ATFjuv_heatmap_deriv_", scenario_name, ".pdf")
  
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


# ---------------------------------------------------------------------------- #
# Consumption (flows) and diet diagnostics - WORK IN PROGRESS####
# ---------------------------------------------------------------------------- #

all_results <- readRDS("results/diagnostics/derivates_59vM04_ssps.rds" )
#null_F5 is just the scene_bioen_best without a rsim run and the other params set up by the F_clim_sim_scene()
flow_mat <- all_results$null_F5$diets$flow
diet_mat <- all_results$null_F5$diets$diet

get_predator_diet_plot(flow_mat, diet_mat, predator = "arrowtooth_flounder_adult", 
                       scenario_name = "null")


