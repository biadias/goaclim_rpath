# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 26 January 2026
#
# R/GOA_02_fitting_setup.R
# Purpose: Script for setting up the fitting scenarios for the GOA model.
# 
# Rpath_fitting/R/GOA_fitting_setup_a_la_carte.r to set up base scenario objects.
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

source("R/GOA_00_fitting_setup_parameters.R")
source("R/GOA_01_fitting_optimization.R")

scenarios_to_run <- list(
  "Base" = scene_base,
  "PrimProd" = scene_primprod,
  "Bioen" = scene_bioen,
  "Full" = scene_full
)

# ---------------------------------------------------------------------------- #
# Run 0: Baseline (non-fitted) ####
# ---------------------------------------------------------------------------- #

results_nofit <- list()

message("Running baseline (non-fitted) scenarios")

for(scen_name in names(scenarios_to_run)) {
  current_scene <- scenarios_to_run[[scen_name]]
  
  run_base <- rsim.fit.run(NA, NA, NA,
                         scene = current_scene,
                         run_method = "AB",
                         run_years = hind_years,
                         verbose = TRUE)
                         

  nll_base <- rsim.fit.run(NA, NA, NA,
                           scene = current_scene,
                           run_method = "AB",
                           run_years = hind_years,
                           verbose = FALSE)
  
  results_nofit[[scen_name]] <- list(
    scen_name = paste0(scen_name, "_NoFit"),
    opt_object = NULL,         # No optimization object
    new_scene = current_scene, # Scene hasn't changed
    final_run = run_base,
    nll = nll_base,
    aic = rfit_aic(current_scene, NA, nll_base), # Ensure rfit_aic handles NA vars
    like_table = rsim.fit.table(current_scene, run_base),
    time = c(0,0,0)            # Zero time taken
  )
}


# ---------------------------------------------------------------------------- #
# Run 1: 63 Parameters Fitting (Vulnerabilities) ####
# ---------------------------------------------------------------------------- #

fit_results_63par <- list()
for (scen_name in names(scenarios_to_run)) {
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_63par"),
    species_vec = par63_species,
    vartype_vec = par63_vartype,
    start_vec = par63_start,
    hind_years = hind_years,
    use_parallel = TRUE,
    penalty_weight = 10,
    nll_details = FALSE
  )
  
  fit_results_63par[[scen_name]] <- fit_result
}

# Save results
saveRDS(fit_results_63par, file = "GOA/GOA_fit_results_63par.rds")

# ---------------------------------------------------------------------------- #
# Run 2: 61 Parameters Fitting (Vulnerabilities + M02) ####
# ---------------------------------------------------------------------------- #

results_61M02par <- list()
for (scen_name in names(scenarios_to_run)) {
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_61M02par"),
    species_vec = par61M02_species,
    vartype_vec = par61M02_vartype,
    start_vec = par61M02_start,
    hind_years = hind_years,
    use_parallel = TRUE
  )
  results_61M02par[[scen_name]] <- fit_result
}

# Save results
saveRDS(results_61M02par, file = "GOA/GOA_fit_results_61M02par.rds")

# ---------------------------------------------------------------------------- #
# Run 3: 59 Parameters Fitting (Vulnerabilities + M03) ####
# ---------------------------------------------------------------------------- #
results_59M03par <- list()
for (scen_name in names(scenarios_to_run)) {
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_59M03par"),
    species_vec = par59M03_species,
    vartype_vec = par59M03_vartype,
    start_vec = par59M03_start,
    hind_years = hind_years,
    use_parallel = TRUE
  )
  results_59M03par[[scen_name]] <- fit_result
}
# Save results
saveRDS(results_59M03par, file = "GOA/GOA_fit_results_59M03par.rds")
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Model Comparison ####
# ---------------------------------------------------------------------------- #
extract_stats <- function(res) {
  if (is.null(res))
    return(NULL)
  data.frame(
    Model = res$name,
    Num_Parameters = length(res$opt_object$par),
    AICc = res$aicc[[4]],
    NLL = res$nll,
    NLLp= res$opt_object$value,
    Time_Min = round(res$time[3] / 60, 2),
    converged = res$opt_object$convergence == 0,
    stringsAsFactors = FALSE
  )
}

extract_stats_nofit <- function(res) {
  if (is.null(res))
    return(NULL)
  data.frame(
    Model = res$scen_name,
    Num_Parameters = NA,
    AICc = res$aic[[4]],
    NLL = res$nll,
    NLLp= NA,
    Time_Min = NA,
    converged = NA,
    stringsAsFactors = FALSE
  )
}
stats_63par <-  do.call(rbind, lapply(fit_results_63par, extract_stats)) 
stats_61M02par <- do.call(rbind, lapply(results_61M02par, extract_stats))
stats_59M03par <- do.call(rbind, lapply(results_59M03par, extract_stats))
stats_nofit <- do.call(rbind, lapply(results_nofit, extract_stats_nofit))


aicc_table <- rbind(
  stats_63par,
  stats_61M02par,
  stats_59M03par,
  stats_nofit
)

min_aicc <- min(aicc_table$AICc, na.rm = TRUE)
aicc_table$Delta_AICc <- aicc_table$AICc - min_aicc
aicc_table <- aicc_table[order(aicc_table$AICc), ]

print(aicc_table)

write.csv(aicc_table, "GOA/wgoa_data_rpath_fitting/GOA_Model_AIC_Ranking_20260129.csv", row.names = FALSE)

# Combine all tables (work in progress)


find_model_result <- function(target_name) {
  all_results <- c(results_nofit, fit_results_63par, results_61M02par, results_59M03par)
  return(all_results[[target_name]]) # Requires names to be assigned to the combined list
}
combined_like_tables <- do.call(rbind, lapply(rownames(aicc_table), function(nm) {
  res <- find_model_result(nm)
  if (!is.null(res)) {
    lt <- res$like_table
    lt$Model <- nm
    return(lt)
  } else {
    return(NULL)
  }
}))
