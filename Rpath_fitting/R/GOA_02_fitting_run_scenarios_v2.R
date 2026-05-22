# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 26 January 2026, updated 21 April 2026
#
# R/GOA_02_fitting_setup.R
# Purpose: Script for setting up the fitting scenarios for the GOA model.
# 
# Rpath_fitting/R/GOA_fitting_setup_a_la_carte.r to set up base scenario objects.
# This script was modified to add the search of vv sensitive groups using the 
# sensitivity analysis results. I noticed that the base model and Bioen had different most sensitive groups organization
# that is due to the constraints of the Bioenergetics parameters that are guiding the fit,
# therefore the optim does not require the same search of VV sensitive groups as the base model.
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #


#source("Rpath_fitting/R/GOA_00_fitting_setup_parameters.R")
source("Rpath_fitting/R/GOA_fitting_setup_a_la_carte.R")
source("Rpath_fitting/R/GOA_01_fitting_optimization.R")
source("Rpath_fitting/R/vv_fit_loop.R")
source("Rpath_fitting/R/fitting_aic.R")

scenarios_to_run <- list(
  #"Base" = scene_base,
  #"PrimProd" = scene_primprod,
  #"Bioen" = scene_bioen,
  #"Full" = scene_full,
  # New: same scenarios with cod recruitment forcing added
  "Base_pcod"     = scene_base_pcod,
  "PrimProd_pcod" = scene_primprod_pcod,
  "Bioen_pcod"    = scene_bioen_pcod,
  "Full_pcod"     = scene_full_pcod
)

#scenarios_to_run <- list("Bioen" = scene_bioen)

# ---------------------------------------------------------------------------- #
# PHASE 1: DYNAMIC SENSITIVITY SEARCH
# This replaces your static species lists
# ---------------------------------------------------------------------------- #
message("Running sensitivity searches to identify top groups per scenario...")
scen_sensitivities <- list()

for (scen_name in names(scenarios_to_run)) {
  curr_scen <- scenarios_to_run[[scen_name]]
  # Get baseline likelihood for this specific scenario
  r_base <- rsim.run(curr_scen, method = "AB", years = hind_years)
  b_like <- rsim.fit.table(curr_scen, r_base)
  
  # Identify sensitivities
  vv_loop_res <- vv_fit_loop(bal, curr_scen, hind_years, run.base.nll, b_like)
  # Store the summary table (sorted by dnll)
  scen_sensitivities[[scen_name]] <- vv_loop_res[[1]][order(vv_loop_res[[1]]$dnll), ]
}


# ---------------------------------------------------------------------------- #
# HELPER: Generate Parameter Vectors on the fly
# ---------------------------------------------------------------------------- #
#get_top_params <- function(vv_table, mzero_list = NULL, total_df = 63) {
#  n_m0 <- length(mzero_list)
#  # Formula: n_m0 + (2 * n_vv) - 1 = total_df
#  n_vv_species <- floor((total_df - n_m0 + 1) / 2)
#
#  # Select top species from sensitivity table
#  top_groups <- head(row.names(vv_table), n_vv_species)
#
#  # Ensure 'offal' is the one with only preyvul (put at the end)
#  if ("offal" %in% top_groups) {
#    top_groups <- c(setdiff(top_groups, "offal"), "offal")
#  }
#
#  spec <- c(mzero_list, top_groups, top_groups[1:(length(top_groups)-1)])
#  type <- c(rep("mzero", n_m0), rep("preyvul", length(top_groups)),
#            rep("predvul", length(top_groups)-1))
#
#  return(list(species = spec, vartype = type))
#}

get_top_params <- function(vv_table,
                           bal,
                           mzero_list = NULL,
                           total_df = 63) {
  # M0 parameters (1 df each)
  current_df <- length(mzero_list)
  
  species_prey <- c()
  species_pred <- c()
  
  # Iterate through ranked sensitive species
  # vv_table should be sorted by dnll (most negative first)
  ranked_species <- row.names(vv_table)
  
  for (sp in ranked_species) {
    # This check is to assign predvul and prey vul
    # Type 0 = Consumer (2 params), Type 1/2 = Producer/Detritus (1 param)
    sp_type <- bal$type[sp]
    params_needed <- ifelse(sp_type == 0, 2, 1)
    
    # Check to see if we still have space on our 63df budget
    if ((current_df + params_needed) <= total_df) {
      species_prey <- c(species_prey, sp)
      
      # predvul to consumers only (Type 0)
      if (sp_type == 0) {
        species_pred <- c(species_pred, sp)
      }
      
      current_df <- current_df + params_needed
    } else {
      # No more adding parameters, we hit our ceiling
      break
    }
  }
  
  # M0s first
  final_species <- mzero_list
  final_vartype <- rep("mzero", length(mzero_list))
  
  # Add all selected PreyVuls
  final_species <- c(final_species, species_prey)
  final_vartype <- c(final_vartype, rep("preyvul", length(species_prey)))
  
  # Add all applicable PredVuls
  final_species <- c(final_species, species_pred)
  final_vartype <- c(final_vartype, rep("predvul", length(species_pred)))
  
  return(list(
    species = final_species,
    vartype = final_vartype,
    df_used = current_df
  ))
}



# ---------------------------------------------------------------------------- #
# Run 0: Baseline (non-fitted) ####
# ---------------------------------------------------------------------------- #

results_nofit <- list()

message("Running baseline (non-fitted) scenarios")

for (scen_name in names(scenarios_to_run)) {
  current_scene <- scenarios_to_run[[scen_name]]
  
  run_base <- rsim.fit.run(
    NA,
    NA,
    NA,
    scene = current_scene,
    run_method = "AB",
    run_years = hind_years,
    verbose = TRUE
  )
  
  
  nll_base <- rsim.fit.run(
    NA,
    NA,
    NA,
    scene = current_scene,
    run_method = "AB",
    run_years = hind_years,
    verbose = FALSE
  )
  
  results_nofit[[scen_name]] <- list(
    scen_name = paste0(scen_name, "_NoFit"),
    opt_object = NULL,
    # No optimization object
    new_scene = current_scene,
    # Scene hasn't changed
    final_run = run_base,
    nll = nll_base,
    aic = rfit_aic(current_scene, NA, nll_base),
    # Ensure rfit_aic handles NA vars
    like_table = rsim.fit.table(current_scene, run_base),
    time = c(0, 0, 0)            # Zero time taken
  )
}

saveRDS(results_nofit, file = "Rpath_fitting/GOA/GOA_fit_results_nofit_v4.rds")

# ---------------------------------------------------------------------------- #
# Run 1: 63 Parameters Fitting (Vulnerabilities) ####
# ---------------------------------------------------------------------------- #

fit_results_63par_codrec <- list()
for (scen_name in names(scenarios_to_run)) {
  p <- get_top_params(scen_sensitivities[[scen_name]], bal = bal, mzero_list = NULL)
  
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_63par"),
    species_vec = p$species,
    vartype_vec = p$vartype,
    start_vec = rep(0, length(p$species)),
    hind_years = hind_years,
    use_parallel = TRUE,
    penalty_weight_vul = 10,
    nll_details = FALSE
  )
  
  fit_results_63par_codrec[[scen_name]] <- fit_result
}

# Save results

#saveRDS(fit_results_63par, file = "GOA/GOA_fit_results_63par.rds")

#saveRDS(fit_results_63par, file = "Rpath_fitting/GOA/GOA_fit_results_63par_v3.rds")
saveRDS(fit_results_63par_codrec, file = "Rpath_fitting/GOA/GOA_fit_results_63par_codrec_v4.rds")

#results_63par <- readRDS("Rpath_fitting/GOA/GOA_fit_results_63par.rds")
fit_results_63par_codrec_v1 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_63par_codrec_v3.rds")
#fit_results_63par_codrec_v4 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_63par_codrec_v4.rds")
# ---------------------------------------------------------------------------- #
# Run 2: 61 Parameters Fitting (Vulnerabilities + M02) ####
# ---------------------------------------------------------------------------- #
mzero_groups_2     <- c("walleye_pollock_adult", "pacific_cod_adult") 

results_61M02par_codrec <- list()
for (scen_name in names(scenarios_to_run)) {
  p <- get_top_params(scen_sensitivities[[scen_name]], bal = bal, mzero_list = mzero_groups_2)
  current_scene <- scenarios_to_run[[scen_name]]
  starts <- c(scenarios_to_run[[scen_name]]$params$MzeroMort[mzero_groups_2], rep(0, length(p$species) -
                                                                                    2))
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_61M02par"),
    species_vec = p$species,
    vartype_vec = p$vartype,
    start_vec = starts,
    hind_years = hind_years,
    penalty_weight_vul = 10,
    use_parallel = TRUE
  )
  results_61M02par_codrec[[scen_name]] <- fit_result
}

# Save results
saveRDS(results_61M02par_codrec, file = "Rpath_fitting/GOA/GOA_fit_results_61M02par_codrec_v4.rds")
#results_61M02par <- readRDS("Rpath_fitting/GOA/GOA_fit_results_61M02par_v3.rds")
results_61M02par_codrec_v1 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_61M02par_codrec_v3.rds")

# ---------------------------------------------------------------------------- #
# Run 3: 59 Parameters Fitting (Vulnerabilities + M03) ####
# ---------------------------------------------------------------------------- #
mzero_groups_3     <- c("walleye_pollock_adult",
                        "pacific_cod_adult",
                        "arrowtooth_flounder_adult")
results_59M03par_codrec <- list()
for (scen_name in names(scenarios_to_run)) {
  p <- get_top_params(scen_sensitivities[[scen_name]], bal = bal, mzero_list = mzero_groups_3)
  starts <- c(scenarios_to_run[[scen_name]]$params$MzeroMort[mzero_groups_3], rep(0, length(p$species) -
                                                                                    3))
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_59M03par"),
    species_vec = p$species,
    vartype_vec = p$vartype,
    start_vec = starts,
    hind_years = hind_years,
    penalty_weight_vul = 10,
    use_parallel = TRUE
  )
  results_59M03par_codrec[[scen_name]] <- fit_result
}
# Save results
saveRDS(results_59M03par_codrec, file = "Rpath_fitting/GOA/GOA_fit_results_59M03par_codrec_v4.rds")

results_59M03par <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M03par_v3.rds")
results_59M03par_codrec_v1 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M03par_codrec_v3.rds")

# ---------------------------------------------------------------------------- #
# Run 4: 59 Parameters Fitting (Vulnerabilities + M04) ####
# Best model ####
# ---------------------------------------------------------------------------- #

mzero_groups_4     <- c(
  "walleye_pollock_adult",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult",
  "deep_water_flatfish"
)
results_59M04par_codrec <- list()
for (scen_name in names(scenarios_to_run)) {
  p <- get_top_params(scen_sensitivities[[scen_name]], bal = bal, mzero_list = mzero_groups_4)
  starts <- c(scenarios_to_run[[scen_name]]$params$MzeroMort[mzero_groups_4], rep(0, length(p$species) -
                                                                                    4))
  
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_59M04par"),
    species_vec = p$species,
    vartype_vec = p$vartype,
    start_vec = starts,
    hind_years = hind_years,
    penalty_weight_vul = 10,
    #penalty_weight_mzero = 10,
    use_parallel = TRUE
  )
  results_59M04par_codrec[[scen_name]] <- fit_result
}
# Save results
saveRDS(results_59M04par_codrec, file = "Rpath_fitting/GOA/GOA_fit_results_59M04par_codrec_v4.rds")
results_59M04par <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par_v3.rds")
results_59M04par_codrec_v1 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par_codrec_v3.rds")
# ---------------------------------------------------------------------------- #
# Run 5: 57 Parameters Fitting (Vulnerabilities + M05) ####
# ---------------------------------------------------------------------------- #



mzero_groups_5     <- c(
  "walleye_pollock_adult",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult",
  "pacific_ocean_perch_adult",
  "deep_water_flatfish"
)


results_57M05par_codrec <- list()
for (scen_name in names(scenarios_to_run)) {
  current_scene <- scenarios_to_run[[scen_name]]
  p <- get_top_params(scen_sensitivities[[scen_name]], bal = bal, mzero_list = mzero_groups_5)
  starts <- c(scenarios_to_run[[scen_name]]$params$MzeroMort[mzero_groups_5], rep(0, length(p$species) -
                                                                                    5))
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_57M05par"),
    species_vec = p$species,
    vartype_vec = p$vartype,
    start_vec = starts,
    hind_years = hind_years,
    penalty_weight_vul = 10,
    use_parallel = TRUE
  )
  results_57M05par_codrec[[scen_name]] <- fit_result
}
# Save results
#saveRDS(results_57M05par_codrec, file = "Rpath_fitting/GOA/GOA_fit_results_57M05par_codrec_v3.rds")

# ---------------------------------------------------------------------------- #
# Run 6: 59 Parameters Fitting (Vulnerabilities + M03 no_pollock) ####
# Best model ####
# ---------------------------------------------------------------------------- #


mzero_groups_3     <- c(
  #"walleye_pollock_adult",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult",
  "deep_water_flatfish"
)
results_59M03par_codrec_noM0poll <- list()
for (scen_name in names(scenarios_to_run)) {
  p <- get_top_params(scen_sensitivities[[scen_name]], bal = bal, mzero_list = mzero_groups_3)
  starts <- c(scenarios_to_run[[scen_name]]$params$MzeroMort[mzero_groups_3], rep(0, length(p$species) -
                                                                                    3))
  
  current_scene <- scenarios_to_run[[scen_name]]
  
  fit_result <- run_rpath_optim(
    scene_obj = current_scene,
    scene_name = paste0(scen_name, "_59M03par_noM0poll"),
    species_vec = p$species,
    vartype_vec = p$vartype,
    start_vec = starts,
    hind_years = hind_years,
    penalty_weight_vul = 10,
    #penalty_weight_mzero = 10,
    use_parallel = TRUE
  )
  results_59M03par_codrec_noM0poll[[scen_name]] <- fit_result
}
# Save results
saveRDS(results_59M03par_codrec_noM0poll, file = "Rpath_fitting/GOA/GOA_fit_results_59M03par_codrec_noM0poll_v4.rds")
results_59M03par_codrec_noM0poll_v1 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M03par_codrec_noM0poll_v3.rds")
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Model Comparison ####
# ---------------------------------------------------------------------------- #

# Helper: parse a model name string into its component dimensions.
# Names follow the pattern: {Interaction}[_pcod]_{ParamSet}
# e.g. "Base_pcod_59M04par", "Bioen_63par", "Full_pcod_59M03par_noM0poll"
parse_model_name <- function(model_name) {
  # Strip trailing _v1 suffix before parsing (added externally to ensure uniqueness)
  nm <- sub("_v1$", "", model_name)
  # Try pcod variant first (more specific)
  m <- regmatches(nm, regexec("^(Base|PrimProd|Bioen|Full)_pcod_(.+)$", nm))[[1]]
  if (length(m) == 3)
    return(data.frame(Interaction = m[2], Has_Cod_Rec = TRUE,
                      Param_Set   = m[3], stringsAsFactors = FALSE))
  # Fallback: no pcod
  m <- regmatches(nm, regexec("^(Base|PrimProd|Bioen|Full)_(.+)$", nm))[[1]]
  if (length(m) == 3)
    return(data.frame(Interaction = m[2], Has_Cod_Rec = FALSE,
                      Param_Set   = m[3], stringsAsFactors = FALSE))
  data.frame(Interaction = NA_character_, Has_Cod_Rec = NA,
             Param_Set   = NA_character_, stringsAsFactors = FALSE)
}

# Extract stats from a fitted model result.
# pollock_ts: number of pollock time series used (1 = current v4, 2 = v1/v3 RDS).
# model_suffix: appended to the Model name to ensure uniqueness across TS versions.
extract_stats <- function(res, pollock_ts = NA_integer_, model_suffix = "") {
  if (is.null(res)) return(NULL)
  n_p1    <- if (!is.null(res$opt_object_s1)) length(res$opt_object_s1$par) else 0L
  conv_s1 <- if (!is.null(res$opt_object_s1)) res$opt_object_s1$convergence == 0 else NA
  n_p2    <- if (!is.null(res$opt_object_s2)) length(res$opt_object_s2$par) else 0L
  conv_s2 <- if (!is.null(res$opt_object_s2)) res$opt_object_s2$convergence == 0 else NA
  model_id <- paste0(res$name, model_suffix)
  parsed   <- parse_model_name(model_id)
  data.frame(
    Model          = model_id,
    Interaction    = parsed$Interaction,
    Has_Cod_Rec    = parsed$Has_Cod_Rec,
    Pollock_TS     = pollock_ts,
    Param_Set      = parsed$Param_Set,
    Num_Parameters = n_p1 + n_p2,
    AICc           = res$aicc[[2]],
    NLL            = res$nll,
    NLLp           = res$nll_penalized,
    Time_Min       = round(res$time[3] / 60, 2),
    converged_s1   = conv_s1,
    converged_s2   = conv_s2,
    stringsAsFactors = FALSE
  )
}

# Extract stats from a non-fitted (baseline) result.
extract_stats_nofit <- function(res, pollock_ts = NA_integer_, model_suffix = "") {
  if (is.null(res)) return(NULL)
  model_id <- paste0(res$scen_name, model_suffix)
  parsed   <- parse_model_name(model_id)
  data.frame(
    Model          = model_id,
    Interaction    = parsed$Interaction,
    Has_Cod_Rec    = parsed$Has_Cod_Rec,
    Pollock_TS     = pollock_ts,
    Param_Set      = "NoFit",
    Num_Parameters = NA_integer_,
    AICc           = res$aic[[2]],
    NLL            = res$nll,
    NLLp           = NA_real_,
    Time_Min       = NA_real_,
    converged_s1   = NA,
    converged_s2   = NA,
    stringsAsFactors = FALSE
  )
}

# --- Extract stats per run block --- #

# No cod recruitment (old runs, no pollock TS distinction)
stats_63par    <- do.call(rbind, lapply(results_63par,    function(r) extract_stats(r)))
stats_61M02par <- do.call(rbind, lapply(results_61M02par, function(r) extract_stats(r)))
stats_59M03par <- do.call(rbind, lapply(results_59M03par, function(r) extract_stats(r)))
stats_59M04par <- do.call(rbind, lapply(results_59M04par, function(r) extract_stats(r)))
#stats_57M05par <- do.call(rbind, lapply(results_57M05par, function(r) extract_stats(r)))

# With cod recruitment – 1 pollock time series shelikof (current v4 runs)
stats_63par_codrec           <- do.call(rbind, lapply(fit_results_63par_codrec,         function(r) extract_stats(r, pollock_ts = 1L)))
stats_61M02par_codrec        <- do.call(rbind, lapply(results_61M02par_codrec,          function(r) extract_stats(r, pollock_ts = 1L)))
stats_59M03par_codrec        <- do.call(rbind, lapply(results_59M03par_codrec,          function(r) extract_stats(r, pollock_ts = 1L)))
stats_59M04par_codrec        <- do.call(rbind, lapply(results_59M04par_codrec,          function(r) extract_stats(r, pollock_ts = 1L)))
stats_59_m03_nopollM0_codrec <- do.call(rbind, lapply(results_59M03par_codrec_noM0poll, function(r) extract_stats(r, pollock_ts = 1L)))
#stats_57M05par_codrec        <- do.call(rbind, lapply(results_57M05par_codrec,          function(r) extract_stats(r, pollock_ts = 1L)))

# With cod recruitment – 2 pollock time series race/shelikof (v1 / v3 RDS)
# Suffix "_v1" is appended to Model names to ensure uniqueness in the table.
stats_63par_codrec_v1           <- do.call(rbind, lapply(fit_results_63par_codrec_v1,         function(r) extract_stats(r, pollock_ts = 2L, model_suffix = "_v1")))
stats_61M02par_codrec_v1        <- do.call(rbind, lapply(results_61M02par_codrec_v1,          function(r) extract_stats(r, pollock_ts = 2L, model_suffix = "_v1")))
stats_59M03par_codrec_v1        <- do.call(rbind, lapply(results_59M03par_codrec_v1,          function(r) extract_stats(r, pollock_ts = 2L, model_suffix = "_v1")))
stats_59M04par_codrec_v1        <- do.call(rbind, lapply(results_59M04par_codrec_v1,          function(r) extract_stats(r, pollock_ts = 2L, model_suffix = "_v1")))
stats_59_m03_nopollM0_codrec_v1 <- do.call(rbind, lapply(results_59M03par_codrec_noM0poll_v1, function(r) extract_stats(r, pollock_ts = 2L, model_suffix = "_v1")))

# Baseline (no fitting) – based on current v4 scenarios (1 pollock TS)
stats_nofit <- do.call(rbind, lapply(results_nofit, function(r) extract_stats_nofit(r, pollock_ts = 1L)))

# --- Build combined AICc table --- #
aicc_table <- rbind(
  stats_nofit,
  # No cod recruitment
  stats_63par,
  stats_61M02par,
  stats_59M03par,
  stats_59M04par,
  #stats_57M05par,
  # Cod recruitment, 1 pollock TS (current v4)
  stats_63par_codrec,
  stats_61M02par_codrec,
  stats_59M03par_codrec,
  stats_59M04par_codrec,
  stats_59_m03_nopollM0_codrec,
  #stats_57M05par_codrec,
  # Cod recruitment, 2 pollock TS (v1 / v3 RDS)
  stats_63par_codrec_v1,
  stats_61M02par_codrec_v1,
  stats_59M03par_codrec_v1,
  stats_59M04par_codrec_v1,
  stats_59_m03_nopollM0_codrec_v1
)

min_aicc              <- min(aicc_table$AICc, na.rm = TRUE)
aicc_table$Delta_AICc <- aicc_table$AICc - min_aicc
aicc_table            <- aicc_table[order(aicc_table$AICc), ]

cat("=== Full model ranking ===\n")
print(aicc_table)

# --- Interaction-level summary (fitted models only) --- #
# Best AICc per combination of Interaction x Has_Cod_Rec x Pollock_TS.
# This lets you directly compare which environmental driver setup wins,
# and whether 1 vs 2 pollock TS makes a consistent difference.
fitted_only <- aicc_table[aicc_table$Param_Set != "NoFit", ]

# Use a string group key to avoid NA being silently dropped by split()
fitted_only$grp_key <- paste(fitted_only$Interaction,
                              fitted_only$Has_Cod_Rec,
                              ifelse(is.na(fitted_only$Pollock_TS), "NoCodRec", fitted_only$Pollock_TS),
                              sep = ".")

interaction_summary <- do.call(rbind, lapply(
  split(fitted_only, fitted_only$grp_key),
  function(grp) {
    if (nrow(grp) == 0) return(NULL)
    best           <- grp[which.min(grp$AICc), ]
    best$N_Configs <- nrow(grp)   # how many param-set configs were evaluated
    best
  }
))
interaction_summary$grp_key    <- NULL   # remove helper column
fitted_only$grp_key            <- NULL
interaction_summary <- interaction_summary[order(interaction_summary$AICc), ]
interaction_summary$Delta_AICc <- interaction_summary$AICc - min_aicc

cat("\n=== Best model per interaction type (fitted only) ===\n")
print(interaction_summary[, c("Interaction", "Has_Cod_Rec", "Pollock_TS", "Param_Set",
                               "Num_Parameters", "AICc", "Delta_AICc", "NLL",
                               "N_Configs", "converged_s1", "converged_s2")])

# --- Save outputs --- #
write.csv(
  aicc_table,
  "Rpath_fitting/GOA/wgoa_data_rpath_fitting/GOA_Model_AIC_Ranking_20260514.csv",
  row.names = FALSE
)
write.csv(
  interaction_summary[, c("Interaction", "Has_Cod_Rec", "Pollock_TS", "Param_Set",
                           "Num_Parameters", "AICc", "Delta_AICc", "NLL",
                           "N_Configs", "converged_s1", "converged_s2")],
  "Rpath_fitting/GOA/wgoa_data_rpath_fitting/GOA_Interaction_Best_AIC_20260514.csv",
  row.names = FALSE
)

# --- Retrieve any model result by name --- #
# For v1 models, strip the "_v1" suffix before looking up in the v1 result lists.
find_model_result <- function(target_name) {
  is_v1    <- grepl("_v1$", target_name)
  base_name <- sub("_v1$", "", target_name)

  pool <- if (is_v1) {
    c(fit_results_63par_codrec_v1,
      results_61M02par_codrec_v1,
      results_59M03par_codrec_v1,
      results_59M04par_codrec_v1,
      results_59M03par_codrec_noM0poll_v1)
  } else {
    c(results_nofit,
      results_63par,
      results_61M02par,
      results_59M03par,
      results_59M04par,
      fit_results_63par_codrec,
      results_61M02par_codrec,
      results_59M03par_codrec,
      results_59M04par_codrec,
      results_59M03par_codrec_noM0poll)
  }

  idx <- which(sapply(pool, function(r) {
    nm <- if (!is.null(r$name)) r$name else r$scen_name
    identical(nm, base_name)
  }))
  if (length(idx) == 0) return(NULL)
  pool[[idx[1]]]
}

# Build a combined likelihood table across all models (work in progress)
combined_like_tables <- do.call(rbind, lapply(aicc_table$Model, function(nm) {
  res <- find_model_result(nm)
  if (!is.null(res) && !is.null(res$like_table)) {
    lt       <- res$like_table
    lt$Model <- nm
    return(lt)
  }
  return(NULL)
}))
