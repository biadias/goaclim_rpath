# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 03 March 2026
#
# code/01_load_best_model.R
# Purpose: Loading the best scenario to the set up scene_bioen
# ---------------------------------------------------------------------------- #

# Getting the top VV sensitive groups for the best model
source("code/00_setup_forecast.R")
source("Rpath_fitting/R/vv_fit_loop.R")

#scenarios_to_run <- list(
#  "Base" = scene_base,
#  "PrimProd" = scene_primprod,
#  "Bioen" = scene_bioen,
#  "Full" = scene_full
#)

scenarios_to_run <- list("Bioen" = scene_bioen)

# ---------------------------------------------------------------------------- #
# PHASE 1: DYNAMIC SENSITIVITY SEARCH
# This replaces your static species lists
# ---------------------------------------------------------------------------- #
#scen_sensitivities <- list()
#
#for(scen_name in names(scenarios_to_run)) {
#  scenarios_to_run[[scen_name]]$fitting$Biomass <- scenarios_to_run[[scen_name]]$fitting$Biomass %>%
#    filter(Year %in% hind_years)
#  
#  scenarios_to_run[[scen_name]]$fitting$Catch <- scenarios_to_run[[scen_name]]$fitting$Catch %>%
#    filter(Year %in% hind_years)
#  
#  curr_scen <- scenarios_to_run[[scen_name]]
#  # Get baseline likelihood for this specific scenario
#  r_base <- rsim.run(curr_scen, method = "AB", years = hind_years)
#  b_like <- rsim.fit.table(curr_scen, r_base)
#  # Identify sensitivities
#  vv_loop_res <- vv_fit_loop(bal, curr_scen, hind_years, run.base.nll, b_like)
#  # Store the summary table (sorted by dnll)
#  scen_sensitivities[[scen_name]] <- vv_loop_res[[1]][order(vv_loop_res[[1]]$dnll), ]
#}
#
#saveRDS(scen_sensitivities, "results/diagnostics/vv_scen_sensitivities_Bioen59M04.rds") 

scen_sensitivities <- readRDS("results/diagnostics/vv_scen_sensitivities_Bioen59M04.rds")
# ---------------------------------------------------------------------------- #
# HELPER: Generate Parameter Vectors on the fly
# ---------------------------------------------------------------------------- #


get_top_params <- function(vv_table, bal, mzero_list = NULL, total_df = 63) {
  # 1. Start with M0 parameters (1 df each)
  current_df <- length(mzero_list)
  
  species_prey <- c()
  species_pred <- c()
  
  # 2. Iterate through ranked sensitive species
  # vv_table should be sorted by dnll (most negative first)
  ranked_species <- row.names(vv_table)
  
  for (sp in ranked_species) {
    # Check group type
    # Type 0 = Consumer (2 params), Type 1/2 = Producer/Detritus (1 param)
    sp_type <- bal$type[sp]
    params_needed <- ifelse(sp_type == 0, 2, 1)
    
    # 3. Check if we have room in the 63 df budget
    if ((current_df + params_needed) <= total_df) {
      species_prey <- c(species_prey, sp)
      
      # Only add to predvul list if it's a consumer (Type 0)
      if (sp_type == 0) {
        species_pred <- c(species_pred, sp)
      }
      
      current_df <- current_df + params_needed
    } else {
      # Stop adding once we hit the 63 limit
      break
    }
  }
  
  # 4. Construct final vectors
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
# 1. Load the best model scenario ####
# ---------------------------------------------------------------------------- #
result_59_M04 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par_v3.rds")
mzero_groups_4     <- c("walleye_pollock_adult", "pacific_cod_adult",
                        "arrowtooth_flounder_adult", "deep_water_flatfish" )
p <- get_top_params(scen_sensitivities$Bioen, bal=bal, mzero_list = mzero_groups_4)



# create scenario object to apply the fit vulnerabilities
combined_pars <- c(
  result_59_M04$Bioen$opt_object_s2$par, 
  result_59_M04$Bioen$opt_object_s1$par
)

scene_bioen_best <- scene_bioen


scene_bioen_best$params <- rsim.fit.apply(
  values = combined_pars, # Uses the fitted vector from the RDS object
  species = p$species,  # Your combined species vector (M0 + prey + pred)
  vartype = p$vartype,  # Your combined vartype vector
  scene.params = scene_bioen$params
)

