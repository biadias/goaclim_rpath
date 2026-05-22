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

scenarios_to_run <- list("Bioen_pcod" = scene_bioen_pcod)

# ---------------------------------------------------------------------------- #
# PHASE 1: DYNAMIC SENSITIVITY SEARCH
# This replaces your static species lists
# ---------------------------------------------------------------------------- #

sens_file <- "results/diagnostics/vv_scen_sensitivities_Bioen_pcod_59M04.rds"

if (file.exists(sens_file)) {
  scen_sensitivities <- readRDS(sens_file)
} else{
  scen_sensitivities <- list()
  #
  for (scen_name in names(scenarios_to_run)) {
    scenarios_to_run[[scen_name]]$fitting$Biomass <- scenarios_to_run[[scen_name]]$fitting$Biomass %>%
      filter(Year %in% hind_years)
    
    scenarios_to_run[[scen_name]]$fitting$Catch <- scenarios_to_run[[scen_name]]$fitting$Catch %>%
      filter(Year %in% hind_years)
    
    curr_scen <- scenarios_to_run[[scen_name]]
    # Get baseline likelihood for this specific scenario
    r_base <- rsim.run(curr_scen, method = "AB", years = hind_years)
    b_like <- rsim.fit.table(curr_scen, r_base)
    # Identify sensitivities
    vv_loop_res <- vv_fit_loop(bal, curr_scen, hind_years, run.base.nll, b_like)
    # Store the summary table (sorted by dnll)
    scen_sensitivities[[scen_name]] <- vv_loop_res[[1]][order(vv_loop_res[[1]]$dnll), ]
  }
}
#saveRDS(scen_sensitivities, "results/diagnostics/vv_scen_sensitivities_Bioen_pcod_59M04.rds") 


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
result_59_M03 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M03par_codrec_noM0poll_v4.rds")
                         #GOA_fit_results_59M04par_codrec_v4.rds")
#result_59_M04 <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par.rds")

mzero_groups_3     <- c( "pacific_cod_adult",
                        "arrowtooth_flounder_adult", 
                        "deep_water_flatfish" )
p <- get_top_params(scen_sensitivities$Bioen_pcod, bal=bal, mzero_list = mzero_groups_3)



# create scenario object to apply the fit vulnerabilities
combined_pars <- c(
  result_59_M03$Bioen_pcod$opt_object_s2$par, 
  result_59_M03$Bioen_pcod$opt_object_s1$par
)

scene_bioen_best <- scene_bioen_pcod


scene_bioen_best$params <- rsim.fit.apply(
  values = combined_pars, # Uses the fitted vector from the RDS object
  species = p$species,  # Your combined species vector (M0 + prey + pred)
  vartype = p$vartype,  # Your combined vartype vector
  scene.params = scene_bioen_pcod$params
)

#----------------------------------------------------------------------------- #
#Pollock handling time DD ####

#
#pollock_idx <- which(scene_full_pcod$params$spname == "walleye_pollock_adult")
#
#cat("DD length:", length(scene_full_pcod$params$DD),
#    "vs NumPredPreyLinks:", scene_full_pcod$params$NumPredPreyLinks, "\n\n")
#
## DD values for links where pollock is PREY (predators eating pollock)
#prey_links  <- which(scene_full_pcod$params$PreyFrom == pollock_idx)
#pred_names  <- scene_full_pcod$params$spname[scene_full_pcod$params$PreyTo[prey_links]]
#cat("DD — pollock as PREY (predators eating pollock):\n")
#print(data.frame(predator = pred_names,
#                 DD = scene_full_pcod$params$DD[prey_links],
#                 link_idx = prey_links))
#
## DD values for links where pollock is PREDATOR (pollock eating prey)
#pred_links <- which(scene_full_pcod$params$PreyTo == pollock_idx)
#prey_names <- scene_full_pcod$params$spname[scene_full_pcod$params$PreyFrom[pred_links]]
#cat("\nDD — pollock as PREDATOR (pollock eating prey):\n")
#print(data.frame(prey = prey_names,
#                 DD = scene_full_pcod$params$DD[pred_links],
#                 link_idx = pred_links))
#
## Pollock handling params$DD
#all_years  <- 1991:2099
#proj_yrs   <- 2021:2099
#proj_yr_idx <- which(all_years %in% proj_yrs)  # rows 31:109
#cat("Projection year row indices:", range(proj_yr_idx), "\n\n")
#
#scene_f0_base <- scene_bioen_best
#all_gears <- colnames(scene_f0_base$fishing$ForcedFRate)
#scene_f0_base$fishing$ForcedFRate[proj_yr_idx, all_gears] <- 0
#
#cat(
#  "ForcedFRate projection sum (should be 0):",
#  sum(scene_f0_base$fishing$ForcedFRate[proj_yr_idx, ]),
#  "\n\n"
#)
#
#dd_values   <- c(10000, 1000, 100, 10, 2)
#
#dd_f0_runs <- purrr::map(dd_values, \(dd) {
#  s <- scene_f0_base
#  s$params$DD[pred_links] <- dd
#  rsim.run(s, method = "AB", years = all_years)
#})
#names(dd_f0_runs) <- paste0("DD_", dd_values)
#
#purrr::map_dfr(names(dd_f0_runs), \(nm) {
#  tibble(
#    year     = all_years,
#    pollock  = dd_f0_runs[[nm]]$annual_Biomass[, "walleye_pollock_adult"],
#    scenario = nm
#  )
#}) |>
#  filter(year >= 2021) |> 
#  ggplot(aes(year, pollock, colour = scenario)) +
#  geom_line() +
#  geom_hline(yintercept = bal$Biomass["walleye_pollock_adult"],
#             linetype = "dashed",
#             alpha = 0.5) +
#  labs(title = "F0 projection: pollock sensitivity to DD (pollock as predator)", y = "Biomass (t/km²)", x = NULL)
