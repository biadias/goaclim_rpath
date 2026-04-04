#' ---------------------------------------------------------------------------- #
#' AUTHORS: Bia Dias
#' ORIGINAL CODE: George (Andy) Whitehouse
#' AFFILIATIONS: CICOES University of Washington
#' E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#' DATE: 25 February 2026
#' @description
#'
#' code/Function_F_clim_sim_scene.R
#' Purpose: Function to set-up a scenario object to be used in solving for target F.
#' Adapted for the GOA fit model (from GOA_fit_results_59M04par.rds).
#' ---------------------------------------------------------------------------- #
#' @param scene the scenario object to be used.
#' @param ssp Share socioeconomic pathway. The ssp scenario to be used for the projection. Options: \code{"126"}, \code{"254"}, \code{"585"}, or \code{"persist"}.
#' @param cons whether to apply the bioenergetics consumption multiplier (TRUE/FALSE). Default \code{TRUE}.
#' @param resp whether to apply the bioenergetics respiration multiplier (TRUE/FALSE). Default \code{TRUE}.
#' @param buf whether to apply the primary production forcing (bottom up forcing-buf) (TRUE/FALSE). Default \code{TRUE}.
#' @param pcod_rec whether to apply the Pacific cod recruitment forcing (TRUE/FALSE). Default \code{TRUE}.
#' @param pcod_rec_method Character. Which thermal curve to use (Laurel and Rogers 2020). Options: \code{"cauchy"} or \code{"gaussian"}. Default \code{"cauchy"}.
#' @param bioen_sp Character vector. List of species to apply bioenergetics
#' @param tdc_hind_bt Data frame. Hindcast of bottom temperature consumption multipliers.
#' @param tdr_hind_bt Data frame. Hindcast of bottom temperature respiration multipliers.
#' @param managed_sp a vector of managed species that I want to change the Fs.
#' @param f_equil F equilibrium for non_managed stocks
#' @param f_zero F zero for the given species in Fzero_sp
#' @param f_ref_yrs years to calculate the mean F rate
#' @param zero_fishing_sp vector of species to set Fzero
#' @param climate_dir Character. Path to the csv files.
#' @param hind_yrs Numeric vector. Hindcast period.\code{1991:2020}.
#' @param proj_yrs Numeric vector. Projection/forecast period. \code{2021:2099}.
#' @param hind_data_start_yr Numeric. Starting year of the hindcast data. Default \code{1991}.
#' @param climate_data_start_yr Numeric. Starting year of the climate data. Default \code{1991}.
#' @param verbose Logical. Whether to print progress messages. Default \code{TRUE}.
#'
#' @return Modified rsim scenario object with updated forcing matrices and parameters for the specified SSP and bioenergetic modifiers.
#' @export
#'


F_clim_sim_scene <- function(scene,
                             ssp,
                             cons = TRUE,
                             resp = TRUE,
                             buf = TRUE,
                             pcod_rec = TRUE,
                             pcod_rec_method = "cauchy",
                             bioen_sp,
                             tdc_hind_bt,
                             tdr_hind_bt,
                             managed_sp,
                             f_equil,
                             f_zero,
                             f_ref_yrs = 2016:2020,
                             zero_fishing_sp = NULL,
                             climate_dir = "data/climate",
                             hind_yrs = 1991:2020,
                             proj_yrs = 2021:2099,
                             hind_data_start_yr = 1991,
                             climate_data_start_yr = 1980,
                             verbose = TRUE) {
  ssp <- as.character(ssp) # Ensure ssp is character for consistent handling
  ssp <- match.arg(ssp,
                   choices = c("126", "245", "585", "persist"),
                   several.ok = FALSE) # Validate ssp input
  
  stopifnot(
    "scene has to be provided" = !missing(scene),
    "bioen_sp must be provided" = !missing(bioen_sp),
    "cons, resp, and buf must be logical" = is.logical(cons) &&
      is.logical(resp) && is.logical(buf) && is.logical(pcod_rec)
  )
  
  if (verbose)
    message(sprintf("setting up climate scenario with ssp: %s", ssp))
  
  n_hind_months <- length(hind_yrs) * 12
  n_proj_months <- length(proj_yrs) * 12
  
  ts_hindcast <- 1:n_hind_months
  ts_projection <- (n_hind_months + 1):(n_hind_months + n_proj_months)
  
  ts_mean_ref <- (n_hind_months - 128 + 1):n_hind_months #133 months for 1991-2001, 120 months for 2001-2020, 128 months for 1991-2020 Change it if needed!
  
  #Calculate subsetting indices for the input data files
  #Hind data (tdc/tdr) offset (1991 start, 1991 hind =0 offset) This is just because I am using the full ROMS output with the Historical period (no physical obs data)
  hind_offset <- (min(hind_yrs) - hind_data_start_yr) * 12
  data_idx_hind <- (hind_offset + 1):(hind_offset + n_hind_months) # 1:360 (1991:2020)
  
  proj_offset <- (min(proj_yrs) - climate_data_start_yr) * 12
  data_idx_proj <- (proj_offset + 1):(proj_offset + n_proj_months)     # 493:1440 (2021:2099)
  
  # -------------------------------------------------------------------------- #
  # Bottom-up forcing (buf) - primary production forcing We don't need right now
  # -------------------------------------------------------------------------- #
  
  if (buf == FALSE) {
    if (verbose)
      message("Not applying primary production forcing")
    # If not applying buf, ensure that the primary producer groups are set to 1 in the forcing matrices
    scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] <- 1
    scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] <- 1
  } else{
    if (ssp == "persist") {
      if (verbose)
        message(
          "Applying primary production forcing for 'persist' scenario by averaging the hindcast data."
        ) # change this if different
      scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] <- mean(scene$forcing$ForcedSearch[ts_mean_ref, "large_phytoplankton"], na.rm = TRUE)
      scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] <- mean(scene$forcing$ForcedSearch[ts_mean_ref, "small_phytoplankton"], na.rm = TRUE)
    } else {
      if (verbose)
        message(
          sprintf(
            "Applying primary production forcing for SSP %s scenario using climate data.",
            ssp
          )
        )
      climate_file <- file.path(climate_dir,
                                paste0("ssp", ssp, "_wide_WGOA_B_1000.csv")) # I need to make this file
      
      if (!file.exists(climate_file))
        stop(sprintf("Climate data file not found: %s", climate_file))
      climate_proj <- read.csv(climate_file, row.names = NULL)
      
      scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] <- climate_proj[data_idx_proj, "large_phytoplankton_ano"]
      scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] <- climate_proj[data_idx_proj, "small_phytoplankton_ano"]
    }
    #zero trap, in order for the biomass not fall under 0
    epsilon <- 1e-15
    scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] <- ifelse(scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] < epsilon,
                                                                               epsilon,
                                                                               scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"])
    scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] <- ifelse(scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] < epsilon,
                                                                               epsilon,
                                                                               scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"])
  }
  
  # ---------------------------------------------------------------------- #
  # Bioen forcing - consumption and respiration multipliers
  # ---------------------------------------------------------------------- #
  bioen_sp_noceph <- bioen_sp[!bioen_sp %in% c("octopus", "squids")]
  if (ssp != "persist" && (cons == TRUE || resp == TRUE)) {
    if (verbose)
      message("Applying bioenergetic forcing for SSP scenario using climate data.")
    bioen_proj <- bioen_params(ssp)
  }
  #consumption multiplier
  if (cons) {
    if (verbose)
      message("Applying consumption multipliers.")
    for (i in bioen_sp_noceph) {
      #Hindcast
      scene$forcing$ForcedSearch[ts_hindcast, i] <- tdc_hind_bt[data_idx_hind, i]
      #Projection
      
      if (ssp == "persist") {
        scene$forcing$ForcedSearch[ts_projection, i] <- mean(scene$forcing$ForcedSearch[ts_mean_ref, i])
      } else {
        scene$forcing$ForcedSearch[ts_projection, i] <- bioen_proj[[1]][, i]
      }
    }
    
  }
  #respiration multiplier
  if (resp) {
    if (verbose)
      message("Applying respiration multipliers.")
    for (i in bioen_sp_noceph) {
      #Hindcast
      scene$forcing$ForcedActresp[ts_hindcast, i] <- tdr_hind_bt[data_idx_hind, i]
      
      #Projection
      if (ssp == "persist") {
        scene$forcing$ForcedActresp[ts_projection, i] <- mean(scene$forcing$ForcedActresp[ts_mean_ref, i])
      } else {
        scene$forcing$ForcedActresp[ts_projection, i] <- bioen_proj[[2]][, i]
      }
    }
  }
  if (verbose)
    message("Climate scenario setup complete.")
  
  # -------------------------------------------------------------------------- #
  # Pacific Cod Recruitment Forcing (Temperature Proxy)
  # -------------------------------------------------------------------------- #
  if (pcod_rec) {
    if (verbose)
      message(sprintf("Applying Pacific cod recruitment forcing based on temperature proxy using %s method", pcod_rec_method))
    
    laurel_rogers_survival_proxy <- function(temp, method="cauchy"){
      if(method == "cauchy"){
        return(1 / (1 + ((temp - 4.192) / 2.2125)^2))
        # Note: The parameters for the Cauchy function (location = 4.192, scale = 2.2125)
      } else if(method == "gaussian"){
        return(exp(-0.5 *((temp-4.5)/1.5)^2))
        # Note: The parameters for the Gaussian function (mean = 4.5, sd = 1.5) 
      } else {
        stop("Invalid method. Choose 'cauchy' or 'gaussian'.")
      }
    }
    
    ssp_clean <- gsub("ssp","", ssp)
    # get the raw temperature projection output from ROMSNPZ scenario
    climate_file <- paste0("Rpath_fitting/GOA/wgoa_data_rpath_fitting/ssp", ssp_clean,"_wide_WGOA_temp_1000",".csv")
    climate_proj <- read.csv(climate_file, row.names = NULL)
    bot_temp     <- as.matrix(climate_proj$btemp[493:1440])
    row.names(bot_temp) <- climate_proj$tstep[493:1440]
    
    
    # Load the monthly temperature data for the projection period
    ssp_clean <- gsub("ssp","", ssp)
    temp_file <- file.path(climate_dir, paste0("ssp", ssp_clean, "_wide_WGOA_temp_1000.csv"))
    
    if (!file.exists(temp_file))
      stop(sprintf("Temperature data file not found: %s", temp_file))
    
    temp_proj <- read.csv(temp_file, row.names = NULL)
    temp_proj_vector <- temp_proj[data_idx_proj, "btemp"] # 'btemp' is the bottom temperature
    
    # Apply the recruitment forcing based on the chosen method
    if (pcod_rec_method == "cauchy") {
      rec_multiplier <- 1 / (1 + ((temp_proj_vector - 4.192) / 2.2125)^2)
    } else if (pcod_rec_method == "gaussian") {
      rec_multiplier <- exp(-0.5 * ((temp_proj_vector - 4.5) / 1.5)^2)
    } else {
      stop("Invalid pcod_rec_method. Choose 'cauchy' or 'gaussian'.")
    }
    
    # Set the recruitment forcing in the scenario
    scene$forcing$ForcedRecs[ts_projection, "pacific_cod_adult"] <- rec_multiplier
  }
  
  
  
  # -------------------------------------------------------------------------- #
  # Fishing scenarios
  # -------------------------------------------------------------------------- #
  
  if(verbose)
    message(sprintf(
      "Calculating mean F from %s-%s...",
      min(f_ref_yrs),
      max(f_ref_yrs)
    ))
  run.hind <- rsim.run(scene, method = "AB", years = hind_yrs)
  ref_idx <- which(hind_yrs %in% f_ref_yrs)
  
  #F=Catch (dicards+landings)/Biomass
  catch_ref <- run.hind$annual_Catch[ref_idx, managed_sp, drop = FALSE]
  bio_ref <- run.hind$annual_Biomass[ref_idx, managed_sp, drop = FALSE]
  f_mean_ref <- apply(catch_ref / bio_ref, 2, mean, na.rm = TRUE)
  
  if (verbose)
    message("applying F_equil and scene fishing rates")
  
  max_rows <- nrow(scene$fishing$ForcedFRate)
  if(max_rows<500){
    target_indices <- which(all_years %in% proj_yrs)
  } else {
    target_indices <- ts_projection
  }
  
  valid_equil_names <- intersect(names(f_equil), colnames(scene$fishing$ForcedFRate))
  background_sp <- setdiff(valid_equil_names, managed_sp)
  
  for (m in target_indices) {
    if(length(background_sp)>0){
    #F equilibrium for everybody
    scene$fishing$ForcedFRate[m, background_sp] <- f_equil[background_sp]
    } 
    
    #Mean F for managed species
    scene$fishing$ForcedFRate[m, managed_sp] <- f_mean_ref[managed_sp]
    
    if (!is.null(zero_fishing_sp)) {
      scene$fishing$ForcedFRate[m, zero_fishing_sp] <- f_zero[zero_fishing_sp]
    }
  }
  
  if (verbose)
    message("All done!")
  return(scene)
  
}
