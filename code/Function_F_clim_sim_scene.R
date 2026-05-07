#' ---------------------------------------------------------------------------- #
#' AUTHORS: Bia Dias
#' ORIGINAL CODE: George (Andy) Whitehouse
#' AFFILIATIONS: CICOES University of Washington
#' E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#' DATE: 25 February 2026
#' @description
#'
#' code/Function_F_clim_sim_scene.R
#' Purpose: Function to set-up a scenario object.
#' Adapted for the GOA fit model (from GOA_fit_results_59M04par.rds).
#' ---------------------------------------------------------------------------- #
#' @param scene the scenario object to be used.
#' @param ssp Share socioeconomic pathway. The ssp scenario to be used for the projection. Options: \code{"126"}, \code{"245"}, \code{"585"}, or \code{"persist"}.
#' @param cons whether to apply the bioenergetics consumption multiplier (TRUE/FALSE). Default \code{TRUE}.
#' @param resp whether to apply the bioenergetics respiration multiplier (TRUE/FALSE). Default \code{TRUE}.
#' @param buf whether to apply the primary production forcing (bottom up forcing-buf) (TRUE/FALSE). Default \code{TRUE}.
#' @param pcod_rec whether to apply the Pacific cod recruitment forcing (TRUE/FALSE). Default \code{FALSE}.
#' @param pcod_rec_method Character. Which thermal curve to use (Laurel and Rogers 2020). Options: \code{"cauchy"} or \code{"gaussian"}. Default \code{"cauchy"}.
#' @param bioen_sp Character vector. List of species to apply bioenergetics
#' @param tdc_hind Data frame. Hindcast of bottom temperature consumption multipliers. It can be sst or bottom temp depending on the fg niche
#' @param tdr_hind Data frame. Hindcast of bottom temperature respiration multipliers. It can be sst or bottom temp depending on the fg niche
#' @param managed_sp a vector of managed species that I want to change the Fs.
#' @param f_equil F equilibrium for non_managed stocks
#' @param f_zero F zero for the given species in zero_fishing_sp
#' @param f_ref_yrs years to calculate the mean F rate
#' @param f_scenario Character. Options: "mean" (historical average for managed), "equil" (equilibrium for all species), "zero" (off for managed), or "zero_all" (off for ALL species in the model). Default \code{"mean"}. 
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
                             pcod_rec = FALSE,
                             pcod_rec_method = "cauchy",
                             bioen_sp,
                             tdc_hind,
                             tdr_hind,
                             managed_sp,
                             f_equil,
                             f_zero,
                             f_ref_yrs = 2016:2020,
                             f_scenario = "mean",
                             zero_fishing_sp = NULL,
                             climate_dir = "data/climate/",
                             hind_yrs = 1991:2020,
                             proj_yrs = 2021:2099,
                             hind_data_start_yr = 1991,
                             climate_data_start_yr = 1980,
                             verbose = TRUE) {
  ssp <- as.character(ssp) # Ensure ssp is character for consistent handling
  ssp <- match.arg(ssp,
                   choices = c("126", "245", "585", "persist"),
                   several.ok = FALSE) # Validate ssp input
  all_years <- c(hind_yrs, proj_yrs)
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
  # Bottom-up forcing (buf) - primary production forcing ####
  # Anomaly reference: monthly means from ssp126 over hind_yrs (1991:2020),
  # matching the reference used in the fitting scene setup.
  # Source file: Long_WGOA_NPZ_PP_monthly_B_added.csv (covers 1980:2099, all SSPs).
  # Change to the function 2026/05/05
  # -------------------------------------------------------------------------- #
  
  if (buf == FALSE) {
    if (verbose)
      message("Not applying primary production forcing")
    
    scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] <- 1
    scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] <- 1
  } else{
    pp_file <- "Rpath_fitting/GOA/wgoa_data_rpath_fitting/Long_WGOA_NPZ_PP_monthly_B_added.csv"
    if(!file.exists(pp_file))
      stop(sprintf("Primary production file not found: %s", pp_file))
    pp_raw <- read.csv(pp_file)
    
    # Hindcast monthly means (ssp126, hind_yrs) — anomaly reference for all SSPs
    pp_ref <- pp_raw %>%
      dplyr::filter(simulation == "ssp126", year %in% hind_yrs) %>% 
      dplyr::group_by(varname, month) %>% 
      dplyr::summarise(Pmean =mean(P_tkm2_mon, na.rm=TRUE), .groups = "drop")
    
    #Helper: compute monthly P_anom series for one varname and SSP, filtered to
    #the requestd years, using the hind_yrs ssp126 monthly mean as the reference
    
    pp_anom_series <- function(varname_sel, years_sel, ssp_sel){
      pp_raw %>% 
        dplyr::filter(simulation == paste0("ssp", ssp_sel),
                      varname    == varname_sel,
                      year      %in% years_sel) %>% 
        dplyr::arrange(year, month) %>% 
        dplyr::left_join(
          dplyr::filter(pp_ref, varname == varname_sel),
                   by = c("varname", "month")
                  ) %>% 
        dplyr::mutate(P_anom = P_tkm2_mon/Pmean) %>% 
        dplyr::pull(P_anom)
    }
    
    
    if (ssp == "persist") {
      if (verbose)
        message("buf: persist — repeating 1991-2020 monthly PP climatology for projection")

      # Repeat the 12-month seasonal cycle of hindcast means across projection years.
      # By construction mean(P_anom) = 1.0 per month, but the seasonal shape is preserved.
      for (grp in list(c("prod_PhL", "large_phytoplankton"),
                       c("prod_PhS", "small_phytoplankton"))) {
        monthly_clim <- pp_anom_series(grp[1], hind_yrs, "126")  %>% 
          matrix(ncol = 12, byrow = TRUE) %>%
          colMeans()                               # 12 monthly climatological means
        scene$forcing$ForcedSearch[ts_projection, grp[2]] <-
          rep(monthly_clim, length(proj_yrs))      # repeat seasonal cycle
      }

    } else {
      if (verbose)
        message(sprintf("buf: applying SSP %s primary production forcing", ssp))

      scene$forcing$ForcedSearch[ts_projection, "large_phytoplankton"] <-
        pp_anom_series("prod_PhL", proj_yrs, ssp)
      scene$forcing$ForcedSearch[ts_projection, "small_phytoplankton"] <-
        pp_anom_series("prod_PhS", proj_yrs, ssp)
    }

    # Zero trap — prevent phytoplankton forcing from hitting exactly zero
    epsilon <- 1e-15
    for (grp in c("large_phytoplankton", "small_phytoplankton")) {
      v <- scene$forcing$ForcedSearch[ts_projection, grp]
      scene$forcing$ForcedSearch[ts_projection, grp] <- pmax(v, epsilon)
    }
  }
  
  # ---------------------------------------------------------------------- #
  # Bioen forcing - consumption and respiration multipliers ####
  # ---------------------------------------------------------------------- #
  bioen_sp_noceph <- bioen_sp[!bioen_sp %in% c("octopus", "squid")]
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
      scene$forcing$ForcedSearch[ts_hindcast, i] <- tdc_hind[data_idx_hind, i]
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
      scene$forcing$ForcedActresp[ts_hindcast, i] <- tdr_hind[data_idx_hind, i]
      
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
# if (pcod_rec) {
#   if (verbose)
#     message(sprintf("Applying Pacific cod recruitment forcing based on temperature proxy using %s method", pcod_rec_method))
#   
#   laurel_rogers_survival_proxy <- function(temp, method="cauchy"){
#     if(method == "cauchy"){
#       return(0.453 / (1 + ((temp - 4.192) / 2.125)^2)) #change to 1 later
#       # Note: The parameters for the Cauchy function (location = 4.192, scale = 2.2125)
#     } else if(method == "gaussian"){
#       return(exp(-0.5 *((temp-4.5)/1.5)^2))
#       # Note: The parameters for the Gaussian function (mean = 4.5, sd = 1.5) 
#     } else {
#       stop("Invalid method. Choose 'cauchy' or 'gaussian'.")
#     }
#   }
#   
#   if(ssp=="persist"){
#     temp_ssp <- "126"
#   } else {
#     temp_ssp <- gsub("ssp","", ssp)
#   }
# 
#   # get the raw temperature projection output from ROMSNPZ scenario
#   temp_file <- paste0(climate_dir,"ssp", temp_ssp,"_wide_WGOA_temp_1000",".csv")
#   
#   if (!file.exists(temp_file)) stop(sprintf("Temp data file not found: %s", temp_file))
#   monthly_temps <- read.csv(temp_file, row.names = NULL)
#   
#   total_years <- length(hind_yrs)+length(proj_yrs)
#   total_months <- total_years * 12 #1:1308 here and 133:1440 tstep
#
#   temp_offset <- (min(hind_yrs) - climate_data_start_yr) * 12
#   data_idx_all <- (temp_offset + 1):(temp_offset + total_months) # 1:1308 (1991:2099)
#   
#   #climate_proj <- read.csv(climate_file, row.names = NULL)
#   #bot_temp     <- as.matrix(climate_proj$btemp[133:1440])
#   #row.names(bot_temp) <- climate_proj$tstep[133:1440]
#   
#   actual_temps_vector <- monthly_temps$btemp[data_idx_all]
#   
#   if(ssp=="persist"){
#     hind_temps <- actual_temps_vector[1:n_hind_months]
#     monthly_mean_temps <- tapply(hind_temps, rep(1:12, length(hind_yrs)), mean)
#     proj_temps <- rep(monthly_mean_temps, length(proj_yrs))
#     actual_temps_vector[(n_hind_months + 1):total_months] <- proj_temps
#   }
#   
#   
#   
#   force_pattern <- rep(1, total_months) # change to 1
#    total_years <- total_months / 12
#   
#   for (yr in 1:total_years) {
#     for (m in 1:12) {
#       month_idx <- (yr - 1) * 12 + m
#       if (m %in% c(2, 3, 4)) {
#         current_temp <- actual_temps_vector[month_idx]
#         force_pattern[month_idx] <- laurel_rogers_survival_proxy(current_temp, method = pcod_rec_method)
#       }
#     }
#   }
#   
#   
#   # Apply to scenario object
#   scene$forcing$ForcedRecs[, "pacific_cod_adult"] <- force_pattern
# }
  
  # -------------------------------------------------------------------------- #
  # Pacific Cod Recruitment Forcing (Temperature Proxy) ####
  # -------------------------------------------------------------------------- #
  if (pcod_rec) {
    if (verbose)
      message(sprintf("Applying Pacific cod recruitment forcing based on temperature proxy using %s method", pcod_rec_method))
    
    laurel_rogers_survival_proxy <- function(temp, method="cauchy"){
      if(method == "cauchy"){
        return(0.453 / (1 + ((temp - 4.192) / 2.125)^2)) #change to 1 later
        # Note: The parameters for the Cauchy function (location = 4.192, scale = 2.2125)
      } else if(method == "gaussian"){
        return(exp(-0.5 *((temp-4.5)/1.5)^2))
        # Note: The parameters for the Gaussian function (mean = 4.5, sd = 1.5) 
      } else {
        stop("Invalid method. Choose 'cauchy' or 'gaussian'.")
      }
    }
    
    if(ssp=="persist"){
      temp_ssp <- "126"
    } else {
      temp_ssp <- gsub("ssp","", ssp)
    }
    
    # get the raw temperature projection output from ROMSNPZ scenario
    temp_file <- paste0(climate_dir,"ssp", temp_ssp,"_wide_WGOA_temp_1000",".csv")
    
    if (!file.exists(temp_file)) stop(sprintf("Temp data file not found: %s", temp_file))
    monthly_temps <- read.csv(temp_file, row.names = NULL)
    
    total_years <- length(hind_yrs)+length(proj_yrs)
    total_months <- total_years * 12 #1:1308 here and 133:1440 tstep
    
    temp_offset <- (min(hind_yrs) - climate_data_start_yr) * 12
    data_idx_all <- (temp_offset + 1):(temp_offset + total_months) # 1:1308 (1991:2099)
    
    btemp_vec <- monthly_temps$btemp[data_idx_all] 
    
    # Compute annual Feb-Apr mean temperature per year #
    # Returns one value per year (length = total_years)
    feb_apr_idx <- function(year_idx) {
      base <- (year_idx - 1) * 12
      base + c(2, 3, 4)  # February, March, April
    }
    annual_FebApr_temp <- sapply(seq_len(total_years), function(yr) {
      mean(btemp_vec[feb_apr_idx(yr)], na.rm = TRUE)
    })
    
    # Reference: hindcast climatological Feb-Apr mean #
    # This is the temperature at which the multiplier = 1.0 by construction.
    n_hind_years   <- length(hind_yrs)
    ref_FebApr_T   <- mean(annual_FebApr_temp[1:n_hind_years], na.rm = TRUE)
    ref_proxy      <- laurel_rogers_survival_proxy(ref_FebApr_T,
                                                   method = pcod_rec_method)
    
    if (verbose) {
      message(sprintf(
        "  Hindcast Feb-Apr btemp climatology: %.3f C (proxy = %.4f)",
        ref_FebApr_T, ref_proxy))
    }
    
    # Build annual multiplier (normalized) — projection years only #
    annual_proxy      <- laurel_rogers_survival_proxy(annual_FebApr_temp,
                                                      method = pcod_rec_method)
    annual_multiplier <- annual_proxy / ref_proxy

    # For 'persist': projection years hold at climatological mean (= 1.0) #
    if (ssp == "persist") {
      annual_multiplier[(n_hind_years + 1):total_years] <- 1.0
    }

    # Expand projection years only to monthly #
    proj_multiplier <- annual_multiplier[(n_hind_years + 1):total_years]
    force_pattern   <- rep(proj_multiplier, each = 12)

    # Diagnostics #
    if (verbose) {
      message(sprintf(
        "  Projection multiplier range: [%.3f, %.3f], mean: %.3f",
        min(force_pattern), max(force_pattern), mean(force_pattern)))
      message(sprintf(
        "  End-of-century (last 20 yrs) mean multiplier: %.3f",
        mean(tail(force_pattern, 20 * 12))))
    }

    # Apply to projection rows only; hindcast keeps calibrated fitted values #
    scene$forcing$ForcedRecs[ts_projection, "pacific_cod_adult"] <- force_pattern
  }
  
   
  
  # -------------------------------------------------------------------------- #
  # Fishing scenarios ####
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
  
  if (verbose) {
    message("Reference F (catch/biomass) by managed species, ",
            min(f_ref_yrs), "-", max(f_ref_yrs), ":")
    print(round(f_mean_ref, 3))
    message("For comparison, F_equil (Ecopath baseline) for these species:")
    print(round(f_equil[managed_sp], 3))
  }
  
  max_rows <- nrow(scene$fishing$ForcedFRate)
  if(max_rows<500){
    target_indices <- which(all_years %in% proj_yrs) # yearly rows
  } else {
    target_indices <- ts_projection #monthly rows
  }
  
  valid_equil_names <- intersect(names(f_equil), colnames(scene$fishing$ForcedFRate))
  background_sp <- setdiff(valid_equil_names, managed_sp)
  all_sp <- colnames(scene$fishing$ForcedFRate)
  
  
  
# if(f_scenario =="zero_all"){
#   scene$fishing$ForcedFRate[target_indices, all_sp] <- 0
# } else {
#   if(length(background_sp) > 0) {
#   scene$fishing$ForcedFRate[target_indices, background_sp] <- matrix(
#     f_equil[background_sp], 
#     nrow = length(target_indices), 
#     ncol = length(background_sp), 
#     byrow = TRUE
#   )
# }
# 
#   if(f_scenario == "mean") {
#     scene$fishing$ForcedFRate[target_indices, managed_sp] <- matrix(
#       f_mean_ref[managed_sp], 
#       nrow = length(target_indices), 
#       ncol = length(managed_sp), 
#       byrow = TRUE
#     )
#   } else if(f_scenario == "zero") {
#     scene$fishing$ForcedFRate[target_indices, managed_sp] <- 0
#   } else {
#     stop("Invalid f_scenario. Choose 'mean', 'zero', or 'zero_all'.")
#   }
#   
#   if (!is.null(zero_fishing_sp)) {
#     scene$fishing$ForcedFRate[target_indices, zero_fishing_sp] <- matrix(
#       f_zero[zero_fishing_sp],
#       nrow = length(target_indices),
#       ncol = length(zero_fishing_sp),
#       byrow = TRUE
#     )
#   }
# }

  # ---- ALWAYS apply F_equil to background species in projection ----------- #
  # This is the symmetric baseline that BOTH F0 and Fmean share. Without this,
  # the two scenarios differ in their treatment of the WHOLE food web, not
  # just managed species, and you can't isolate the effect of managed-species
  # fishing.
  
  if (length(background_sp) > 0) {
    scene$fishing$ForcedFRate[target_indices, background_sp] <- matrix(
      f_equil[background_sp],
      nrow  = length(target_indices),
      ncol  = length(background_sp),
      byrow = TRUE
    )
  }
  
  # ---- Apply the chosen managed-species F scenario ------------------------- #
  if (f_scenario == "mean") {
    if (verbose) message("Applying f_mean_ref to managed species.")
    scene$fishing$ForcedFRate[target_indices, managed_sp] <- matrix(
      f_mean_ref[managed_sp],
      nrow  = length(target_indices),
      ncol  = length(managed_sp),
      byrow = TRUE
    )
    
  } else if (f_scenario == "zero") {
    if (verbose) message("Applying F = 0 to managed species (background unchanged).")
    scene$fishing$ForcedFRate[target_indices, managed_sp] <- 0
    
  } else if (f_scenario == "equil") {
    if (verbose) message("Applying F_equil to managed species (status-quo Ecopath).")
    scene$fishing$ForcedFRate[target_indices, managed_sp] <- matrix(
      f_equil[managed_sp],
      nrow  = length(target_indices),
      ncol  = length(managed_sp),
      byrow = TRUE
    )
    
  } else if (f_scenario == "zero_all") {
    #  Use only when you
    # actually want to release fishing on every species in the model.
    if (verbose) message("Applying F = 0 to ALL species (whole-ecosystem release).")
    all_sp <- colnames(scene$fishing$ForcedFRate)
    scene$fishing$ForcedFRate[target_indices, all_sp] <- 0
    
  } else {
    stop("Invalid f_scenario. Choose 'mean', 'zero', 'equil', or 'zero_all'.")
  }
  
  # ---- Optional override list --------------------------------------------- #
  if (!is.null(zero_fishing_sp)) {
    scene$fishing$ForcedFRate[target_indices, zero_fishing_sp] <- 0
  }
 
  if (verbose)
   message("All done!")
 return(scene) 

  }  
  
  



  
  
  