# ---------------------------------------------------------------------------- #
# wgoa_add_pcod_rec_to_scene.R
#
# Adds Pacific cod recruitment forcing to a fitting scene using the
# Laurel-Rogers (2020) survival proxy applied to annual Feb-Apr bottom
# temperature. This is applied to scene$forcing$ForcedRecs
#
# Forcing logic:
#   - Annual Feb-Apr mean btemp computed per year
#   - Laurel-Rogers proxy applied to each annual temperature
#   - Multiplier = proxy(year_T) / proxy(hindcast_mean_FebApr_T)
#   - Mean multiplier across hindcast = 1.0 (anomaly forcing)
#   - Single annual multiplier expanded to all 12 months of that year
#
# This is the SAME formula used in F_clim_sim_scene.R for projection, ensuring
# the projection forcing is internally consistent with the fit.
#
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Laurel & Rogers (2020) survival proxy
# ---------------------------------------------------------------------------- #
laurel_rogers_survival_proxy <- function(temp,  method = "cauchy") {
  if (method == "cauchy") {
    return(1 / (1 + ((temp - 4.192) / 2.125)^2)) # since we are on anomaly space I don't need to put the ceiling value.
  } else if (method == "gaussian") {
    return(exp(-0.5 * ((temp - 4.5) / 1.5)^2))
  } else {
    stop("Invalid method. Choose 'cauchy' or 'gaussian'.")
  }
}

# ---------------------------------------------------------------------------- #
# Add cod recruitment forcing to a fitting scene
#
# scene        : Rpath fitting scene (covering hind_years only)
# hind_years   : vector of hindcast years (e.g., 1991:2020)
# temp_file    : path to ROMS hindcast monthly temperature CSV
#                (the same one used for tdc_hind/tdr_hind)
# rec_method   : "cauchy" or "gaussian"
# verbose      : if TRUE, prints multiplier diagnostics
# ---------------------------------------------------------------------------- #
add_pcod_rec_forcing <- function(scene,
                                 hind_years,
                                 temp_file = "Rpath_fitting/GOA/wgoa_data_rpath_fitting/Long_WGOA_temp_monthly_1000.csv",
                                 rec_method = "cauchy",
                                 verbose = TRUE) {
  
  if (verbose)
    message(sprintf("Adding pcod recruitment forcing (%s method)", rec_method))
  
  # ---- Load ROMS hindcast temperatures (ssp126 hindcast) ----------------- #
  if (!file.exists(temp_file))
    stop(sprintf("Temp file not found: %s", temp_file))
  
  roms_hind <- read.csv(temp_file) %>%
    dplyr::filter(simulation == "ssp126",
                  year %in% hind_years) %>%
    dplyr::select(year, month, depthclass, area_weighted_temp) %>%
    tidyr::pivot_wider(names_from = depthclass,
                       values_from = area_weighted_temp) %>%
    dplyr::rename(temp_b5 = Bottom) %>%
    dplyr::arrange(year, month)
  
  expected_rows <- length(hind_years) * 12
  if (nrow(roms_hind) != expected_rows)
    stop(sprintf("ROMS hindcast has %d rows but expected %d (years %d-%d).",
                 nrow(roms_hind), expected_rows,
                 min(hind_years), max(hind_years)))
  
  # ---- Compute annual Feb-Apr mean btemp -------------------------------- #
  feb_apr <- roms_hind %>%
    dplyr::filter(month %in% c(2, 3, 4)) %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(FebApr_btemp = mean(temp_b5, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(year)
  
  if (nrow(feb_apr) != length(hind_years))
    stop("Feb-Apr aggregation produced wrong number of years.")
  
  # ---- Reference: hindcast climatological Feb-Apr mean ------------------ #
  ref_FebApr_T <- mean(feb_apr$FebApr_btemp, na.rm = TRUE)
  ref_proxy    <- laurel_rogers_survival_proxy(ref_FebApr_T, method = rec_method)
  
  # ---- Annual multiplier (anomaly forcing, mean = 1.0) ------------------ #
  annual_proxy      <- laurel_rogers_survival_proxy(feb_apr$FebApr_btemp,
                                                    method = rec_method)
  annual_multiplier <- annual_proxy / ref_proxy
  
  # ---- Expand to monthly (12 months per year) --------------------------- #
  force_pattern <- rep(annual_multiplier, each = 12)
  end_hind      <- length(hind_years) * 12

  # The scene may span beyond hind_years (e.g. into a forecast period).
  # Apply forcing only to hindcast rows; remaining rows are left unchanged.
  if (end_hind > nrow(scene$forcing$ForcedRecs))
    stop(sprintf("Hindcast forcing (%d months) exceeds scene length (%d rows).",
                 end_hind, nrow(scene$forcing$ForcedRecs)))

  # ---- Diagnostics ------------------------------------------------------ #
  if (verbose) {
    message(sprintf("  Hindcast Feb-Apr btemp climatology: %.3f C",
                    ref_FebApr_T))
    message(sprintf("  Multiplier range: [%.3f, %.3f], mean: %.3f",
                    min(annual_multiplier), max(annual_multiplier),
                    mean(annual_multiplier)))
    # Identify warmest and coldest years
    yr_max <- feb_apr$year[which.max(feb_apr$FebApr_btemp)]
    yr_min <- feb_apr$year[which.min(feb_apr$FebApr_btemp)]
    message(sprintf("  Warmest Feb-Apr: %d (%.2f C, mult=%.3f)",
                    yr_max, max(feb_apr$FebApr_btemp),
                    annual_multiplier[which.max(feb_apr$FebApr_btemp)]))
    message(sprintf("  Coldest Feb-Apr: %d (%.2f C, mult=%.3f)",
                    yr_min, min(feb_apr$FebApr_btemp),
                    annual_multiplier[which.min(feb_apr$FebApr_btemp)]))
  }
  
  # ---- Apply to scene (hindcast rows only) ------------------------------ #
  scene$forcing$ForcedRecs[1:end_hind, "pacific_cod_adult"] <- force_pattern
  
  return(scene)
}