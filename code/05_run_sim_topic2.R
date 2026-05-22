# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 03 March 2026
#
# code/05_run_sim_topic2.R
# Purpose: Results (GOA_fit_results_59M03par_codrec_noM0poll_v4.rds) 
# for the GFDL SSP 126, SSP 425, SSP 585 scenarios, as well as a persistence scenario.
# ---------------------------------------------------------------------------- #

#devtools::install_github("NOAA-EDAB/Rpath@ssb_output_diag")
#library(Rpath)
source("code/02_run_sim_pcod.R")

# ---------------------------------------------------------------------------- #
# 1. Model lists ####
# ---------------------------------------------------------------------------- #
# format model outputs
# species: pollock, pcod, and arrowtooth
T2_species      <- c("POL", "COD", "ATF")
T2_adults       <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
T2_full_species <- c("walleye_pollock_adult", "walleye_pollock_juvenile",
                     "pacific_cod_adult",     "pacific_cod_juvenile",
                     "arrowtooth_flounder_adult", "arrowtooth_flounder_juvenile")
# December vector
tdec <- seq(12,1308,12)

list_ssps_pp <- list(
  "F0 persist Pcod rec PP" = fore_pers_buf_pcod_f0,
  "F0 126 Pcod rec PP" = forecast_126_buf_pcod_f0,
  "F0 245 Pcod rec PP" = forecast_245_buf_pcod_f0,
  "F0 585 Pcod rec PP" = forecast_585_buf_pcod_f0,
  "Fmean persist Pcod rec PP" = fore_pers_buf_pcod_fmean,
  "Fmean 126 Pcod rec PP" = forecast_126_buf_pcod_fmean,
  "Fmean 245 Pcod rec PP" = forecast_245_buf_pcod_fmean,
  "Fmean 585 Pcod rec PP" = forecast_585_buf_pcod_fmean
)

list_ssps_pcod <- list(
  "F0 persist Pcod rec" = fore_pers_pcod_f0,
  "F0 126 Pcod rec" = forecast_126_pcod_f0,
  "F0 245 Pcod rec" = forecast_245_pcod_f0,
  "F0 585 Pcod rec" = forecast_585_pcod_f0,
  "Fmean persist Pcod rec" = fore_pers_pcod_fmean,
  "Fmean 126 Pcod rec" = forecast_126_pcod_fmean,
  "Fmean 245 Pcod rec" = forecast_245_pcod_fmean,
  "Fmean 585 Pcod rec" = forecast_585_pcod_fmean
)

list_ssps <- list(
  "F0 persist Pcod" = fore_pers_f0,
  "F0 126 Pcod" = forecast_126_f0,
  "F0 245 Pcod" = forecast_245_f0,
  "F0 585 Pcod" = forecast_585_f0,
  "Fmean persist Pcod" = fore_pers_fmean,
  "Fmean 126 Pcod" = forecast_126_fmean,
  "Fmean 245 Pcod" = forecast_245_fmean,
  "Fmean 585 Pcod" = forecast_585_fmean
)

# ---------------------------------------------------------------------------- #
# 2. Helper function: extract SSB, total biomass, and annual catch ####
# ---------------------------------------------------------------------------- #
extract_model_data <- function(model, model_name, scenario_type) {
  # Parse Climate from model name
  climate <- dplyr::case_when(
    grepl("persist", model_name) ~ "persistence",
    grepl("126",     model_name) ~ "ssp126",
    grepl("245",     model_name) ~ "ssp245",
    grepl("585",     model_name) ~ "ssp585"
  )
  # Parse fishing rate from model name
  frate <- ifelse(grepl("^F0", model_name), "F0", "Fmean")

  # SSB
  SSB_long <- as.data.frame(cbind(1991:2099, model$out_SSB[tdec, T2_adults] * 234769)) |>
    setNames(c("year", T2_species)) |>
    pivot_longer(cols = 2:4, names_to = "Species", values_to = "SSB")

  # Total biomass (juv + adu per species)
  totB_long <- as.data.frame(cbind(
    1991:2099,
    (model$out_Biomass[tdec, "walleye_pollock_juvenile"]      + model$out_Biomass[tdec, "walleye_pollock_adult"])      * 234769,
    (model$out_Biomass[tdec, "pacific_cod_juvenile"]          + model$out_Biomass[tdec, "pacific_cod_adult"])          * 234769,
    (model$out_Biomass[tdec, "arrowtooth_flounder_juvenile"]  + model$out_Biomass[tdec, "arrowtooth_flounder_adult"])  * 234769
  )) |>
    setNames(c("year", T2_species)) |>
    pivot_longer(cols = 2:4, names_to = "Species", values_to = "totB")

  # Annual catch
  annC_long <- as.data.frame(cbind(1991:2099, model$annual_Catch[, T2_adults] * 234769)) |>
    setNames(c("year", T2_species)) |>
    pivot_longer(cols = 2:4, names_to = "Species", values_to = "Catch")

  # Combine and tag
  result <- cbind(totB_long, SSB = SSB_long[["SSB"]], Catch = annC_long[["Catch"]])
  result$Climate  <- climate
  result$Frate    <- frate
  result$Region   <- "WGOA"
  result$Model    <- "Rpath"
  result$Scenario <- scenario_type

  return(result)
}

# ---------------------------------------------------------------------------- #
# 3. Loop over all model lists and combine ####
# ---------------------------------------------------------------------------- #
all_model_lists <- list(
  base        = list_ssps,
  pcod_rec    = list_ssps_pcod,
  pcod_rec_pp = list_ssps_pp
)

all_model_data <- do.call(rbind, lapply(names(all_model_lists), function(scenario_type) {
  model_list <- all_model_lists[[scenario_type]]
  do.call(rbind, lapply(names(model_list), function(model_name) {
    extract_model_data(model_list[[model_name]], model_name, scenario_type)
  }))
}))

# ---------------------------------------------------------------------------- #
# 4. Export ####
# ---------------------------------------------------------------------------- #
write.csv(all_model_data, "results/topic2_all_model_data.csv", row.names = FALSE)

# ---------------------------------------------------------------------------- #
# 5. Sanity check plot  ####
# ---------------------------------------------------------------------------- #

ggplot(all_model_data,
       aes(
         x = year,
         y = totB / 1e6,
         color = Climate,
         linetype = Frate
       )) +
  geom_line(alpha = 0.8) +
  facet_grid(Species ~ Scenario, scales = "free_y") +
  labs(
    x = "Year",
    y = "Total Biomass (million mt)",
    color = "Climate",
    linetype = "F rate"
  )

# plot SSB, Catch, B

all_model_data %>%
  pivot_longer(
    cols = c(totB, SSB, Catch),
    names_to = "Metric",
    values_to = "value"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("totB", "SSB", "Catch"))) %>%
  ggplot(aes(
    x = year,
    y = value / 1e6,
    color = Metric,
    group = interaction(Climate, Frate, Metric)
  )) +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  facet_grid(Species ~ Scenario, scales = "free_y") +
  labs(x = "Year", y = "Million mt", color = "Metric")
