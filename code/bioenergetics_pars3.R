#------------------------------------------------------------------------------#
#AUTHORS: Bia Dias
#ORIGINAL AUTHORS: Andy Whitehouse code developed for ACLIM2 
#ORIGINAL CODE: bioen_pars2.R
#AFFILIATIONS: CICOES University of Washington/ Alaska Fisheries Science Center
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#
#Bioenergetic params for GOACLIM
#------------------------------------------------------------------------------#
#Bioenergetics

library(tidyverse)
library(here)
library(viridis)

#------------------------------------------------------------------------------#
# 1) Load Data ####
#------------------------------------------------------------------------------#
bioen_pars <- read.csv("Rpath_fitting/GOA/wgoa_bioenergetics_code/wgoa_bioen.csv", header=TRUE, sep=',', dec='.', row.names=1)
bioen_sp <- row.names(bioen_pars)
# habitat niche preference
fg_pref <- read.csv("data/bioenergetics/fg-env-preference-parameters_WGOA.csv", header=TRUE)

# Helper function to clean names for matching
clean_name <- function(x) {
  tolower(gsub(" ", "_", x))
}

# Flagdem loop. Flag Demersal is 1, Pelagic is 0. 
# We will use this to determine which temp to use for each species.

niche_lookup <- fg_pref %>%
  mutate(clean_ewe = clean_name(EwE_name)) %>%
  select(clean_ewe, flagdem) %>% 
  deframe() 

# Survey Temperatures (Bottom and SST)
hind_temp_raw <- read.csv("data/bioenergetics/species_weighted_temp_WGOA.csv", header=TRUE, 
                          sep=',', dec='.', row.names=1)

# Process Bottom Temps
hind_bt <- hind_temp_raw %>%
  mutate(race_clean = clean_name(race_group)) %>%
  select(year, race_clean, agg_weighted_bot_temp) %>%
  pivot_wider(names_from = race_clean, values_from = agg_weighted_bot_temp) %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  column_to_rownames(var="year")

# Process Surface Temps
hind_st <- hind_temp_raw %>%
  mutate(race_clean = clean_name(race_group)) %>%
  select(year, race_clean, agg_weighted_surf_temp) %>%
  pivot_wider(names_from = race_clean, values_from = agg_weighted_surf_temp) %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  column_to_rownames(var="year")

#------------------------------------------------------------------------------#
# 2) Determine Habitat Niche ####
#------------------------------------------------------------------------------#
mean_ref_temps <- setNames(numeric(length(bioen_sp)), bioen_sp)

for (sp in bioen_sp) {
  # 1. Determine if the current name (or its adult proxy) is in the data
  target_name <- if (sp %in% colnames(hind_bt)) sp else gsub("_juvenile", "_adult", sp)
  
  if (target_name %in% colnames(hind_bt)) {
    # 2. Get the habitat niche for the base species
    habitat_type <- niche_lookup[sp]
    
    # 3. Calculate directly from the source dataframes
    if (!is.na(habitat_type) && habitat_type == 0) {
      mean_ref_temps[sp] <- round(mean(hind_st[[target_name]]), 2)
    } else {
      mean_ref_temps[sp] <- round(mean(hind_bt[[target_name]]), 2)
    }
  } else {
    # 4. Global fallback if neither sp nor adult_proxy is found
    mean_ref_temps[sp] <- 5.0 
  }
}




# ---------------------------------------------------------------------------- #
# 3) Calculate Scaled Consumption (Kitchell) ####
# ---------------------------------------------------------------------------- #

xx <-
  (((
    log(bioen_pars$Q10) * (bioen_pars$Tmax - bioen_pars$Topt)^2
  )) / 400) *
  (1 + (1 + sqrt(40 / (
    log(bioen_pars$Q10) * (bioen_pars$Tmax - bioen_pars$Topt + 2)
  )))^2)
names(xx) <- bioen_sp
# Kitchell equation: rc is proportion of max consumption at a given temp
# temps to evaluate over
Ctemp <- seq(-2, 30, 0.01)
rc <- matrix(nrow = length(Ctemp), ncol = length(bioen_sp))
colnames(rc) <- bioen_sp
row.names(rc) <- Ctemp

#Kitchell equation: rc is the proportion of max consumption at a given temp. 

for (i in bioen_sp) {
  for (j in Ctemp) {
    rc[j, i] <- (((bioen_pars[i, 'Tmax'] - Ctemp[j]) / (bioen_pars[i, 'Tmax'] - bioen_pars[i, 'Topt']))^xx[i]) *
      exp(xx[i] * (1 - ((bioen_pars[i, 'Tmax'] - Ctemp[j]) / (bioen_pars[i, 'Tmax'] - bioen_pars[i, 'Topt'])
      )))
  }
}

# Alternatively, we can use species-specific biomass-weighted mean bottom temps
# from the survey.
# mean temperatures are the means of the sp.-specific survey time series
# Scale rc to the reference temperature (Surface for Pelagic, Bottom for Benthic)
rownames(rc) <- sprintf("%.2f", round(as.numeric(rownames(rc)), digits=2))

rc_scaled <- rc
for(i in bioen_sp){
  ref_t_str <- sprintf("%.2f", mean_ref_temps[i])
  rc_scaled[, i] <- rc[, i] / rc[ref_t_str, i]
}
rc_scaled[is.nan(rc_scaled)] <- 1e-08

# ---------------------------------------------------------------------------- #
# Kitchell curve plots ####
# ---------------------------------------------------------------------------- #
#par(mfrow = c(2, 1))
## rc unscaled
#plot(
#  Ctemp,
#  rc[, "arrowtooth_flounder_adult"],
#  type = 'n',
#  xlab = "Celcius",
#  ylab = "rc (proportion max consumption)",
#  ylim = c(0, 1.1)
#)
#for (i in 1:18) {
#  lines(Ctemp, rc[, i], col = viridis(18)[i], lwd = 2)
#}
## rc scaled to mean bottom temp 1991–1994
#plot(
#  Ctemp,
#  rc_scaled[, "arrowtooth_flounder_adult"],
#  type = 'n',
#  xlab = "Celcius",
#  ylab = "rc_scaled_b (proportion max consumption)",
#  ylim = c(0, max(rc_scaled, na.rm = TRUE)),
#  main = "rc_scaled to mean bottom temp_ATF"
#)
#for (i in 1:18) {
#  lines(Ctemp, rc_scaled[, i], col = viridis(18)[i], lwd = 2)
#  abline(v = mean_ref_temps[i], lty = 2, col = "gray50")
#}
#abline(h=1)
# # # rc scaled to mean surface temp 1991–1994
# # plot(
# #   Ctemp,
# #   rc_scaled_s[, "arrowtooth_adu"],
# #   type = 'n',
# #   xlab = "Celcius",
# #   ylab = "rc_scaled_s (proportion max consumption)",
# #   ylim = c(0, max(rc_scaled_s, na.rm = TRUE)),
# #   main = "rc_scaled to mean surface temp"
# # )
# # for (i in 1:26) {
# #   lines(Ctemp, rc_scaled_s[, i], col = viridis(26)[i], lwd = 2)
# #   abline(v = mean_sur_temps, lty = 2, col = "gray50")
# # }
# # abline(h=1)
#  
#  
#  # # rc scaled to bottom temp plots by species
#  par(mfrow = c(4, 8))
#  for (i in 1:31) {
#    plot(
#      Ctemp,
#      rc_scaled_b[, i],
#      type = 'l',
#      lwd = 2,
#      main = bioen_sp[i],
#      ylim = c(0, max(rc_scaled_b[, i], na.rm = TRUE)),
#      ylab = "rc_scaled (BT)",
#      xlab = "temp (C)"
#    )
#    abline(v = mean_bot_temps[i], lty = 2, col = "gray50")
#  }
# # rc scaled to surface temp plots by species
# par(mfrow = c(4, 8))
# for (i in 1:31) {
#   plot(
#     Ctemp,
#     rc_scaled_s[, i],
#     type = 'l',
#     lwd = 2,
#     main = bioen_sp[i],
#     ylim = c(0, max(rc_scaled_s[, i], na.rm = TRUE)),
#     ylab = "rc_scaled (ST)",
#     xlab = "temp (C)"
#   )
#   abline(v = mean_sur_temps[i], lty = 2, col = "gray50")
# }
# 
# 
# # By species both plots on the same graph
# par(mfrow = c(4, 8))
# for (i in 1:31) {
#   plot(
#     Ctemp,
#     rc_scaled_b[, i],
#     type = 'l',
#     lwd = 2,
#     main = bioen_sp[i],
#     ylim = c(0, max(rc_scaled_b[, i], na.rm = TRUE)),
#     ylab = "rc_scaled",
#     xlab = "temp (C)",
#     col = "blue"
#   )
#   abline(v = mean_bot_temps[i], lty = 2, col = "blue")
#   lines(Ctemp, rc_scaled_s[, i], lwd=2, col = "red", lty=2)
#   abline(v = mean_sur_temps[i], lty = 2, col = "red")
#   abline(h=1)
# }
# par(mfrow=c(1,1))

# Forced search consumption modifier 
# get species-specific consumption modifiers
# temp time series from GOACLIM hindcast
#**hindcast**: representing the final years of the spinup forced with observed oceanographic conditions to better represent historical conditions (1990 to 2020),
#**projection**: GFDL-ESM2M downscaled projection (2015 to 2099)
#**historical**: representing model spinup (1980 to 2014).




# ---------------------------------------------------------------------------- #
# 4) Calculate Scaled Respiration (Modified Arrhenius)  ####
# Modified Arrhenius equation from Blanchard et al. (2012)
# ---------------------------------------------------------------------------- #
Ktemp <- Ctemp + 273.15 # Kelvin
k <- 8.62 * 10^(-5)
E <- 0.63
c1 <- 25.55
tau <- exp(c1 - (E / (k * Ktemp)))

tau_scaled <- matrix(nrow = length(tau), ncol = length(bioen_sp), 
                     dimnames = list(sprintf("%.2f", Ctemp), bioen_sp))
for (i in bioen_sp) {
  ref_k <- mean_ref_temps[i] + 273.15
  ref_tau <- exp(c1 - (E / (k * ref_k)))
  tau_scaled[, i] <- tau / ref_tau
}
# ---------------------------------------------------------------------------- #
# 5) Calculate ForcedActResp (Modifier) ####
# ---------------------------------------------------------------------------- #
rc_rows <- row.names(rc_scaled)
# TotCons scaled by niche temp
TotCons <- matrix(nrow = length(Ctemp), ncol = length(bioen_sp), 
                  dimnames = list(rc_rows, bioen_sp))
for (i in bioen_sp) {
  TotCons[, i] <- rc_scaled[, i] * w.bal$QB[i] * w.bal$Biomass[i]
}
# Total respiration at sp.-specific reference means
TotResp_mean <- vector(mode = "numeric", length = length(bioen_sp))
names(TotResp_mean) <- bioen_sp
for (i in bioen_sp) {
  TotResp_mean[i] <- TotCons[sprintf("%.2f", mean_ref_temps[i]), i] * scene0$params$ActiveRespFrac[i]
}

# Total respiration curve
TotResp <- matrix(nrow = length(Ctemp), ncol = length(bioen_sp), dimnames = list(rc_rows, bioen_sp))
for (i in bioen_sp) {
  TotResp[, i] <- tau_scaled[, i] * TotResp_mean[i]
}

# ActiveRespFrac by temperature
ActiveRespFrac <- TotResp / TotCons

# Annual respiration modifier (ForcedActResp)
ForcedActResp <- matrix(nrow = length(Ctemp), ncol = length(bioen_sp), dimnames = list(rc_rows, bioen_sp))
for(i in bioen_sp){
  ForcedActResp[,i] <- ActiveRespFrac[,i] / scene0$params$ActiveRespFrac[i]
}
ForcedActResp[is.nan(ForcedActResp)] <- 1e04

# Plot
# par(mfrow = c(4, 8))
# for (i in 1:31) {
#   plot(
#     Ctemp,
#     ForcedActResp[, i],
#     type = 'l',
#     lwd = 2,
#     main = bioen_sp[i],
#     ylim = c(0, 2),
#     ylab = "ForcedActResp_b",
#     xlab = "temp (C)",
#     col = "blue"
#   )
#   abline(v = bioen_pars[i,5], lty = 2, col = "blue")
#   lines(Ctemp, ForcedActResp[, i], lwd=2, col = "red", lty=2)
#   abline(v = bioen_pars[i,4], lty = 2, col = "red")
# }
# par(mfrow=c(1,1))
 # ---------------------------------------------------------------------------- #
 # 6) Hindcast Forcing Matrices (tdc and tdr) ####
 # ---------------------------------------------------------------------------- #
 # Process ROMS temperatures
 roms_hind_temp <- read.csv("Rpath_fitting/GOA/wgoa_data_rpath_fitting/Long_WGOA_temp_monthly_1000.csv", 
                            header=TRUE, sep=',', dec='.') %>% 
   filter(simulation == "ssp126", year > 1990 & year < 2021) %>% 
   select(year, month, depthclass, area_weighted_temp) %>% 
   pivot_wider(names_from = depthclass, values_from = area_weighted_temp) %>% 
   rename(temp_b5 = Bottom, temp_s5 = Surface)

 # make sure the scale is in two decimal places (e.g. 2.20 instead of 2.2) to match the row names of rc_scaled and ForcedActResp
 row.names(ForcedActResp) <- sprintf("%.2f", as.numeric(row.names(ForcedActResp)))
 
 # Consumption Modifier (tdc) and Respiration Modifier (tdr)
 tdc_hind <- matrix(nrow = nrow(roms_hind_temp), ncol = length(bioen_sp), dimnames = list(1:nrow(roms_hind_temp), bioen_sp))
 tdr_hind <- tdc_hind
 
 for (i in 1:nrow(roms_hind_temp)) {
   for (j in bioen_sp) {
     # Habitat selection
     habitat_type <- niche_lookup[j]
     
     # Pick Surface (s5) for Pelagic (0), else Bottom (b5)
     raw_temp <- ifelse(!is.na(habitat_type) && habitat_type == 0, 
                        roms_hind_temp$temp_s5[i], 
                        roms_hind_temp$temp_b5[i])
     
     # Handle NA temperatures (Critical check)
     if (is.na(raw_temp)) {
       # If temp is missing, you can use a default mean or the previous step
       # Using 5.0 as a fallback example
       raw_temp <- 5.0 
     }
     
     # Clamp and Format
     target_temp <- max(-2, min(30, raw_temp))
     temp_str <- sprintf("%.2f", round(target_temp, digits = 2))
     
     # 5. Safe Lookup
     if (temp_str %in% rownames(rc_scaled) && j %in% colnames(rc_scaled)) {
       tdc_hind[i, j] <- rc_scaled[temp_str, j]
       tdr_hind[i, j] <- ForcedActResp[temp_str, j]
     } else {
       # This will tell you exactly which species or temp string is failing
       stop(paste0("Lookup failed at row ", i, ": Species '", j, "' or Temp '", temp_str, "' not found."))
     }
   }
 }
 
 
 