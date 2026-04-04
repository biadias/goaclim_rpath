# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 02 April 2026
#
# code/05_run_sim_cod_rec.R
# Purpose: This script appliers the recruitment function to Pcod according to 
# Laurel and Rogers 2020 paper hatching success
#
# From https://doi.org/10.1139/cjfas-2019-0238 Laurel and Rogers 2020
# Our study suggests that poor hatch success due to a reduction in thermal habitat 
# suitability may have contributed to low recruitment during the recent marine heatwave in the GOA.
# ---------------------------------------------------------------------------- #

source("code/02_run_sim.R")

scene_bioen_cod <- scene_bioen_best


laurel_rogers_survival_proxy <- function(temp, method="cauchy"){
  if(method == "cauchy"){
    survival_multiplier <- 1 / (1 + ((temp - 4.192) / 2.2125)^2)
    # Note: The parameters for the Cauchy function (location = 4.192, scale = 2.2125)
  } else if(method == "gaussian"){
    survival_multiplier <- exp(-0.5 *((temp-4.5)/1.5)^2) 
    # Note: The parameters for the Gaussian function (mean = 4.5, sd = 1.5) 
  } else {
    stop("Invalid method. Choose 'cauchy' or 'gaussian'.")
  }
  return(survival_multiplier)
}



total_years <- length(all_years)
total_months <- total_years * 12

monthly_temps_ssp126 <- read.csv("data/climate/ssp126_wide_WGOA_temp_1000.csv")
monthly_temps_ssp245 <- read.csv("data/climate/ssp245_wide_WGOA_temp_1000.csv")
monthly_temps_ssp585 <- read.csv("data/climate/ssp585_wide_WGOA_temp_1000.csv")

actual_temps_vector <- monthly_temps_ssp126$btemp


force_pattern <- rep(1, total_months)


for (yr in 1:total_years) {
  for (m in 1:12) {
    month_idx <- (yr - 1) * 12 + m
    
    
    if(m %in% c(2,3,4)){
      current_temp <- actual_temps_vector[month_idx]
      force_pattern[month_idx] <- laurel_rogers_survival_proxy(current_temp, method="cauchy")
    }else{
        force_pattern[month_idx] <- 1
      }
  }
}

scene_bioen_cod$forcing$ForcedRecs[,"pacific_cod_adult"] <- force_pattern

runtest <- rsim.run(scene_bioen_cod, "AB", all_years)
run_test_nocod <- rsim.run(scene_bioen_best, "AB", all_years)

par(mfcol=c(3,2))
plot(run_test_nocod$out_Biomass[,"pacific_cod_adult"],         type="l", main="PCOD_adu_nc", ylab="Biomass")
plot(run_test_nocod$out_Biomass[,"arrowtooth_flounder_adult"], type="l", main="ATF_adu_nc", ylab="Biomass")
plot(run_test_nocod$out_Biomass[,"walleye_pollock_adult"],     type="l", main="WPOL_adu_nc", ylab="Biomass")
plot(runtest$out_Biomass[,"pacific_cod_adult"],                type="l", main="PCOD_adu", ylab="Biomass")
plot(runtest$out_Biomass[,"arrowtooth_flounder_adult"],        type="l", main="ATF_adu", ylab="Biomass")
plot(runtest$out_Biomass[,"walleye_pollock_adult"],            type="l", main="WPOL_adu", ylab="Biomass")
plot(runtest$out_Biomass[,"pacific_cod_juvenile"], type="l", main="Juvenile Pacific Cod", ylab="Biomass")

par(mfcol=c(1,2))
plot(runtest$out_Biomass[,"pacific_cod_juvenile"], type="l", main="Juvenile Pacific Cod", ylab="Biomass")
plot(runtest$out_Biomass[,"pacific_cod_adult"],                type="l", main="PCOD_adu", ylab="Biomass")


par(mfcol=c(1,2))
plot(run_test_nocod$out_Biomass[,"pacific_cod_juvenile"], type="l", main="Juvenile Pacific Cod", ylab="Biomass")
plot(run_test_nocod$out_Biomass[,"pacific_cod_adult"],                type="l", main="PCOD_adu", ylab="Biomass")
