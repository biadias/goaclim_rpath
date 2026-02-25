#------------------------------------------------------------------------------#
#AUTHORS: Bia Dias
#ORIGINAL AUTHORS: Andy Whitehouse code developed for ACLIM2 
#ORIGINAL CODE: bioen_pars2.R
#AFFILIATIONS: CICOES University of Washington/ Alaska Fisheries Science Center
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#
#Bioenergetic params for GOACLIM
#------------------------------------------------------------------------------#

library(tidyverse)
library(here)
library(viridis)
library(viridisLite)



bioen_pars <- read.csv("GOA/wgoa_bioenergetics_code/WGOA_bioen.csv", header=TRUE, sep=',', 
                       dec='.', row.names=1) #this file is valid for both EGOA and WGOA
# X parameter in Kitchell equation
bioen_sp <- row.names(bioen_pars)
# bioen_sp <- bioen_sp[ !bioen_sp == "octopus" & !bioen_sp == "squids"]
# shortrakers have not been caught in the shelf BTS. Applying rougheye mean temps
# as a proxy
#bioen_pars["shortraker_rock",4:5] <- bioen_pars["rougheye_rock",4:5]
# species-specific temp time series from BT survey
hind_bt <- read.csv(
  "GOA/wgoa_bioenergetics_code/species_weighted_temp_WGOA.csv",
  header = TRUE,
  sep = ',',
  dec = '.',
  row.names = 1
) %>%
  select(c(year, race_group, agg_weighted_bot_temp)) %>%
  pivot_wider(names_from = race_group, values_from = agg_weighted_bot_temp)

#hind_st <- read.csv(
#  "WGOA_source_data/species_weighted_temp_WGOA.csv",
#  header = TRUE,
#  sep = ',',
#  dec = '.',
#  row.names = 1
#) %>%
#  select(c(year, race_group, agg_weighted_surf_temp)) %>%
#  pivot_wider(names_from = race_group, values_from = agg_weighted_surf_temp)

# replace missing (NA) with mean
hind_bt <- hind_bt %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>% 
  column_to_rownames(var="year")

#for(i in 1:(dim(hind_bt)[2])){
#  hind_bt[,i][is.na(hind_bt[,i])] <- mean(hind_bt[,i], na.rm=TRUE)
#}
#for(i in 1:(dim(hind_st)[2])){
#  hind_st[,i][is.na(hind_st[,i])] <- mean(hind_st[,i], na.rm=TRUE)
#}
mean_bot_temps <- round(colMeans(hind_bt),digits=2)
#mean_sur_temps <- round(colMeans(hind_st),digits=2)
#BD do the same for the other functional groups that have proxy species ####
# adding shortraker_rock using rougheye temps as proxy
#mean_bot_temps["shortraker_rock"] <- mean_bot_temps["rougheye_rock"]
#mean_sur_temps["shortraker_rock"] <- mean_sur_temps["rougheye_rock"]

# indicator species for multi-species functional groups
# oth_flatfish is represented by rex_sole
# oth_rockfish is represented by dusky_rockfish

# xx is a term in the Kitchell equation
xx <-
  (((
    log(bioen_pars$Q10) * (bioen_pars$Tmax - bioen_pars$Topt) ^ 2
  )) / 400) *
  (1 + (1 + sqrt(40 / (
    log(bioen_pars$Q10) * (bioen_pars$Tmax - bioen_pars$Topt + 2)
  ))) ^ 2)
names(xx) <- bioen_sp
# Kitchell equation: rc is proportion of max consumption at a given temp
# temps to evaluate over
Ctemp <- seq(-2, 30, 0.01)
rc <- matrix(nrow = length(Ctemp), ncol = length(bioen_sp))
colnames(rc) <- bioen_sp
row.names(rc) <- Ctemp
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
rc_scaled_b <- matrix(nrow = dim(rc)[1], ncol = dim(rc)[2])
colnames(rc_scaled_b) <- bioen_sp
names(mean_bot_temps) <- bioen_sp
for(i in 1:(length(bioen_sp))){
  rc_scaled_b[,bioen_sp[i]] <- rc[,bioen_sp[i]]/rc[as.character(mean_bot_temps[bioen_sp[i]]),bioen_sp[i]]
}
# have to name the rows after creating the matrix, and in this fashion or it doesn't
# work. It eliminates trailing zeroes which throws everything off in later loops.
rc_rows <- sprintf("%.2f", round(Ctemp, digits=2))
row.names(rc_scaled_b) <- rc_rows
# mean biomass-weighted species-specific surface temp 1991-1994
#rc_scaled_s <- matrix(nrow = dim(rc)[1], ncol = dim(rc)[2])
#colnames(rc_scaled_s) <- bioen_sp
#names(mean_sur_temps) <- bioen_sp
#for(i in 1:(length(bioen_sp))){
#  rc_scaled_s[,bioen_sp[i]] <- rc[,bioen_sp[i]]/rc[as.character(mean_sur_temps[bioen_sp[i]]),bioen_sp[i]]
#}
#row.names(rc_scaled_s) <- rc_rows
# replace NaN with a really small value
rc_scaled_b[is.nan(rc_scaled_b)] <- 1e-08
#rc_scaled_s[is.nan(rc_scaled_s)] <- 1e-08

# Kitchell curve plots -------------------------------------------------------
#  par(mfrow = c(2, 1))
#  # rc unscaled
#  plot(
#    Ctemp,
#    rc[, "arrowtooth_flounder_adult"],
#    type = 'n',
#    xlab = "Celcius",
#    ylab = "rc (proportion max consumption)",
#    ylim = c(0, 1.1)
#  )
#  for (i in 1:18) {
#    lines(Ctemp, rc[, i], col = viridis(18)[i], lwd = 2)
#  }
#  # rc scaled to mean bottom temp 1991–1994
#  plot(
#    Ctemp,
#    rc_scaled_b[, "arrowtooth_flounder_adult"],
#    type = 'n',
#    xlab = "Celcius",
#    ylab = "rc_scaled_b (proportion max consumption)",
#    ylim = c(0, max(rc_scaled_b, na.rm = TRUE)),
#    main = "rc_scaled to mean bottom temp"
#  )
#  for (i in 1:18) {
#    lines(Ctemp, rc_scaled_b[, i], col = viridis(18)[i], lwd = 2)
#    abline(v = mean_bot_temps[i], lty = 2, col = "gray50")
#  }
#  abline(h=1)
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

# Forced search consumption modifier ####
# get species-specific consumption modifiers
# temp time series from GOACLIM hindcast
#**hindcast**: representing the final years of the spinup forced with observed oceanographic conditions to better represent historical conditions (1990 to 2020),
#**projection**: GFDL-ESM2M downscaled projection (2015 to 2099)
#**historical**: representing model spinup (1980 to 2014).

 
roms_hind_temp <- read.csv("GOA/wgoa_data_rpath_fitting/Long_WGOA_temp_monthly_1000.csv", 
                           header=TRUE, sep=',', dec='.') %>% 
  filter(simulation=="ssp126", year>=1990 & year<2021) %>% 
  select(c( year, month, depthclass,area_weighted_temp)) %>% 
  pivot_wider(names_from = depthclass, values_from = area_weighted_temp) %>% 
  rename(temp_b5=Bottom,  temp_s5=Surface)
roms_hind_temp$tstep <- 1:372

roms_hind_npz <- read.csv("GOA/wgoa_bioenergetics_code/Long_WGOA_B_summary_month1000_v2_corrected.csv", 
                           header=TRUE, sep=',', dec='.') %>% 
  filter(simulation=="ssp126", year>=1990 & year<2021) %>% 
  select(c( year, month, varname,biomass_tonnes_km2)) %>% 
  pivot_wider(names_from = varname, values_from = biomass_tonnes_km2) %>% 
  rename(cop=Cop, eup=Eup, mzl=MZL, mzs=MZS, nca=NCa, phl=PhL, phs=PhS)  
  #select(cop, eup, mzl, phl, phs) #to match aclim_rpath
roms_hind_npz$tstep <- 1:372

goaclim_hind_raw <- roms_hind_temp %>% left_join(roms_hind_npz, by= join_by(tstep,year, month)) %>% 
  select(tstep, year, month, cop, eup, mzl, phl, phs, temp_b5, temp_s5)


goaclim_hind_bt <- as.matrix(goaclim_hind_raw$temp_b5)
row.names(goaclim_hind_bt) <- goaclim_hind_raw$tstep
row.names(goaclim_hind_bt) <- goaclim_hind_raw$tstep
goaclim_hind_st <- as.matrix(goaclim_hind_raw$temp_s5)
row.names(goaclim_hind_st) <- goaclim_hind_raw$tstep
# for some reason hind_st[52] gets rounded off to "-0.00" which is a problem in
# the consumption modifier below
goaclim_hind_st[52] <- 0.00
# bottom temperature
tdc_hind_bt <- matrix(nrow=dim(goaclim_hind_bt)[1], ncol=length(bioen_sp))
colnames(tdc_hind_bt) <- bioen_sp 
for(i in 1:dim(goaclim_hind_bt)[1]){
  for(j in bioen_sp){
    tdc_hind_bt[i,j] <- rc_scaled_b[sprintf("%.2f", round(goaclim_hind_bt[i], digits=2)),j]
  }
}
row.names(tdc_hind_bt) <- row.names(goaclim_hind_bt)
## surface temperature
#tdc_hind_st <- matrix(nrow=dim(goaclim_hind_st)[1], ncol=length(bioen_sp))
#colnames(tdc_hind_st) <- bioen_sp # w/o shortraker rockfish—never caught on shelf BT survey
#for(i in 1:dim(goaclim_hind_st)[1]){
#  for(j in bioen_sp){
#    tdc_hind_st[i,j] <- rc_scaled_s[sprintf("%.2f", round(goaclim_hind_st[i], digits=2)),j]
#  }
#}
#row.names(tdc_hind_st) <- row.names(goaclim_hind_st)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Respiration
# Modified Arrhenius equation from Blanchard et al. (2012)

# temps to evaluate (Kelvin = C + 273.15)
Ktemp <- seq(271, 303, 0.01) # from -2 to 30 Celcius
k  <- 8.62 * 10 ^ (-5) # Boltzmann constant
E  <- 0.63             # activation energy
c1 <- 25.55            # constant
tau <- matrix(nrow = length(Ktemp), ncol = 1)
row.names(tau) <- Ktemp
for (i in 1:(length(Ktemp))) {
  tau[i] <- exp(c1 - (E / (k * Ktemp[i])))
}

# tau scaled by sp.-specific mean bottom temps over survey time series
tau_scaled_b <- matrix(nrow=length(tau), ncol=dim(bioen_pars)[1])
for (i in 1:(dim(bioen_pars)[1])) {
  tau_scaled_b[,i] <- tau/(tau[as.character(mean_bot_temps[i] + 273.15),])
}
colnames(tau_scaled_b) <- bioen_sp
row.names(tau_scaled_b) <- Ktemp
# tau scaled
# Using species-specific biomass-weighted mean surface temp from survey
#tau_scaled_s <- matrix(nrow=length(tau), ncol=dim(bioen_pars)[1])
#for (i in 1:(dim(bioen_pars)[1])) {
#  tau_scaled_s[,i] <- tau/(tau[as.character(mean_sur_temps[i] + 273.15),])
#}
#colnames(tau_scaled_s) <- bioen_sp
#row.names(tau_scaled_s) <- Ktemp

#------------------------------------------------------------------------------#
# Total Consumption scaled by bottom temp (TotCons = rc_scaled * bal$QB * bal$Biomass)
TotCons_b <- matrix(nrow = dim(rc_scaled_b)[1], ncol = dim(rc_scaled_b)[2])
colnames(TotCons_b) <- bioen_sp
for (i in bioen_sp) {
  TotCons_b[, i] <- rc_scaled_b[, i] * w.bal$QB[i] * w.bal$Biomass[i]
}
# row.names(TotCons_b) <- Ctemp
row.names(TotCons_b) <- rc_rows
# Total Consumption scaled by surface temp (TotCons = rc_scaled * bal$QB * bal$Biomass)
#TotCons_s <- matrix(nrow = dim(rc_scaled_s)[1], ncol = dim(rc_scaled_s)[2])
#colnames(TotCons_s) <- bioen_sp
#for (i in bioen_sp) {
#  TotCons_s[, i] <- rc_scaled_s[, i] * w.bal$QB[i] * w.bal$Biomass[i]
#}
## row.names(TotCons_s) <- Ctemp
#row.names(TotCons_s) <- rc_rows


#------------------------------------------------------------------------------#
# Total respiration at sp.-specific mean temps from survey time series
# bottom temp
TotResp_btmean <- vector(mode = "numeric", length = length(bioen_sp))
names(TotResp_btmean) <- bioen_sp
for (i in bioen_sp) {
  TotResp_btmean[i] <- TotCons_b[sprintf("%.2f", mean_bot_temps[i]), i] * scene0$params$ActiveRespFrac[i] # proportion of consumption lost to respirartion
} #IMPORTANT chose scenario #####


## surface temp
#TotResp_stmean <- vector(mode = "numeric", length = length(bioen_sp))
#names(TotResp_stmean) <- bioen_sp
#for (i in bioen_sp) {
#  TotResp_stmean[i] <- TotCons_s[sprintf("%.2f", mean_sur_temps[i]), i] * scene$params$ActiveRespFrac[i]
#}

#------------------------------------------------------------------------------#
# Total respiration by temperature curve
# bottom temperature
TotResp_b <- matrix(nrow = dim(tau_scaled_b)[1], ncol = length(bioen_sp))
colnames(TotResp_b) <- bioen_sp
row.names(TotResp_b) <- row.names(tau_scaled_b)
for (i in bioen_sp) {
  TotResp_b[, i] <- tau_scaled_b[,i] * (TotResp_btmean[i] / tau_scaled_b[as.character(mean_bot_temps[i] + 273.15), i]) #transformation to Kelvin
}
## surface temperature
#TotResp_s <-  matrix(nrow = dim(tau_scaled_s)[1], ncol = length(bioen_sp))
#colnames(TotResp_s) <- bioen_sp
#row.names(TotResp_s) <- row.names(tau_scaled_s)
#for (i in bioen_sp) {
#  TotResp_s[, i] <- tau_scaled_s[,i] * (TotResp_stmean[i] / tau_scaled_s[as.character(mean_sur_temps[i] + 273.15), i])
#}

#------------------------------------------------------------------------------#
# ActiveRespFrac by temperature
# bottom temperature
ActiveRespFrac_b <- matrix(nrow=dim(TotResp_b)[1], ncol=dim(TotResp_b)[2])
colnames(ActiveRespFrac_b) <- bioen_sp
for(i in bioen_sp){
  ActiveRespFrac_b[,i] <- TotResp_b[,i]/TotCons_b[,i]
}
row.names(ActiveRespFrac_b) <- rc_rows
## surface temperature
#ActiveRespFrac_s <- matrix(nrow=dim(TotResp_s)[1], ncol=dim(TotResp_s)[2])
#colnames(ActiveRespFrac_s) <- bioen_sp
#for(i in bioen_sp){
#  ActiveRespFrac_s[,i] <- TotResp_s[,i]/TotCons_s[,i]
#}
#row.names(ActiveRespFrac_s) <- rc_rows

# plot ActiveRespFrac
# par(mfrow = c(4, 8))
# for (i in 1:31) {
#   plot(
#     Ctemp,
#     ActiveRespFrac_b[, i],
#     type = 'l',
#     lwd = 2,
#     main = bioen_sp[i],
#     ylim = c(0, 1),
#     ylab = "ActiveRespFrac",
#     xlab = "temp (C)",
#     col = "blue"
#   )
#   abline(v = bioen_pars[i,5], lty = 2, col = "blue")
#   lines(Ctemp, ActiveRespFrac_s[, i], lwd=2, col = "red", lty=2)
#   abline(v = bioen_pars[i,4], lty = 2, col = "red")
# }
# par(mfrow=c(1,1))

#------------------------------------------------------------------------------#
# Annual respiration modifier (ForcedActResp)
# bottom temperature
ForcedActResp_b <- matrix(nrow=dim(ActiveRespFrac_b)[1], ncol=dim(ActiveRespFrac_b)[2])
colnames(ForcedActResp_b) <- bioen_sp
for(i in bioen_sp){
  ForcedActResp_b[,i] <- ActiveRespFrac_b[,i]/scene0$params$ActiveRespFrac[i]
}
row.names(ForcedActResp_b) <- rc_rows
## surface temperature
#ForcedActResp_s <- matrix(nrow=dim(ActiveRespFrac_s)[1], ncol=dim(ActiveRespFrac_s)[2])
#colnames(ForcedActResp_s) <- bioen_sp
#for(i in bioen_sp){
#  ForcedActResp_s[,i] <- ActiveRespFrac_s[,i]/scene$params$ActiveRespFrac[i]
#}
#row.names(ForcedActResp_s) <- rc_rows

# replace NaN with a really large value
ForcedActResp_b[is.nan(ForcedActResp_b)] <- 1e04
#ForcedActResp_s[is.nan(ForcedActResp_s)] <- 1e04

# Plot
# par(mfrow = c(4, 8))
# for (i in 1:31) {
#   plot(
#     Ctemp,
#     ForcedActResp_b[, i],
#     type = 'l',
#     lwd = 2,
#     main = bioen_sp[i],
#     ylim = c(0, 2),
#     ylab = "ForcedActResp_b",
#     xlab = "temp (C)",
#     col = "blue"
#   )
#   abline(v = bioen_pars[i,5], lty = 2, col = "blue")
#   lines(Ctemp, ForcedActResp_s[, i], lwd=2, col = "red", lty=2)
#   abline(v = bioen_pars[i,4], lty = 2, col = "red")
# }
# par(mfrow=c(1,1))


#------------------------------------------------------------------------------#
# ForcedActResp modifier
# get species-specific respiration modifiers
# bottom temperature 
tdr_hind_bt <- matrix(nrow=dim(goaclim_hind_bt)[1], ncol=length(bioen_sp))
colnames(tdr_hind_bt) <- bioen_sp 
for(i in 1:dim(goaclim_hind_bt)[1]){
  for(j in bioen_sp){
    tdr_hind_bt[i,j] <- ForcedActResp_b[sprintf("%.2f", round(goaclim_hind_bt[i], digits=2)),j]
  }
}
row.names(tdr_hind_bt) <- row.names(goaclim_hind_bt)

#tdr_hind_bt is the final file ####

## surface temperature (hind_st from conusmption modifier above)
#tdr_hind_st <- matrix(nrow=dim(goaclim_hind_st)[1], ncol=length(bioen_sp))
#colnames(tdr_hind_st) <- bioen_sp 
#for(i in 1:dim(aclim_hind_st)[1]){
#  for(j in bioen_sp){
#    tdr_hind_st[i,j] <- ForcedActResp_s[sprintf("%.2f", round(aclim_hind_st[i], digits=2)),j]
#  }
#}
#row.names(tdr_hind_st) <- row.names(aclim_hind_st)


