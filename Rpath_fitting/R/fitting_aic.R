# ---------------------------------------------------------------------------- #
# AUTHORS: George (Andy) Whitehouse
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: gaw@uw.edu
# UPDATED: 3 November 2025
#
# Purpose: Function to calculate AIC and AICc for fitting with Rpath.
# ---------------------------------------------------------------------------- #

rfit_aic <- function(scenario, vartype, nll) {
  # the scenario object used in fitting
  scene <- scenario
  # number biomass observations
  nbio_obs <- dim(scene$fitting$Biomass)[1]
  # number catch observations
  ncatch_obs <- dim(scene$fitting$Catch)[1]
  # total number of observations
  Kobs <- nbio_obs + ncatch_obs
  print(c("no. observations", " ", Kobs))
  
  # number of parameters being fit
  nparams <- length(vartype)
  print(c("no. params", " ", nparams))
  # negative log-likelihood (e.g., as output from rsim.fit.run())
  nll <- nll
  
  # Conventional AIC assumes obs are independent
  # AIC
  fit_aic_obs <- 2*nll + 2*nparams
  # AICc
  fit_aicc_obs <- fit_aic_obs + (2*nparams*(nparams + 1))/(Kobs - nparams - 1)
  
  # AIC treating each time series as an observation due to lack of
  # independence (i.e., due to autocorrelation)
  # number of biomass time series
  nbiots <- length(rsim.fit.list.bio.series(scene)$all)
  # number of catch time series
  ncatchts <- length(rsim.fit.list.catch.series(scene)$groups)
  # total number of time series
  Kts <- nbiots + ncatchts
  print(c("no. time series", " ", Kts))
  # AIC
  fit_aic_ts <- 2*nll + 2*nparams
  # AICc
  fit_aicc_ts <- fit_aic_ts + (2*nparams*(nparams + 1))/(Kts - nparams - 1)
  
  AIC_out <- c(fit_aic_obs, fit_aicc_obs, fit_aic_ts, fit_aicc_ts)
  names(AIC_out) <- c("aic_obs", "aicc_obs", "aic_time_series", "aicc_time_series")
  return(AIC_out)
}
