# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 26 January 2026
#
# R/01_fitting_optimization.R
# Purpose: Script for fitting of GOA Rpath for GOACLIM
# Rpath_fitting/R/GOA_fitting_setup_a_la_carte.r to set up base scenario objects.
# ---------------------------------------------------------------------------- #
# General instructions:
# First, fit each group separately for vulnerabilities (get realistic results?).
# Record the change in NLL and qualitatively describe change in fit of the group
# being fit.
# Second, fit the model by only fitting the most sensitive groups, plus fitting
# the preyvuls for detritus pools, and VVs for zooplankton and microbes. And do
# not fit more params than time series.
# ---------------------------------------------------------------------------- #

run_rpath_optim <- function(scene_obj,
                            scene_name,
                            species_vec,
                            vartype_vec,
                            start_vec,
                            hind_years,
                            use_parallel = TRUE,
                            penalty_weight_vul = 0,
                            penalty_weight_mzero = 0,
                            nll_details = TRUE) {
  require(parallel)
  require(optimParallel)
  
  message(paste("Starting optimization for scenario:", scene_name))
  
  #Sequential fitting, first vulnerability parameters, then M0s fits. Change 2/3/2026
  # Vulnerabilities
  idx_vul <- which(vartype_vec %in% c("predvul", "preyvul"))
  #M0s
  idx_mzero <- which(vartype_vec %in% c("mzero"))
  
  
  if (length(idx_vul) == 0 && length(idx_mzero) == 0) {
    stop("No parameters to fit. Check species_vec and vartype_vec inputs.")
  }
  
  new_scene <- scene_obj
  final_scene <- scene_obj
  opt_result_s1 <- NULL
  opt_result_s2 <- NULL
  final_opt_obj <-  NULL
  
  
  if (use_parallel) {
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    setDefaultCluster(cl)
    parallel::clusterSetRNGStream(cl, 12345)
    
    clusterEvalQ(cl, {
      library(Rpath)
      source("R/GOA_fitting_setup_a_la_carte.R")
      source("R/fitting_aic.R")
    })
    
    clusterExport(
      cl,
      varlist = c(
        "rsim.fit.run",
        "rsim.fit.apply",
        "rsim.fit.obj",
        "rsim.run",
        "scene_obj",
        "rfit_aic",
        "hind_years"
      ),
      envir = environment()
    )
  }
    ptm <- proc.time() 
    
    if (length(idx_vul) > 0) {
    message(paste("Fitting", length(idx_vul), "vulnerability parameters first."))
  s1_species <- species_vec[idx_vul]
  s1_vartype <- vartype_vec[idx_vul]
  s1_start   <- start_vec[idx_vul]
  
  
  if (use_parallel) {
    opt_result_s1 <- optimParallel::optimParallel(
      par = s1_start, #start_vec,
      fn = rsim.fit.run,
      method = "L-BFGS-B",
      lower = -4,
      upper = 4,
      hessian = FALSE,
      control = list(
        maxit = 1000,
        trace = 2,
        factr = 1e7),
      species = s1_species, #species_vec,
      vartype = s1_vartype, #vartype_vec,
      scene = scene_obj,
      run_method = "AB",
      run_years = hind_years,
      penalty_weight = penalty_weight_vul,
      nll_details = FALSE
    )
    
  } else {
    opt_result_s1 <- optim(
      par = s1_start, #start_vec,
      fn = rsim.fit.run,
      method = "L-BFGS-B",
      lower = -4,
      upper = 4,
      hessian = FALSE,
      control = list(
        maxit = 1000,
        trace = 2,
        factr = 1e7),
      species = s1_species, #species_vec,
      vartype = s1_vartype, #vartype_vec,
      scene = scene_obj,
      run_method = "AB",
      run_years = hind_years,
      penalty_weight = penalty_weight_vul,
      nll_details = FALSE
    )
  }
  message(paste("Optimization for vul"))#scene_name, "completed in", round(time_taken[3]/60,2), "minutes."))
  
stage1_scene <- scene_obj  
stage1_scene$params <- rsim.fit.apply(
    opt_result_s1$par,
    s1_species, #species_vec,
    s1_vartype, #vartype_vec,
    scene_obj$params)
  
  #CHANGE ####
  current_scene <- stage1_scene
  
  #final_opt_obj <- opt_result_s1
  #final_scene <- new_scene
  
} else{
  message("No vulnerability parameters to fit, skipping to M0 fitting.")
  current_scene <- scene_obj
}
  # M0s ####
    
if (length(idx_mzero) > 0) {
  message(paste(
    "Starting M0 optimization for scenario:",
    length(idx_mzero),
    "parameters."
  ))
  
  s2_species <- species_vec[idx_mzero]
  s2_vartype <- vartype_vec[idx_mzero]
  s2_start   <- start_vec[idx_mzero]
  
  if (use_parallel) {
    opt_result_s2 <- optimParallel::optimParallel(
      par = s2_start, #start_vec,
      fn = rsim.fit.run,
      method = "L-BFGS-B",
      lower = 0.001,
      upper = 1.5,
      hessian = FALSE,
      control = list(
        maxit = 1000,
        trace = 2,
        factr = 1e7
      ),
      species = s2_species, #species_vec,
      vartype = s2_vartype, #vartype_vec,
      scene = current_scene, #stage1 scene vv
      run_method = "AB",
      run_years = hind_years,
      penalty_weight = penalty_weight_mzero,
      nll_details = FALSE
    )
  } else {
    opt_result_s2 <- optim(
      par = s2_start, #start_vec,
      fn = rsim.fit.run,
      method = "L-BFGS-B",
      lower = 0.001,
      upper = 1.5,
      hessian = FALSE,
      control = list(
        maxit = 1000,
        trace = 2,
        factr = 1e7
      ),
      species = s2_species, #species_vec,
      vartype = s2_vartype, #vartype_vec,
      scene = current_scene, #stage1 scene vv
      run_method = "AB",
      run_years = hind_years,
      penalty_weight = penalty_weight_mzero,
      nll_details = FALSE
    )
  }
  
  message(paste("Optimization for M0 NLL"))#scene_name, "completed in", round(time_taken[3]/60,2), "minutes."))
  
#  final_scene <- new_scene
  final_scene <- current_scene
  final_scene$params <- rsim.fit.apply(
    opt_result_s2$par,
    s2_species, #species_vec,
    s2_vartype, #vartype_vec,
    current_scene$params
  )
  
  
#  final_opt_obj <- opt_result_s2
  
} else{
  message("No M0 parameters to fit, skipping M0 fitting.")
  final_scene <- current_scene
  #final_opt_obj <- opt_result_s1
  #if(is.null(final_opt_obj)) {
  #  final_opt_obj <- NULL
  #}
}
if (use_parallel) stopCluster(cl)

time_taken <- proc.time() - ptm
message(paste("Total optimization time", round(time_taken[3] / 60, 2), "minutes."))

#diagnostics and final run with best parameters



final_run <- rsim.fit.run(
  NA,#opt_result_s2$par,
  NA,#s2_species,
  NA,#s2_vartype,
  scene = final_scene,
  run_method = "AB",
  run_years = hind_years,
  verbose = TRUE,
  nll_details = TRUE
)

nll_value <- rsim.fit.run(
  NA,
  NA,
  NA,
  scene = final_scene,
  run_method = "AB",
  run_years = hind_years,
  verbose = FALSE,
  nll_details = FALSE,
  penalty_weight = 0
)

penalty_s1 <- 0
if (!is.null(opt_result_s1) && penalty_weight_vul>0){
  
  nll_s1_penalized <- rsim.fit.run(
    opt_result_s1$par,
    s1_species,
    s1_vartype,
    scene = new_scene,
    run_method = "AB",
    run_years = hind_years,
    verbose = FALSE,
    nll_details = FALSE,
    penalty_weight = penalty_weight_vul
  )
  
  nll_s1_unpenalized <- rsim.fit.run(
    opt_result_s1$par,
    s1_species,
    s1_vartype,
    scene = new_scene,
    run_method = "AB",
    run_years = hind_years,
    verbose = FALSE,
    nll_details = FALSE,
    penalty_weight = 0
  
  )
  penalty_s1 <- nll_s1_penalized - nll_s1_unpenalized
}
                  
penalty_s2 <- 0
if (!is.null(opt_result_s2) && penalty_weight_mzero>0){
  penalty_s2 <- opt_result_s2$value-nll_value
}

total_nll_penalized <- nll_value+penalty_s1+penalty_s2
  
  
#AICc
aic_val <- rfit_aic(final_scene, vartype_vec, nll_value)
likelihood_table <- rsim.fit.table(scene_obj, final_run)

return(
  list(
    name = scene_name,
    opt_object_s1 = opt_result_s1,
    opt_object_s2 = opt_result_s2,
    final_scene = final_scene, #new_scene=final_scene,
    final_run = final_run,
    nll = nll_value,
    nll_penalized = total_nll_penalized,
    penalties = c(vul= penalty_s1, mzero = penalty_s2),
    aicc = aic_val,
    like_table = likelihood_table,
    time = time_taken
  ))
}
