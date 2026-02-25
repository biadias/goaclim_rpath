# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, George (Andy) Whitehouse and Bridget Ferriss
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
                            use_parallel=TRUE,
                            penalty_weight=0,
                            nll_details=TRUE){
  require(parallel)
  require(optimParallel)
  
  message(paste("Starting optimization for scenario:", scene_name))
  
  if(use_parallel){
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
        "hind_years",
        "species_vec",
        "vartype_vec"
      ), envir = environment()
    )
    
    ptm <- proc.time()
    opt_result <- optimParallel::optimParallel(par = start_vec,
                                              fn = rsim.fit.run,
                                              method = "L-BFGS-B",
                                              lower = -4,
                                              upper = 4,
                                              hessian = FALSE,
                                              control = list(maxit=1000, trace = 2, factr = 1e7),
                                              species = species_vec,
                                              vartype = vartype_vec,
                                              scene = scene_obj,
                                              run_method = "AB",
                                              run_years = hind_years,
                                              penalty_weight=penalty_weight,
                                              nll_details=FALSE)
    
    stopCluster(cl)
  } else {
    ptm <- proc.time()
    opt_result <- optim(par = start_vec,
                        fn = rsim.fit.run,
                        method = "L-BFGS-B",
                        lower = -4,
                        upper = 4,
                        hessian = FALSE,
                        control = list(maxit=1000,trace = 2, factr = 1e7),
                        species = species_vec,
                        vartype = vartype_vec,
                        scene = scene_obj,
                        run_method = "AB",
                        run_years = hind_years,
                        penalty_weight=penalty_weight,
                        nll_details=FALSE)
  }
  time_taken <- proc.time() - ptm
  message(paste("Optimization for scenario", scene_name, "completed in", round(time_taken[3]/60,2), "minutes."))
  
  new_scene <- scene_obj
  new_scene$params <- rsim.fit.apply(opt_result$par,
                                    species_vec,
                                    vartype_vec,
                                    scene_obj$params)
  final_run <- rsim.fit.run(opt_result$par,
                            species_vec,
                            vartype_vec,
                            scene = scene_obj,
                            run_method = "AB",
                            run_years = hind_years,
                            verbose = TRUE,
                            nll_details=TRUE)
  
  nll_value <- rsim.fit.run(NA,
                            NA,
                            NA,
                            scene = new_scene,
                            run_method = "AB",
                            run_years = hind_years,
                            verbose = FALSE,
                            nll_details = FALSE)
  
  aic_val <- rfit_aic(new_scene, vartype_vec, nll_value)
  likelihood_table <- rsim.fit.table(scene_obj, final_run)
  
  return(list(
    name=scene_name,
    opt_object = opt_result,
    new_scene = new_scene,
    final_run = final_run,
    nll = nll_value,
    nll_penalized = opt_result$value,
    aicc = aic_val,
    like_table = likelihood_table,
    time=time_taken
  ))
}