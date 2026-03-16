# ---------------------------------------------------------------------------- #
# AUTHORS: George (Andy) Whitehouse
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: gaw@uw.edu
# UPDATED: 15 September 2025
#
# Purpose: Function for first pass fitting of Mzero. 
# Arguments:
# 1) bal: balanced Rpath model
# 2) scene: Rpath scenario object
# 3) fit_years: The years over which the fitting is taking place
# 4) base_nll: The negative log-likelihood of the "base" model or the NLL of the
#              model that you are making comparisons with.
# 5) base_like: A table of functional group NLL's for the "base" model (as
#               produced by rsim.fit.table())
# 6) fit_groups: the group or vector of groups that are going to be fit
# ---------------------------------------------------------------------------- #

mzero_fit_loop <- function(bal, scene, fit_years, base_nll, base_like, fit_groups) {
  
  bal          <- bal        # balanced Ecopath
  scene_mzero  <- scene      # scenario object
  fit_years    <- fit_years  # years over which we are fitting
  fit_groups   <- fit_groups # functional groups to have mzero fit
  
  # output objects
  # mzero table
  mzero_out            <- matrix(NA, ncol = 9, nrow = (length(fit_groups)))
  colnames(mzero_out)  <- c("spnum", "mzero_fit", "mzero_bal", "nll", "dnll",
                            "nllB_group", "nllB_group_base", "dnllB_group",
                            "convergence")
  row.names(mzero_out) <- scene_mzero$params$spname[fit_groups]
  mzero_out[,1]        <- scene_mzero$params$spnum[fit_groups]
  mzero_out            <- as.data.frame(mzero_out)
  # table of functional group changes in NLL
  dnll_Biomass   <- matrix(NA, 
                           ncol = length(fit_groups), 
                           nrow = (bal$NUM_LIVING + bal$NUM_DEAD))
  row.names(dnll_Biomass) <- scene_mzero$params$spname[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1)]
  colnames(dnll_Biomass)  <- paste0("dnll_", 
                                    fit_groups)
  dnll_Biomass            <- as.data.frame(dnll_Biomass)
  
  for(i in 1:dim(mzero_out)[1]) {
    # specify functional group
    print(paste0(row.names(mzero_out)[i], " ", i))
    func_grp     <- row.names(mzero_out)[i]
    
    opt_vartype  <- "mzero" 
    opt_species  <- func_grp
    print(func_grp)
    start_mzero     <- 0.3
    # optimize
    opt_call <- optim(start_mzero, rsim.fit.run, method = "L-BFGS-B",
                      lower   = -4.5, upper = 4.5, hessian = FALSE,
                      control = list(trace = 2, factr = 1e13),
                      species = opt_species, 
                      vartype = opt_vartype,
                      scene   = scene_mzero, 
                      run_method = "AB", run_years = fit_years)
    # store mzero
    mzero_out[i, 2] <- opt_call$par
    print(paste0(row.names(mzero_out)[i], " ", opt_call$par))
    # what was the Ecopath mzero
    mzero_out[i, 3] <- (1 - bal$EE[func_grp])
    # store convergence
    mzero_out[i, "convergence"]   <- as.character(opt_call$convergence == 0)
    print(opt_call$convergence == 0)
    # apply params to a new scene
    scene_apply        <- scene_mzero
    scene_apply$params <- rsim.fit.apply(opt_call$par, 
                                         opt_species, 
                                         opt_vartype, 
                                         scene_mzero$params)
    # get new likelihood
    opt_nll <- rsim.fit.run(NA, NA, NA, 
                            scene = scene_apply, 
                            run_method = "AB", 
                            run_years = hind_years, verbose = F)
    mzero_out[i, 4] <- opt_nll
    # store change in likelihood
    mzero_out[i, 5] <- opt_nll - base_nll
    # get the individual group likelihood
    run.mzero_out <- rsim.fit.run(opt_call$par, 
                                  opt_species, 
                                  opt_vartype, 
                                  scene = scene_mzero, 
                                  run_method = "AB", 
                                  run_years = hind_years, verbose = T) 
    group_nll_mzero_out <- rsim.fit.table(scene_mzero, run.mzero_out)
    mzero_out[func_grp, 6] <- group_nll_mzero_out[func_grp, "Biomass"]
    mzero_out[func_grp, 7] <- base_like[func_grp, "Biomass"]
    # get the change in individual group likelihood
    mzero_out[func_grp, 8] <- group_nll_mzero_out[func_grp, "Biomass"] - 
                                  base.like[func_grp, "Biomass"]
    # apply changes to individual group NLL's to dnll_Biomass table
    dnll_Biomass[,i] <- group_nll_mzero_out[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1), "Biomass"] - 
                                  base.like[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1), "Biomass"]
  }
  mzero_out <- list(mzero_out, dnll_Biomass)
  return(mzero_out)
  
} 
