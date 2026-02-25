# ---------------------------------------------------------------------------- #
# AUTHORS: George (Andy) Whitehouse
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: gaw@uw.edu
# UPDATED: 1 August 2025
#
# Purpose: Function for first pass fitting of vulnerabilities . Run 
# Rpath_fitting/R/EBS_fitting.r to set up base scenario object.
# ---------------------------------------------------------------------------- #

vv_fit_loop <- function(bal, scene, fit_years, base_nll) {
  
  bal       <- bal       # balanced Ecopath
  scene_vv  <- scene     # scenario object
  fit_years <- fit_years # years over which we are fitting
  
  # output objects
  # vulnerability table
  vv_out     <- matrix(NA, ncol = 8, nrow = (bal$NUM_LIVING + bal$NUM_DEAD))
  colnames(vv_out)  <- c("spnum","preyvul","predvul","nll","dnll",
                         "nll_groupB","dnll_groupB","convergence")
  row.names(vv_out) <- scene_vv$params$spname[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1)]
  vv_out[,1] <- scene_vv$params$spnum[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1)]
  vv_out     <- as.data.frame(vv_out)
  # table of group changes in NLL
  dnll_Biomass   <- matrix(NA, ncol = (bal$NUM_LIVING + bal$NUM_DEAD), 
                       nrow = (bal$NUM_LIVING + bal$NUM_DEAD))
  row.names(dnll_Biomass) <- scene_vv$params$spname[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1)]
  colnames(dnll_Biomass)  <- paste0("dnll_", 
                                scene$params$spname[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1)])
  dnll_Biomass   <- as.data.frame(dnll_Biomass)
  
  for(i in 1:dim(vv_out)[1]) {
    # specify functional group
    print(paste0(row.names(vv_out)[i], " ", i))
    func_grp     <- row.names(vv_out)[i]

    opt_vartype  <- c(rep("preyvul", length(func_grp)), 
                      rep("predvul",length(func_grp))) 
    opt_species  <- c(func_grp, func_grp)
    start_vv     <- rep(0,length(opt_species))
    # optimize
    opt_call <- optim(start_vv, rsim.fit.run, method = "L-BFGS-B",
                      lower = -4, upper = 4, hessian = FALSE,
                      control = list(trace = 2, factr = 1e7), #changed from factr=1e11
                      species = opt_species, vartype = opt_vartype,
                      scene = scene_vv, run_method = "AB", run_years = fit_years)
    # store vulnerabilities
    vv_out[i, 2:3] <- opt_call$par
    print(paste0(row.names(vv_out)[i], " ", opt_call$par))
    # store convergence
    vv_out[i, "convergence"]   <- as.character(opt_call$convergence == 0)
    print(opt_call$convergence == 0)
    # apply params to a new scene
    scene_apply        <- scene_vv
    scene_apply$params <- rsim.fit.apply(opt_call$par, opt_species, 
                                         opt_vartype, scene_vv$params)
    # get new likelihood
    opt_nll <- rsim.fit.run(NA, NA, NA, scene = scene_apply, run_method = "AB", 
                            run_years = hind_years, verbose = F)
    vv_out[i, 4] <- opt_nll
    # store change in likelihood
    vv_out[i, 5] <- opt_nll - base_nll
    # get the individual group likelihood
    run.vv_out <- rsim.fit.run(opt_call$par, opt_species, opt_vartype, 
                               scene = scene_vv, run_method = "AB", 
                               run_years = hind_years, verbose = T) 
    group_nll_vv_out <- rsim.fit.table(scene_vv, run.vv_out)
    vv_out[i, 6] <- group_nll_vv_out[i+1, "Biomass"]
    # get the change in individual group likelihood
    vv_out[i, 7] <- group_nll_vv_out[i+1, "Biomass"] - base.like[i+1, "Biomass"]
    # apply changes to individual group NLL's to dnll_Biomass table
    dnll_Biomass[,i] <- group_nll_vv_out[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1), "Biomass"] - 
                                base.like[2:(bal$NUM_LIVING + bal$NUM_DEAD + 1), "Biomass"]
  }
  vv_out <- list(vv_out, dnll_Biomass)
  return(vv_out)

} 
