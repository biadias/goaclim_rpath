# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias, Andy Whitehouse
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 08 July 2026
#
# Function to run climate-enhanced Rpath simulations for ACLIM 3.0.
# This function setup to accept scenario objects generated with
# F_clim_sim_scene_prim_prod.R
#
# Arguments:
# 1) ssps: climate scenario label, e.g. "persist", "126", "245", "585"
# 2) scene: Rsim scenario object
# 3) target_F: target F rates for managed stocks
# 4) target_B: corresponding target biomass for managed stocks
# 5) hcr: harvest control rule index (currently HCR = 1 status quo)
#
# Assumed to be available in the global environment:
#   B_equil, F_equil
#   hind_years, fore_years
#   groundfish, ssl_sp, nonssl_sp, managed_sp, medLH(medium-lived groundfish)
#
# cap: logical. Reserved for a future total-catch cap / apportionment rule
#      (analogous to ATTACH in the EBS model). Set to FALSE for now.
# ---------------------------------------------------------------------------- #


effort_to_F<-function(bal,scene){
  
}

start_biomass <- function(rsim){
  return(rsim$out_Biomass[1, 2:(dim(rsim$out_Biomass)[2])])
}

end_biomass <- function(rsim){
  return(rsim$out_Biomass[dim(rsim$out_Biomass)[1], 2:(dim(rsim$out_Biomass)[2])])
}

end_SSB <- function(rsim){
  return(rsim$out_SSB[dim(rsim$out_SSB)[1], 2:(dim(rsim$out_SSB)[2])])
}

start_catch <- function(rsim){
  return(rsim$out_Catch[1, 2:(dim(rsim$out_Catch)[2])])
}

end_catch <- function(rsim){
  return(rsim$out_Catch[dim(rsim$out_Catch)[1], 2:(dim(rsim$out_Catch)[2])])
}

living_groups   <- function(bal){return(bal$Group[1:bal$NUM_LIVING])}
detrital_groups <- function(bal){return(bal$Group[(bal$NUM_LIVING+1):(bal$NUM_LIVING+bal$NUM_DEAD)])}
gear_groups     <- function(bal){return(bal$Group[(bal$NUM_LIVING+bal$NUM_DEAD+1):(bal$NUM_GROUPS)])}

rsim.plotly<-function(rsim,species){
  ptxt<-NULL
  for (c in species){
    if (is.null(ptxt)){
      ptxt <- paste("p <- plot_ly(x=as.numeric(rownames(rsim$annual_Biomass)), y=rsim$annual_Biomass[,'",c,"'], name='",c,"', type = 'scatter', mode = 'lines') %>%", sep="") 
    }else {ptxt <- paste(ptxt," add_trace(y=rsim$annual_Biomass[,'",c,"'], name='",c,"', type = 'scatter', mode = 'lines') %>%",sep="")}
  }
  ptxt<-substr(ptxt,1,nchar(ptxt)-4)
  eval(parse(text = ptxt))
  return(p)       
}

# ---------------------------------------------------------------------------- #

a3_hcr_clim_sim <- function(ssps, scene, target_F, target_B, cap= FALSE, hcr=1) {
  cap <- cap
  if (cap == TRUE)  {
    cat("cap ON\n")
    flush.console()
  }
  if (cap == FALSE) {
    cat("cap OFF\n")
    flush.console()
  }
  
  F_target <- target_F
  B_target <- target_B[managed_sp, "Btarget_SQ"]
  
  
 
  
  # ========================================================================== #
  # Fisheries projections begin here
  
  # Run hindcast years
  run.hind  <- rsim.run(scene, method = 'AB', years = hind_years)
  goaclim.sim <- run.hind
  
  # Set up for observation error ---------- ---------- ---------- ---------- - #
  # object to store epsilon
  epsilon_mat            <- matrix(nrow = (length(fore_years)), ncol = (length(groundfish)))
  row.names(epsilon_mat) <- fore_years
  colnames(epsilon_mat)  <- groundfish
  # set autocorrelation and standard deviation of biomass estimates
  # long-lived LH params from Wiedenmann et al. (2015)
  phi   <- 0.89
  sigma_obs <- 0.34
  # medium-lived LH params from Wiedenmann et al. (2015) for pollock_adu, atka, and octopus
  phi_med   <- 0.84
  sigma_obs_med <- 0.31
  
  # object to store "stock status"
  ss_mat            <- epsilon_mat
  row.names(ss_mat) <- fore_years
  colnames(ss_mat)  <- groundfish
  # object to store ABC
  abc_mat            <- epsilon_mat
  row.names(abc_mat) <- fore_years
  colnames(abc_mat)  <- groundfish
  # object to store TAC
  # tac_mat            <- matrix(nrow=(length(fore_years)), ncol=(length(groundfish)))
  # row.names(tac_mat) <- fore_years
  # colnames(tac_mat)  <- groundfish
  medLH <- intersect(medLH, groundfish)
  longLH <- groundfish[!groundfish %in% medLH]
  
  # State variable to track alpha across years for HCR 2 (Rebuilding)
  alpha_state <- rep(0.05, length(managed_sp))
  names(alpha_state) <- managed_sp
  
  
  # 2. Forecast Loop
  
  # projections with HCRs and maybe ATTACH like CAP -------------------------- #
  # else {
  for (yr in fore_years) {
    cat(yr, ssps, "cap = ", cap, "hcr = ", hcr, "\n")
    flush.console()
    Ftarget    <- F_equil
    Btarget    <- B_equil

    # A. WGOA Observation Error & Assessment ####
    
    yr_char <- as.character(yr)
    prev_yr <- as.character(yr - 1)
    
    # "Assessment" of stock status with error
    # calculate epsilon
    if (yr == fore_years[1]) {
      # long-lived groundfish
      epsilon_mat[yr_char, longLH] <- rnorm(length(longLH), 0, sigma_obs)
      # medium-lived groundfish
      epsilon_mat[yr_char, medLH] <-rnorm(length(medLH), 0, sigma_obs_med)
    }else {
      # long-lived groundfish
      epsilon_mat[yr_char, longLH] <- phi * epsilon_mat[prev_yr, longLH] + sqrt(1 - phi^2) * 
        rnorm(length(longLH), 0, sigma_obs)
      # medium-lived groundfish
      epsilon_mat[yr_char, medLH] <- phi_med * epsilon_mat[prev_yr, medLH] +
        sqrt(1 - phi_med^2) * rnorm(length(medLH), 0, sigma_obs_med)
    }
    # do the assessment
    # catch_bio: total biomass × obs error — always used for C_ABC (catch = F × total_bio)
    # assessment: SSB substituted for SSB stocks × obs error — used only for Bratio
    # Keeping them separate prevents catch being computed from SSB instead of total biomass.
    catch_bio              <- end_biomass(goaclim.sim)
    assessment             <- catch_bio
    assessment[ssb_stocks] <- end_SSB(goaclim.sim)[ssb_stocks]

    groundfish_obs_error         <- exp(epsilon_mat[yr_char, groundfish] - 0.5 * sigma_obs^2)
    groundfish_obs_error[medLH]  <- exp(epsilon_mat[yr_char, medLH] - 0.5 * sigma_obs_med^2)

    catch_bio[groundfish]  <- catch_bio[groundfish]  * groundfish_obs_error  # for C_ABC
    assessment[groundfish] <- assessment[groundfish] * groundfish_obs_error  # for Bratio
    # store assessed stock status
    ss_mat[yr_char, ] <- assessment[groundfish]
    
    # ---------------------------------------------------------------------- #
    # B. Set Targets & Thresholds based on HCR ####
    # ---------------------------------------------------------------------- #
    
    if(hcr %in% c(1,3)){
      Ftarget[managed_sp] <- if(hcr==1) F_target else F50
      Btarget[managed_sp] <- if(hcr==1) B_target else B50
    }else{
      Ftarget[managed_sp] <- F40[managed_sp]
      Btarget[managed_sp] <- B40[managed_sp]
    }
    # For SSB stocks, override Btarget with the SSB-based reference point
    # to match the SSB-based assessment (line 180). Without this,
    # Bratio = SSB / B40_total, overstating stock status by 1.2-4.7x
    # (worst case: pacific_cod at 4.7x, pollock at 1.4x).
    ssb_in_managed <- intersect(ssb_stocks, managed_sp)
    
    if (hcr == 1) {
      # Assuming HCR 1 uses status quo; make sure "Btarget_SQ_SSB" exists in your target_B, 
      # or default it to B40_SSB
      Btarget[ssb_in_managed] <- target_B[ssb_in_managed, "B40_SSB"] 
    } else if (hcr == 3) {
      Btarget[ssb_in_managed] <- target_B[ssb_in_managed, "B50_SSB"]
    } else {
      Btarget[ssb_in_managed] <- target_B[ssb_in_managed, "B40_SSB"]
    }

    # Define beta cutoffs (directed fishing closures)
    beta_vec <- rep(0.05, length(managed_sp))
    names(beta_vec) <- managed_sp
    
    if(hcr==1){
      #FLAG ####
      beta_vec[ssl_sp] <- 0.50 #steller prey limits (check if that stands for WGOA)
      beta_vec[nonssl_sp] <- 0.05
      
    }else if(hcr==2){
      beta_vec[] <- 0.625
    }else if(hcr==3){
      beta_vec[] <- 0.40
    }else if(hcr %in% c(5,10)){
      beta_vec[] <- 0.50
    }
    
    gamma_val <- if(hcr==5) exp(0.1) else 0.1
    # ---------------------------------------------------------------------- #
    # C. HCR executions ####
    # ---------------------------------------------------------------------- #
    Bratio <- assessment[managed_sp]/Btarget[managed_sp]
    F_ABC <- Ftarget
    
    # Update state-dependent alpha for HCR 2
    if (hcr == 2) {
      for (sp in managed_sp) {
        if (Bratio[sp] < 0.625) {
          alpha_state[sp] <- 0.25
        } else if (Bratio[sp] >= 1.0) {
          alpha_state[sp] <- 0.05
        }
      }
    } else {
      alpha_state[] <- 0.05
    }
    
    # Vectorized / looped calculation of F_ABC based on ACLIM_HCR tiers
    for (sp in managed_sp) {
      b_ratio <- Bratio[sp]
      f_lim   <- Ftarget[sp]
      alpha   <- alpha_state[sp]
      limit   <- max(alpha, beta_vec[sp]) # Strict threshold
      
      if (b_ratio <= limit) {
        # Closed to directed fishing
        F_ABC[sp] <- 0.0
        
      } else if (b_ratio < 1.0) {
        # Rebuilding slope
        F_ABC[sp] <- f_lim * ((b_ratio - alpha) / (1.0 - alpha))
        
      } else {
        # Healthy status (> 1.0) dictates strategy
        if (hcr %in% c(1, 2, 3)) {
          F_ABC[sp] <- f_lim
          
        } else if (hcr == 5) {
          F_ABC[sp] <- f_lim * exp(-gamma_val * (b_ratio - 1))
          
        } else if (hcr == 10) {
          if (b_ratio >= (1 + gamma_val)) {
            F_ABC[sp] <- f_lim / (b_ratio / (1 + gamma_val))
          } else {
            F_ABC[sp] <- f_lim
          }
        }
      }
    }
    
    C_ABC <- F_ABC * catch_bio
    abc_mat[yr_char, ] <- C_ABC[groundfish]
    
    
    # ---------------------------------------------------------------------- #
    # D. Total-catch cap / apportionment (placeholder for WGOA implementation) ####
    # ---------------------------------------------------------------------- #
    
    if (cap == TRUE) {
      # TODO: implement WGOA-specific cap / apportionment rule here.
      # C_ABC should be modified in place before being written to ForcedCatch.
      stop("cap = TRUE is not yet implemented for the WGOA model.")
    }
    
    # Add resulting fishing to current year's CATCH matrix
    # scene <- adjust.fishing(scene, "CATCH", names(C_ABC), yr, value = C_ABC)
    scene$fishing$ForcedCatch[as.character(yr), names(C_ABC)] <- C_ABC
    # print("scene$fishing$ForcedCatch[as.character(yr), names(C_ABC)]")
    # print(scene$fishing$ForcedCatch[as.character(yr), names(C_ABC)])
    # Run an rsim year
    goaclim.sim <- rsim.step(scene, goaclim.sim, method = 'AB', year.end = yr)
  }
  # return(goaclim.sim)
  # }
  
  # outputs ------------------------------------------------------------------ #
  sim_out <- list(goaclim.sim, ss_mat, abc_mat)
  names(sim_out) <- c("sim_out", "ss_mat", "abc_mat")
  
  return(sim_out)
  
  # return(goaclim.sim)
}
