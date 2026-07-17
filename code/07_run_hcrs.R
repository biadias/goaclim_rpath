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

# ---------------------------------------------------------------------------- #
# Read reference point outputs from 06_run_btarget_ftarget_v1.R ####
# ---------------------------------------------------------------------------- #

bftarget_dir <- "data/bftarget"

# Biomass reference points (B0, B40, B40_SSB, Btarget_SQ, Blim, …)
target_bio <- read.csv(
  file.path(bftarget_dir, "B_target_WGOA_GFDL_persist.csv"),
  row.names = 1
)

# Climate-informed Ftarget (one column per scenario; "persist" for now)
Ftarg_matrix <- read.csv(
  file.path(bftarget_dir, "Ftarget_WGOA_GFDL_persist.csv"),
  row.names = 1
)

# Named vector of Ftarget for the persistence scenario
F_target_persist <- Ftarg_matrix[, "persist"]

# ---------------------------------------------------------------------------- #

a3_hcr_clim_sim <- function(ssps, scene, target_F, target_B, cap, hcr) {
  cap <- cap
  if (cap == TRUE)  {
    cat("cap ON\n")
    flush.console()
  }
  if (cap == FALSE) {
    cat("cap OFF\n")
    flush.console()
  }
  
  hcr <- hcr
  print(hcr)
  
  print(ssps)
  
  scene <- scene
  
  F_target <- target_F
  B_target <- target_B[managed_sp, "Btarget_SQ"] #; print(B_target)
  # Last year of projection period
  pyrl <- dim(scene$fishing$ForcedFRate)[1]
  
  # ========================================================================== #
  # Fisheries projections begin here
  
  # Run hindcast years
  run.hind  <- rsim.run(scene, method = 'AB', years = hind_years)
  aclim.sim <- run.hind
  
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
  ss_mat            <- matrix(nrow = (length(fore_years)), ncol = (length(groundfish)))
  row.names(ss_mat) <- fore_years
  colnames(ss_mat)  <- groundfish
  # object to store ABC
  abc_mat            <- matrix(nrow = (length(fore_years)), ncol = (length(groundfish)))
  row.names(abc_mat) <- fore_years
  colnames(abc_mat)  <- groundfish
  # object to store TAC
  # tac_mat            <- matrix(nrow=(length(fore_years)), ncol=(length(groundfish)))
  # row.names(tac_mat) <- fore_years
  # colnames(tac_mat)  <- groundfish
  
  
  # non-ATTACH and non-HCR scenarios ----------------------------------------- #
  # if (attach == "none") {
  #   # F_2015 scenario
  #   if (hcr == "static") {
  #     scene_2015 <- scene
  #     # set F in projection period to F_2015
  #     for (i in 32:109) {
  #       scene_2015$fishing$ForcedFRate[i, 1:77] <- F_2015
  #     }
  #     aclim.sim <- rsim.run(scene_2015, method = 'AB', years = all_years)
  #   }
  #   # F_0 scenario
  #   if (hcr == 0) {
  #     scene_0 <- scene
  #     # set F in projection period to F_2015
  #     for (i in 32:109) {
  #       scene_0$fishing$ForcedFRate[i, 1:77] <- F_zero
  #     }
  #     aclim.sim <- rsim.run(scene_0, method = 'AB', years = all_years)
  #   }
  #   # F_target (i.e., fish at F_target regardless of stock status)
  #   if (hcr == 99) {
  #     scene_99 <- scene
  #     # set F to F_target for managed stocks
  #     for (i in 32:109) {
  #       scene_99$fishing$ForcedFRate[i, managed_sp] <- F_target[managed_sp]
  #     }
  #     aclim.sim <- rsim.run(scene_99, method = 'AB', years = all_years)
  #   }
  # }
  
  # projections with HCRs and maybe ATTACH like CAP -------------------------- #
  # else {
  for (yr in fore_years) {
    cat(yr, ssps, "cap = ", cap, "hcr = ", hcr, "\n")
    flush.console()
    Ftarget    <- F_equil
    Btarget    <- B_equil
    
    # "Assessment" of stock status with error ---------- ---------- --------
    # calculate epsilon
    if (yr == 2022) {
      # long-lived groundfish
      epsilon_mat[as.character(yr), groundfish[!groundfish %in% medLH]] <-
        rnorm(length(groundfish[!groundfish %in% medLH]), 0, sigma_obs)
      # medium-lived groundfish
      epsilon_mat[as.character(yr), medLH] <-
        rnorm(length(medLH), 0, sigma_obs_med)
    }
    else {
      # long-lived groundfish
      epsilon_mat[as.character(yr), groundfish[!groundfish %in% medLH]] <-
        phi * epsilon_mat[as.character(yr - 1), groundfish[!groundfish %in% medLH]] +
        sqrt(1 - phi^2) * rnorm(length(groundfish[!groundfish %in% medLH]), 0, sigma_obs)
      # medium-lived groundfish
      epsilon_mat[as.character(yr), medLH] <-
        phi_med * epsilon_mat[as.character(yr - 1), medLH] +
        sqrt(1 - phi_med^2) * rnorm(length(medLH), 0, sigma_obs_med)
    }
    # do the assessment
    assessment             <- end_biomass(aclim.sim)
    # print("assessment"); print(assessment[crab_sp])
    # assessment[ssb_stocks] <- end_SSB(aclim.sim)[ssb_stocks]  ; print(assessment[ssb_stocks])
    # groundfish_obs_error <- exp(epsilon_mat[as.character(yr), groundfish] - 0.5 * sigma_obs^2) #; print(groundfish_obs_error)
    groundfish_obs_error <- exp(epsilon_mat[as.character(yr), groundfish] - 0.5 * sigma_obs^2) #; print(groundfish_obs_error)
    # medium-lived groundfish
    groundfish_obs_error[medLH] <- exp(epsilon_mat[as.character(yr), medLH] - 0.5 * sigma_obs_med^2) #; print(groundfish_obs_error)
    assessment[groundfish] <- assessment[groundfish] * groundfish_obs_error #; print(assessment[groundfish])
    # store assessed stock status
    ss_mat[as.character(yr), ] <- assessment[groundfish]
    
    # Set params for the different HCRs ---------- ---------- ---------- ---
    if (hcr == 2 & yr > 2022) {
      print(falpha_ssl)
    }
    
    if (hcr == 1) {
      Ftarget[managed_sp] <- F_target # this is F_target for status quo
      Btarget[managed_sp] <- B_target #; print("Btarget"); print(Btarget[crab_sp])# B_target for status quo
      falpha_crab         <- 0.1
      fbeta_crab          <- 0.25 # Beta is 25% of B target (B35)
      falpha_ssl          <- 0.05
      fbeta_ssl           <- 0.50 # B20 (i.e., half of B40)
      falpha_nonssl       <- 0.05
      fbeta_nonssl        <- 0.05
    }
    
    if (hcr == 2) {
      Ftarget[managed_sp] <- F40
      Btarget[managed_sp] <- B40
      fbeta_crab          <- 0.625 # cutoff is B25 (i.e., 62.5% of B40)
      fbeta_ssl           <- 0.625 # cutoff is B25 (i.e., 62.5% of B40)
      fbeta_nonssl        <- 0.625 # cutoff is B25 (i.e., 62.5% of B40)
      
      if (yr == 2022) {
        # Ftarget[managed_sp] <- F40
        # Btarget[managed_sp] <- B40
        falpha_crab     <- rep(0.05, length(crab_sp))
        names(falpha_crab) <- crab_sp
        # fbeta_crab      <- 0.625 # cutoff is B25 (i.e., 62.5% of B40)
        # for(i in 1:length(ssl_sp)) {
        #   if(yr > 2022) {
        #     if(falpha_ssl[ssl_sp[i]] < 0.25) {falpha_ssl[ssl_sp[i]] == 0.05}
        # }
        falpha_ssl     <- rep(0.05, length(ssl_sp))
        names(falpha_ssl) <- ssl_sp #; print(falpha_ssl)
        # fbeta_ssl      <- 0.625 # cutoff is B25 (i.e., 62.5% of B40)
        falpha_nonssl     <- rep(0.05, length(nonssl_sp))
        names(falpha_nonssl) <- nonssl_sp
        # fbeta_nonssl      <- 0.625 # cutoff is B25 (i.e., 62.5% of B40)
      }
      # check if alpha is re-building alpha. if so keep "big" alpha until
      # rebuilt, i.e., Bratio >= 1
      else {
        for (i in 1:length(crab_sp)) {
          ifelse(falpha_crab[i] >  0.05,
                 falpha_crab[i] <- falpha_crab[i],
                 falpha_crab[i] <- 0.05)
        }
        for (i in 1:length(ssl_sp)) {
          ifelse(falpha_ssl[i] >  0.05,
                 falpha_ssl[i] <- falpha_ssl[i],
                 falpha_ssl[i] <- 0.05)
        }
        for (i in 1:length(nonssl_sp)) {
          ifelse(
            falpha_nonssl[i] >  0.05,
            falpha_nonssl[i] <- falpha_nonssl[i],
            falpha_nonssl[i] <- 0.05
          )
        }
        # print(falpha_ssl)
      }
    }
    
    if (hcr == 3) {
      Ftarget[managed_sp] <- F50
      Btarget[managed_sp] <- B50
      falpha_crab         <- 0.05
      fbeta_crab          <- 0.4 # cutoff is B20 (i.e., 40% of B50)
      falpha_ssl          <- 0.05
      fbeta_ssl           <- 0.4 # cutoff is B20 (i.e., 40% of B50)
      falpha_nonssl       <- 0.05
      fbeta_nonssl        <- 0.4 # cutoff is B20 (i.e., 40% of B50)
    }
    
    if (hcr == 5) {
      Ftarget[managed_sp] <- F40
      Btarget[managed_sp] <- B40
      falpha_crab         <- 0.05
      fbeta_crab          <- 0.50
      falpha_ssl          <- 0.05
      fbeta_ssl           <- 0.50 # B20 (i.e., half of B40)
      falpha_nonssl       <- 0.05
      fbeta_nonssl        <- 0.50
      fgamma_ssl          <- exp(0.1)
      fgamma_nonssl       <- exp(0.1)
      fgamma_crab         <- exp(0.1)
    }
    
    if (hcr == 10) {
      Ftarget[managed_sp] <- F40
      Btarget[managed_sp] <- B40
      falpha_crab         <- 0.05
      fbeta_crab          <- 0.50
      falpha_ssl          <- 0.05
      fbeta_ssl           <- 0.50 # B20 (i.e., half of B40)
      falpha_nonssl       <- 0.05
      fbeta_nonssl        <- 0.50
      fgamma_ssl          <- 0.1
      fgamma_nonssl       <- 0.1
      fgamma_crab         <- 0.1
    }
    
    # Create Species-indexed vector F_ABC:
    F_ABC <- Ftarget
    C_ABC <- F_ABC * end_biomass(aclim.sim)
    
    # ---------------------------------------------------------------------- #
    # Institute control rule for crab species ---------- ---------- -------- #
    sp <- crab_sp
    # assessment <- end_biomass(aclim.sim) #* rnorm(length(end_biomass(aclim.sim)), mean=1, sd=0.2)
    Bratio     <- assessment / Btarget
    # Bstatus[tyr,sp] <- Bratio[sp]
    if (hcr == 2) {
      for (i in 1:length(crab_sp)) {
        if (Bratio[crab_sp[i]] <  fbeta_crab) {
          falpha_crab[crab_sp[i]] <- 0.25
        }
        if (Bratio[crab_sp[i]] >= 1)          {
          falpha_crab[crab_sp[i]] <- 0.05
        }
      }
      F_ABC[sp]  <- ifelse((Bratio[sp] <= falpha_crab[sp]) |
                             (Bratio[sp] <= fbeta_crab),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_crab[sp]) / (1.0 - falpha_crab[sp]),
                             Ftarget[sp]
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 1 | hcr == 3) {
      F_ABC[sp]  <- ifelse((Bratio[sp] <= falpha_crab) |
                             (Bratio[sp] <= fbeta_crab),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_crab) / (1.0 - falpha_crab),
                             Ftarget[sp]
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 5) {
      F_ABC[sp] <-  ifelse((Bratio[sp] <= falpha_crab) |
                             (Bratio[sp] <= fbeta_crab),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_crab) / (1.0 - falpha_crab),
                             Ftarget[sp] * (exp(-fgamma_crab * (Bratio[sp] - 1)))
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 10) {
      F_ABC[sp] <-  ifelse((Bratio[sp] <= falpha_crab) |
                             (Bratio[sp] <= fbeta_crab),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_crab) / (1.0 - falpha_crab),
                             ifelse(
                               Bratio[sp] > 1.0 & Bratio[sp] < (1 + fgamma_crab),
                               Ftarget[sp],
                               Ftarget[sp] / (Bratio[sp] * (1 / (
                                 1 + fgamma_crab
                               )))
                             )
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    # print("ABC"); print(C_ABC[crab_sp])
    # ---------------------------------------------------------------------- #
    # Institute control rule for ssl prey species (pollock, p.cod, and atka)
    sp <- ssl_sp
    # assessment <- end_biomass(aclim.sim) #* rnorm(length(end_biomass(aclim.sim)), mean=1, sd=0.2)
    # Bratio     <- assessment / Btarget
    if (hcr == 2) {
      for (i in 1:length(ssl_sp)) {
        if (Bratio[ssl_sp[i]] <  fbeta_ssl) {
          falpha_ssl[ssl_sp[i]] <- 0.25
        }
        if (Bratio[ssl_sp[i]] >= 1)         {
          falpha_ssl[ssl_sp[i]] <- 0.05
        }
      }
      F_ABC[sp]  <- ifelse((Bratio[sp] <= falpha_ssl[sp]) |
                             (Bratio[sp] <= fbeta_ssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_ssl[sp]) / (1.0 - falpha_ssl[sp]),
                             Ftarget[sp]
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 1 | hcr == 3) {
      F_ABC[sp]  <- ifelse((Bratio[sp] <= falpha_ssl) |
                             (Bratio[sp] <= fbeta_ssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_ssl) / (1.0 - falpha_ssl),
                             Ftarget[sp]
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp] #; print(F_ABC[sp])
    }
    
    if (hcr == 5) {
      F_ABC[sp] <-  ifelse((Bratio[sp] <= falpha_ssl) |
                             (Bratio[sp] <= fbeta_ssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_ssl) / (1.0 - falpha_ssl),
                             Ftarget[sp] * (exp(-fgamma_ssl * (Bratio[sp] - 1)))
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 10) {
      F_ABC[sp] <-  ifelse((Bratio[sp] <= falpha_ssl) |
                             (Bratio[sp] <= fbeta_ssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_ssl) / (1.0 - falpha_ssl),
                             ifelse(
                               Bratio[sp] > 1.0 & Bratio[sp] < (1 + fgamma_ssl),
                               Ftarget[sp],
                               Ftarget[sp] / (Bratio[sp] * (1 / (1 + fgamma_ssl)))
                             )
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    # ---------------------------------------------------------------------- #
    # Institute control rule for capped species (Federal Groundfish) --------
    sp <- nonssl_sp
    # assessment <- end_biomass(aclim.sim) #* rnorm(length(end_biomass(aclim.sim)), mean=1, sd=0.2)
    # Bratio     <- assessment / Btarget
    if (hcr == 2) {
      for (i in 1:length(nonssl_sp)) {
        if (Bratio[nonssl_sp[i]] <  fbeta_nonssl) {
          falpha_nonssl[nonssl_sp[i]] <- 0.25
        }
        if (Bratio[nonssl_sp[i]] >= 1)            {
          falpha_nonssl[nonssl_sp[i]] <- 0.05
        }
      }
      F_ABC[sp]  <- ifelse((Bratio[sp] <= falpha_nonssl[sp]) |
                             (Bratio[sp] <= fbeta_nonssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_nonssl[sp]) / (1.0 - falpha_nonssl[sp]),
                             Ftarget[sp]
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 1 | hcr == 3) {
      F_ABC[sp]  <- ifelse((Bratio[sp] <= falpha_nonssl) |
                             (Bratio[sp] <= fbeta_nonssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_nonssl) / (1.0 - falpha_nonssl),
                             Ftarget[sp]
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 5) {
      F_ABC[sp] <-  ifelse((Bratio[sp] <= falpha_nonssl) |
                             (Bratio[sp] <= fbeta_nonssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_nonssl) / (1.0 - falpha_nonssl),
                             Ftarget[sp] * (exp(-fgamma_nonssl * (Bratio[sp] - 1)))
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    if (hcr == 10) {
      F_ABC[sp] <-  ifelse((Bratio[sp] <= falpha_nonssl) |
                             (Bratio[sp] <= fbeta_nonssl),
                           0.0,
                           ifelse(
                             Bratio[sp] < 1.0,
                             Ftarget[sp] * (Bratio[sp] - falpha_nonssl) / (1.0 - falpha_nonssl),
                             ifelse(
                               Bratio[sp] > 1.0 & Bratio[sp] < (1 + fgamma_nonssl),
                               Ftarget[sp],
                               Ftarget[sp] / (Bratio[sp] * (1 / (
                                 1 + fgamma_nonssl
                               )))
                             )
                           )
      )
      C_ABC[sp] <- F_ABC[sp] * assessment[sp]
    }
    
    
    # store ABCs
    abc_mat[as.character(yr), ] <- C_ABC[groundfish]
    # for Ebett and Stephani's pollock no fishing sims:
    # C_ABC["pollock_adu"] <- 0 #; print(C_ABC[groundfish])
    
    # ---------------------------------------------------------------------- #
    # Total-catch cap / apportionment (placeholder for WGOA implementation)
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
    aclim.sim <- rsim.step(scene, aclim.sim, method = 'AB', year.end = yr)
  }
  # return(aclim.sim)
  # }
  
  # outputs ------------------------------------------------------------------ #
  sim_out <- list(aclim.sim, ss_mat, abc_mat)
  names(sim_out) <- c("sim_out", "ss_mat", "abc_mat")
  
  return(sim_out)
  
  # return(aclim.sim)
}