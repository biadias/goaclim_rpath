#------------------------------------------------------------------------------#
#REVIEWED: Bia Dias
#ORIGINAL AUTHORS: Andy Whitehouse
#AFFILIATIONS: CICOES University of Washington/ Alaska Fisheries Science Center
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov and andy.whitehouse@noaa.gov
#
# bioenergetic projections
#------------------------------------------------------------------------------#

source("code/bioenergetics_pars3.R")

bioen_params <- function(ssp) {
  ssp_clean <- gsub("ssp", "", ssp)
  
  # 1. Get the raw temperature projection output from ROMSNPZ scenario
  climate_file <- paste0("Rpath_fitting/GOA/wgoa_data_rpath_fitting/ssp", ssp_clean, "_wide_WGOA_temp_1000.csv")
  climate_proj <- read.csv(climate_file, row.names = NULL)
  
  # Extract BOTH Bottom and Surface temperatures for the projection window
  bot_temp <- as.matrix(climate_proj$btemp[493:1440])
  sur_temp <- as.matrix(climate_proj$stemp[493:1440])
  time_steps <- climate_proj$tstep[493:1440]
  
  n_steps <- length(time_steps)
  
  # 2. Initialize multiplier matrices
  cons_mult <- matrix(nrow = n_steps, ncol = length(bioen_sp_noceph))
  resp_mult <- matrix(nrow = n_steps, ncol = length(bioen_sp_noceph))
  
  colnames(cons_mult) <- bioen_sp_noceph
  colnames(resp_mult) <- bioen_sp_noceph 
  row.names(cons_mult) <- time_steps
  row.names(resp_mult) <- time_steps
  
  # 3. Populate matrices using Niche Logic and Temperature Clamping
  for(i in 1:n_steps) {
    for(j in bioen_sp_noceph) {
      
      # Flexible Habitat selection 
      habitat_type <- niche_lookup[j]
      
      # Determine if pelagic (flagdem == 0) or manually flag known pelagics like capelin
      is_pelagic <- (!is.na(habitat_type) && habitat_type == 0)
      
      # Pick Surface temp for Pelagic, else Bottom temp
      raw_temp <- ifelse(is_pelagic, sur_temp[i], bot_temp[i])
      
      # Handle potential NAs in climate projections
      if (is.na(raw_temp)) raw_temp <- 5.0 
      
      # Clamp to [-2, 30] to match your Kitchell/Arrhenius matrices
      target_temp <- max(-2, min(30, raw_temp))
      temp_str <- sprintf("%.2f", round(target_temp, digits=2))
      
      # Safe Assignment using the unified matrices
      if (temp_str %in% rownames(rc_scaled) && j %in% colnames(rc_scaled)) {
        cons_mult[i, j] <- rc_scaled[temp_str, j]
        resp_mult[i, j] <- ForcedActResp[temp_str, j]
      } else {
        # Neutral fallback if something goes completely wrong
        cons_mult[i, j] <- 1.0
        resp_mult[i, j] <- 1.0
      }
    }
  }
  
  # 4. Outputs
  bioen_parameters <- list(cons_mult, resp_mult)
  names(bioen_parameters) <- c(paste0("tdc_", ssp, "_niche"),
                               paste0("tdr_", ssp, "_niche"))
  
  return(bioen_parameters)
}
