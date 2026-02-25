#------------------------------------------------------------------------------#
#REVIEWED: Bia Dias
#ORIGINAL AUTHORS: Andy Whitehouse
#AFFILIATIONS: CICOES University of Washington/ Alaska Fisheries Science Center
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov and andy.whitehouse@noaa.gov
#
# bioenergetic projections
#------------------------------------------------------------------------------#

source("GOA/wgoa_bioenergetics_code/bioenergetic_pars2.r")

bioen_params <- function(ssp) {
  

  # get the raw temperature projection output from ROMSNPZ scenario
  climate_file <- paste0("GOA/wgoa_data_rpath_fitting/", ssp,"_wide_WGOA_temp_1000",".csv")
  climate_proj <- read.csv(climate_file, row.names = NULL)
  bot_temp     <- as.matrix(climate_proj$btemp[493:1440])
  row.names(bot_temp) <- climate_proj$tstep[493:1440]
  
  # consumption multiplier
  cons_mult    <- paste0("tdc_", ssp, "_bt")
  cons_mult    <- matrix(nrow=dim(bot_temp)[1], ncol=length(bioen_sp_noceph))
  colnames(cons_mult) <- bioen_sp_noceph
  for(i in 1:dim(bot_temp)[1]){
    for(j in bioen_sp_noceph){
      cons_mult[i,j] <- rc_scaled_b[sprintf("%.2f", round(bot_temp[i], digits=2)),j]
    }
  }
  row.names(cons_mult) <- row.names(bot_temp)
  
  # respiration modifiers
  resp_mult    <- paste0("tdr_", ssp, "_bt")
  resp_mult    <- matrix(nrow=dim(bot_temp)[1], ncol=length(bioen_sp_noceph))
  colnames(resp_mult) <- bioen_sp_noceph 
  for(i in 1:dim(bot_temp)[1]){
    for(j in bioen_sp_noceph){
      resp_mult[i,j] <- ForcedActResp_b[sprintf("%.2f", round(bot_temp[i], digits=2)),j]
    }
  }
  row.names(resp_mult) <- row.names(bot_temp)
  
  # outputs
  bioen_parameters <- list(cons_mult, resp_mult)
  names(bioen_parameters) <- c(paste0("tdc_", ssp, "_bt"),
                               paste0("tdr_", ssp, "_bt"))

  return(bioen_parameters)

}
#-------------------------------------------------------------------------------#




