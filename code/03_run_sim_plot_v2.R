# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 03 March 2026
#
# code/02_run_sim.R
# Purpose: Running the F_clim_sim_scene function to set up the scenarios and then 
# running the simulations for the best model scenario (GOA_fit_results_59M04par.rds) 
# for the GFDL SSP 126, SSP 425, SSP 585 scenarios, as well as a persistence scenario.
# ---------------------------------------------------------------------------- #
source("code/02_run_sim_pcod.R") 
source("code/Function_plot_scenario_comparison.R") 


list_ssps_pp <- list( "F0 persist Pcod rec PP" = fore_pers_buf_pcod_f0, "F0 126 Pcod rec PP" = forecast_126_buf_pcod_f0, "F0 245 Pcod rec PP" = forecast_245_buf_pcod_f0, "F0 585 Pcod rec PP" = forecast_585_buf_pcod_f0, "Fmean persist Pcod rec PP"= fore_pers_buf_pcod_fmean, "Fmean 126 Pcod rec PP" = forecast_126_buf_pcod_fmean, "Fmean 245 Pcod rec PP" = forecast_245_buf_pcod_fmean, "Fmean 585 Pcod rec PP" = forecast_585_buf_pcod_fmean )

list_ssps_pcod<- list( "F0 persist Pcod rec" = fore_pers_pcod_f0, "F0 126 Pcod rec" = forecast_126_pcod_f0, "F0 245 Pcod rec" = forecast_245_pcod_f0, "F0 585 Pcod rec" = forecast_585_pcod_f0, "Fmean persist Pcod rec"= fore_pers_pcod_fmean, "Fmean 126 Pcod rec" = forecast_126_pcod_fmean, "Fmean 245 Pcod rec" = forecast_245_pcod_fmean, "Fmean 585 Pcod rec" = forecast_585_pcod_fmean )

list_ssps <- list( "F0 persist Pcod" = fore_pers_f0, "F0 126 Pcod" = forecast_126_f0, "F0 245 Pcod" = forecast_245_f0, "F0 585 Pcod" = forecast_585_f0, "Fmean persist Pcod"= fore_pers_fmean, "Fmean 126 Pcod" = forecast_126_fmean, "Fmean 245 Pcod" = forecast_245_fmean, "Fmean 585 Pcod" = forecast_585_fmean )

# --------------------------------------------------------------------------- #
# PLOT — save 3-page PDF, one page per recruitment/buffer scenario
# --------------------------------------------------------------------------- #

species_list <- c(#"arrowtooth_flounder_adult", 
  "pacific_cod_adult", #"walleye_pollock_adult",
                  #"arrowtooth_flounder_juvenile", 
  "pacific_cod_juvenile")#, "walleye_pollock_juvenile")

cols <- c("persist" = "black", "126" = "#009E73", "245" = "#56B4E9", "585" = "#E69F00")
ltys <- c("F0" = 2, "Fmean" = 1)

scenario_lists <- list(
  "Base"           = list_ssps,
  "Pcod rec"       = list_ssps_pcod,
  "Buf + Pcod rec" = list_ssps_pp
)

get_ssp <- function(nm) {
  for (s in c("persist", "585", "245", "126")) {
    if (grepl(s, nm)) return(s)
  }
  NA_character_
}
get_fishing <- function(nm) {
  if (grepl("^Fmean", nm)) "Fmean" else if (grepl("^F0", nm)) "F0" else NA_character_
}

# Shared ylim per species across all list types
global_max <- setNames(
  sapply(species_list, function(sp) {
    max(sapply(unlist(scenario_lists, recursive = FALSE),
               function(run) max(run$annual_Biomass[, sp], na.rm = TRUE)))
  }),
  species_list
)

out_path <- "C:/Users/biadias/Documents/goaclim_rpath/figures/biomass_trajectories_pcod.pdf"
pdf(out_path, width = 14, height = 6)

for (list_label in names(scenario_lists)) {
  lst <- scenario_lists[[list_label]]
  
  # Extract years from row names
  years <- as.numeric(rownames(lst[[1]]$annual_Biomass))
  
  par(mfrow = c(2, 3),
      mar   = c(4, 4, 3, 1),
      oma   = c(0, 0, 2, 0))
  
  for (sp in species_list) {
    
    plot(years, rep(NA_real_, length(years)),
         type = "n",
         #ylim = c(0, global_max[sp] * 1.1),
         xlab = "Year",
         ylab = "Biomass",
         main = sp)
    
    # Vertical line at last hindcast year
    abline(v = 2020, lty = 2, col = "grey40", lwd = 1.5)
    
    for (run_name in names(lst)) {
      ssp     <- get_ssp(run_name)
      fishing <- get_fishing(run_name)
      
      if (!is.na(ssp) && !is.na(fishing)) {
        lines(years, lst[[run_name]]$annual_Biomass[, sp],
              col = cols[ssp],
              lty = ltys[fishing],
              lwd = 2)
      }
    }
    
    if (sp == species_list[1]) {
      legend("bottomleft",
             legend = names(cols), col = cols,
             lty = 1, lwd = 2,
             title = "SSP", bty = "n", cex = 1.5)
      
      legend("bottomright",
             legend = names(ltys), col = "black",
             lty = ltys, lwd = 2,
             title = "Fishing", bty = "n", cex = 1.5)
    }
  }
  
  mtext(paste("Scenario:", list_label), outer = TRUE, cex = 1.3, font = 2)
}

dev.off()
message("Saved to: ", out_path)
