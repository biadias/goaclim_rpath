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
source("code/02_run_sim.R") 
source("code/Function_plot_scenario_comparison.R") 
# ---------------------------------------------------------------------------- #
# 1. Run the simulations (this generates the actual output matrices)####
# ---------------------------------------------------------------------------- #
#This is saved in the 02_run_sim.R
all_scenes <- readRDS("results/diagnostics/all_scenes_cons_resp_ssps.rds")

forecast_gfdl_persist      <- rsim.run(scene_gfdl_persist, method = "AB", years = all_years)
forecast_gfdl_persist_cons <- rsim.run(all_scenes$persist_cons, method = "AB", years = all_years)
forecast_gfdl_persist_resp <- rsim.run(all_scenes$persist_resp, method = "AB", years = all_years)
forecast_gfdl_persist_none <- rsim.run(all_scenes$persist_none, method = "AB", years = all_years)
forecast_gfdl_126          <- rsim.run(scene_gfdl_126, method = "AB", years = all_years)
forecast_gfdl_126_cons     <- rsim.run(all_scenes$`126_cons`, method = "AB", years = all_years)
forecast_gfdl_126_resp     <- rsim.run(all_scenes$`126_resp`, method = "AB", years = all_years)
forecast_gfdl_126_none     <- rsim.run(all_scenes$`126_none`, method = "AB", years = all_years)
forecast_gfdl_245          <- rsim.run(scene_gfdl_245, method = "AB", years = all_years)
forecast_gfdl_245_cons     <- rsim.run(all_scenes$`245_cons`, method = "AB", years = all_years)
forecast_gfdl_245_resp     <- rsim.run(all_scenes$`245_resp`, method = "AB", years = all_years)
forecast_gfdl_245_none     <- rsim.run(all_scenes$`245_none`, method = "AB", years = all_years)
forecast_gfdl_585          <- rsim.run(scene_gfdl_585, method = "AB", years = all_years)
forecast_gfdl_585_cons     <- rsim.run(all_scenes$`585_cons`, method = "AB", years = all_years)
forecast_gfdl_585_resp     <- rsim.run(all_scenes$`585_resp`, method = "AB", years = all_years)
forecast_gfdl_585_none     <- rsim.run(all_scenes$`585_none`, method = "AB", years = all_years)


# ---------------------------------------------------------------------------- #
# 2. Plot ssps####
# ---------------------------------------------------------------------------- #
# The names here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios <- list(
  "Persist" = forecast_gfdl_persist,
  "SSP 126" = forecast_gfdl_126,
  "SSP 245" = forecast_gfdl_245,
  "SSP 585" = forecast_gfdl_585
)

#plot_spp <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
plot_arrow_prey <- c(
  "euphausiids",
  "pandalid_shrimp",
  "pacific_capelin",
  "pacific_sandlance",
  "walleye_pollock_adult",
  "walleye_pollock_juvenile",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult",
  "arrowtooth_flounder_juvenile"
)
plot_arrow_pred <- c(
  "steller_sea_lion",
  "longnose_skate",
  "pacific_halibut_adult",
  "arrowtooth_flounder_adult",
  "shallow_water_flatfish",
  "pacific_cod_adult"
)

plot_ground <- c(
  "arrowtooth_flounder_adult",
  "pacific_cod_adult",
  "walleye_pollock_adult")

plot_ground_juv <- c(
  "arrowtooth_flounder_juvenile",
  "pacific_cod_juvenile",
  "walleye_pollock_juvenile")

fg1 <- plot.species[1:12]
fg2 <- plot.species[13:24]
fg3 <- plot.species[25:36]
fg4 <- plot.species[37:48]
fg5 <- plot.species[49:60]
fg6 <- plot.species[61:72]
fg6 <- plot.species[73:84]

vector_list <- list(fg1, fg2, fg3, fg4, fg5, fg6)


for(i in 1:length(vector_list)){
  ssp_comps <- list()
  ssp_comps[[i]] <- plot_scenario_comparison(
    sim_list = my_scenarios, 
    species_to_plot = vector_list[[i]], 
    start_year = 1991, 
    variable = "Biomass"  # You can also change this to "Catch"
  )
  print(ssp_comps[[i]])
}

print(ssp_comps)

plot_scenario_comparison(
  sim_list = my_scenarios, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)
# ---------------------------------------------------------------------------- #
# 2. Plot ssp 126 res con ####
# ---------------------------------------------------------------------------- #

# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios_126 <- list(
  "SSP 126" = forecast_gfdl_126,
  "SSP 126 (consumption only)" = forecast_gfdl_126_cons,
  "SSP 126 (respiration only)" = forecast_gfdl_126_resp,
  "SSP 126 (no bioenergetic modifiers)" = forecast_gfdl_126_none

)

#plot_arrow_prey <- c("euphausiids", "pandalid_shrimp", "pacific_capelin","pacific_sandlance",
#                     "walleye_pollock_adult", "walleye_pollock_juvenile", "pacific_cod_adult",
#                     "arrowtooth_flounder_adult", "arrowtooth_flounder_juvenile")
#plot_arrow_pred <- c("steller_sea_lion",
#                     "longnose_skate",  
#                     "pacific_halibut_adult", "arrowtooth_flounder_adult", "shallow_water_flatfish", 
#                     "pacific_cod_adult"
#)

plot_ground <- c("arrowtooth_flounder_adult", "pacific_cod_adult","walleye_pollock_adult" )
ssp126 <- plot_scenario_comparison(
  sim_list = my_scenarios_126, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(ssp126)

# ---------------------------------------------------------------------------- #
# 2. Plot ssp 245 res con####
# ---------------------------------------------------------------------------- #

# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios_245 <- list(
  "SSP 245" = forecast_gfdl_245,
  "SSP 245 (consumption only)" = forecast_gfdl_245_cons,
  "SSP 245 (respiration only)" = forecast_gfdl_245_resp,
  "SSP 245 (no bioenergetic modifiers)" = forecast_gfdl_245_none
  
)

ssp245 <- plot_scenario_comparison(
  sim_list = my_scenarios_245, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(ssp245)

# ---------------------------------------------------------------------------- #
# 2. Plot ssp 585 res con####
# ---------------------------------------------------------------------------- #

# The names you use here (e.g., "SSP 126") will become the labels in the plot legend!
my_scenarios_585 <- list(
  "SSP 585" = forecast_gfdl_585,
  "SSP 585 (consumption only)" = forecast_gfdl_585_cons,
  "SSP 585 (respiration only)" = forecast_gfdl_585_resp,
  "SSP 585 (no bioenergetic modifiers)" = forecast_gfdl_585_none
  
)

ssp585 <- plot_scenario_comparison(
  sim_list = my_scenarios_585, 
  species_to_plot = plot_ground, 
  start_year = 1991, 
  variable = "Biomass"  # You can also change this to "Catch"
)

print(ssp585)




# --------------------------------------------------------------------------- #
# PLOT — one figure per recruitment/buffer scenario
# --------------------------------------------------------------------------- #
list_ssps_pp <- list( "F0 persist Pcod rec PP" = fore_pers_buf_pcod_f0, "F0 126 Pcod rec PP" = forecast_126_buf_pcod_f0, "F0 245 Pcod rec PP" = forecast_245_buf_pcod_f0, "F0 585 Pcod rec PP" = forecast_585_buf_pcod_f0, "Fmean persist Pcod rec PP"= fore_pers_buf_pcod_fmean, "Fmean 126 Pcod rec PP" = forecast_126_buf_pcod_fmean, "Fmean 245 Pcod rec PP" = forecast_245_buf_pcod_fmean, "Fmean 585 Pcod rec PP" = forecast_585_buf_pcod_fmean )

list_ssps_pcod<- list( "F0 persist Pcod rec" = fore_pers_pcod_f0, "F0 126 Pcod rec" = forecast_126_pcod_f0, "F0 245 Pcod rec" = forecast_245_pcod_f0, "F0 585 Pcod rec" = forecast_585_pcod_f0, "Fmean persist Pcod rec"= fore_pers_pcod_fmean, "Fmean 126 Pcod rec" = forecast_126_pcod_fmean, "Fmean 245 Pcod rec" = forecast_245_pcod_fmean, "Fmean 585 Pcod rec" = forecast_585_pcod_fmean )

list_ssps <- list( "F0 persist Pcod" = fore_pers_f0, "F0 126 Pcod" = forecast_126_f0, "F0 245 Pcod" = forecast_245_f0, "F0 585 Pcod" = forecast_585_f0, "Fmean persist Pcod"= fore_pers_fmean, "Fmean 126 Pcod" = forecast_126_fmean, "Fmean 245 Pcod" = forecast_245_fmean, "Fmean 585 Pcod" = forecast_585_fmean )


species_list <- c("arrowtooth_flounder_adult", "pacific_cod_adult", "walleye_pollock_adult",
                  "arrowtooth_flounder_juvenile", "pacific_cod_juvenile", "walleye_pollock_juvenile")

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

# One figure per list type
for (list_label in names(scenario_lists)) {
  lst <- scenario_lists[[list_label]]
  
  dev.new(width = 12, height = 8)
  par(mfrow = c(2, 3),
      mar   = c(4, 4, 3, 1),
      oma   = c(0, 0, 2, 0))
  
  for (sp in species_list) {
    n_years <- nrow(lst[[1]]$annual_Biomass)
    
    plot(seq_len(n_years), rep(NA_real_, n_years),
         type = "n",
         ylim = c(0, global_max[sp] * 1.1),
         xlab = "Year Index",
         ylab = "Biomass",
         main = sp)
    
    for (run_name in names(lst)) {
      ssp     <- get_ssp(run_name)
      fishing <- get_fishing(run_name)
      
      if (!is.na(ssp) && !is.na(fishing)) {
        lines(lst[[run_name]]$annual_Biomass[, sp],
              col = cols[ssp],
              lty = ltys[fishing],
              lwd = 2)
      }
    }
    
    # Legends on first panel only
    if (sp == species_list[1]) {
      legend("bottomleft",
             legend = names(cols), col = cols,
             lty = 1, lwd = 2,
             title = "SSP", bty = "n", cex = 0.75)
      
      legend("bottomright",
             legend = names(ltys), col = "black",
             lty = ltys, lwd = 2,
             title = "Fishing", bty = "n", cex = 0.75)
    }
  }
  
  mtext(paste("Scenario:", list_label), outer = TRUE, cex = 1.2, font = 2)
}

