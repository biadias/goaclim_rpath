#' ---------------------------------------------------------------------------- #
#' AUTHORS: Bia Dias
#' ORIGINAL CODE: George (Andy) Whitehouse
#' AFFILIATIONS: CICOES University of Washington
#' E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#' DATE: 10 March 2026
#' @description
#' 
#' code/Function_plot_sim.R
#' Purpose: Plot Relative Biomass Projections
#' ---------------------------------------------------------------------------- #
#' @param sim_list A named list of rsim model outputs (the result of `rsim.run()`).
#' @param species_to_plot Character vector of species/groups you want to plot.
#' @param start_year Numeric. The starting year of your simulation. Default is 1991.
#' @param variable Character. Which variable to plot. Options are "Biomass" (default) or "Catch".
#'
#' @return A ggplot object comparing scenarios faceted by species.
#' @export
#' 
plot_scenario_comparison <- function(sim_list, species_to_plot, start_year = 1991, variable = "Biomass") {
  
  # Determine which output matrix to extract based on user choice
  var_name <- if(variable == "Biomass") "out_Biomass" else "out_Catch"
  
  # Extract and bind data from all scenarios
  plot_data <- do.call(rbind, lapply(names(sim_list), function(scen_name) {
    
    # Get the specific simulation output
    sim <- sim_list[[scen_name]]
    
    # Check if the output variable exists
    if(!var_name %in% names(sim)) {
      stop(paste("Variable", var_name, "not found in the simulation output."))
    }
    
    # Subset the matrix to just the species we care about
    mat <- sim[[var_name]][, species_to_plot, drop = FALSE]
    
    # Convert to dataframe and add time/scenario info
    df <- as.data.frame(mat)
    df$Month <- 1:nrow(df)
    df$Year <- start_year + (df$Month - 1) / 12  # Convert months to fractional years
    df$Scenario <- scen_name
    
    # Melt data from wide to long format for ggplot
    df_long <- pivot_longer(df, 
                            cols = all_of(species_to_plot), 
                            names_to = "Species", 
                            values_to = "Value")
    
    return(df_long)
  }))
  
  # Ensure the order of species in the plot matches the vector provided
  plot_data$Species <- factor(plot_data$Species, levels = species_to_plot)

  
  # Create the ggplot
  p <- ggplot(plot_data, aes(x = Year, y = Value, color = Scenario)) +
    geom_line(size = 1) +
    # Free Y scales so each species' plot fits its own biomass magnitude
    facet_wrap(~ Species, scales = "free_y") + 
    scale_colour_manual(values = c("SSP 585" = "#d55e00",
                                   "SSP 126" = "#009e73", 
                                   "SSP 245" = "#e69f00", 
                                   "Persist" = "#000000")) +
    # Add a vertical dashed line where the projection begins (end of 2020)
    geom_vline(xintercept = 2021, linetype = "dashed", color = "gray50") + 
    theme_bw() +
    labs(
      #title = paste("Comparison of", variable, "Across Climate Scenarios"),
      #subtitle = "Dashed line indicates start of projection (2021)",
      x = "Year",
      y = paste(variable, "(t/km^2)"),
      color = "Scenario"
    ) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold")
    )
  
  return(p)
}