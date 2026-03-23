#' ---------------------------------------------------------------------------- #
#' AUTHORS: Bia Dias
#' ORIGINAL CODE: Kerim Aydin
#' AFFILIATIONS: CICOES University of Washington
#' E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#' DATE: 17 March 2026
#' @description
#' 
#' code/Function_diagnostics_TEST.R
#' Purpose: functions for the diagnostics of the test suite
#' ---------------------------------------------------------------------------- #
#' comparing vulnerabilities from fitted model
#' 
#' Compares the vulnerability values from a fitted model to the original scenario 
#' parameters.
#' 
#' @param original_scene the baseline scene (not fitted)
#' @param fit_results the fitted results object containing the final VVs
#' @param plot logical, whether to generate a plot comparing original vs fitted VVs
#' @param scene_name Optional character string to include in plot title for context
#' 
#' @return a data frame containing the old vs new VV values for each predator-prey link
#' @importFrom dplyr mutate filter summarize n %>% 
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme_minimal theme element_text labs
#' @importFrom scales rescale_mid label_number squish
#' @export
#' 

get_vv_comparison <- function(original_scene, fit_results, plot = TRUE, scene_name=NULL, limit_type = c("log", "abs")) {
  
  limit_type <- match.arg(limit_type)
  
  old_vv <- original_scene$params$VV
  new_vv <- fit_results$final_scene$params$VV
  prey_idx <- original_scene$params$PreyFrom + 1
  pred_idx <- original_scene$params$PreyTo + 1
  sp_names <- original_scene$params$spname
  
  vv_comp <- data.frame(
    Prey     = sp_names[prey_idx],
    Predator = sp_names[pred_idx],
    Old_VV   = old_vv,
    New_VV   = new_vv,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Rel_Change = (New_VV / Old_VV) - 1,
      Log_Ratio  = log10(New_VV / Old_VV)
    )
  
  boundary_report <- vv_comp %>%
    dplyr::summarise(
      Total_Links    = dplyr::n(),
      Decreased      = sum(Rel_Change < -1e-5),
      Increased      = sum(Rel_Change > 1e-5),
      At_Floor_1.0   = sum(New_VV <= 1.0001),
      Max_Multiplier = max(New_VV / Old_VV, na.rm = TRUE),
      Min_Multiplier = min(New_VV / Old_VV, na.rm = TRUE)
    )
  message("--- Vulnerability Fit Summary ---")
  print(boundary_report)
  
  if(plot){
    plot_data <- vv_comp %>% dplyr::filter(abs(Rel_Change) > 1e-6)
    if(nrow(plot_data)>0){
      max_abs_change <- max(abs(plot_data$Rel_Change), na.rm = TRUE)
      max_log <- max(abs(vv_comp$Log_Ratio), na.rm = TRUE)
      
      plot_limit <- if(limit_type=="log") max_log else max_abs_change
      
      sub_title <- "Relative change (New/Old - 1)"
      if(!is.null(scene_name)) sub_title <- paste(sub_title, "for Scenario:",scene_name)
      
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Predator, y = Prey, fill = Rel_Change)) +
        ggplot2::geom_tile(color = "gray95", linewidth = 0.05) +
        ggplot2::scale_fill_gradientn(
          colours = c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020"),
          limits = c(-plot_limit, plot_limit),
          rescaler = scales::rescale_mid,
          oob= scales::squish,
          labels = scales::label_number(accuracy=0.1)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
                       axis.text.y = element_text(size = 7)) +
        ggplot2::labs(
          title = "Vulnerability Changes",
          subtitle = sub_title,
          fill = "Rel Change"
        )
      print(p)
      
    }
  }
  return(vv_comp)
}



#' Extract Diet and Flow Matrices for a predator
#' This function runs a simulation and extracts the flow and diet proportion matrices 
#' for a specified predator over the period interval (years)
#' 
#' @param scene an rsim scenario object
#' @param years a numeric vector of years to run the simulation
#' @param predator a character string of the predator name to extract the matrices 
#' 
#' @return a list containing two matrices: "flow" (annual flow to the predator) and "diet" (normalized diet proportions)
#' 
#' @export
#' 
#' @examples 
#' \dontrun{get_predator_diet(scene, years = 1990:2020, predator = "arrowtooth_flounder_adult")}
#' 

get_predator_diet <- function(scene, years, predator) {
  drun <- rsim.run(scene, "AB", years)
  from_species <- scene$params$spname[scene$params$PreyFrom +1]
  to_species   <- scene$params$spname[scene$params$PreyTo   +1]
  
  flow_mat <-     drun$annual_Qlink[,which(to_species == predator)]
  colnames(flow_mat) <- from_species[which(to_species == predator)]
  
  diet_mat <- flow_mat / rowSums(flow_mat)
  
  return(list(flow =flow_mat, diet= diet_mat))
}


#' Extract derivatives (derivT) for target species over time
#' 
#' Steps throught the simulation years to extract the derivative values for diagnostics
#' (gains, losses, fishing rates, etc.) for a specified target species.
#' 
#' @param scene an rsim scenario object
#' @param years a numeric vector of years to run the simulation
#' @param target_species a character string of the target species name to extract derivatives
#' 
#' @return a named list of data frames, each containing the derivatives for the target species over time
#' 
#' @export
#' 

get_species_derivatives <- function(scene, years, target_species) {
  output_list <- list()
  for(sp in target_species){
    output_list[[sp]] <- data.frame(matrix(NA, length(years), 13))
    row.names(output_list[[sp]]) <- years
  }
for(i in seq_along(years)){
  current_years <- years[1:i]
  step_run <- rsim.run(scene, method="AB", years= current_years)
  
  deriv.table <- cbind(
    Biomass=step_run$end_state$Biomass,
    as.data.frame(step_run$dyt[c("TotGain","TotLoss","DerivT", "FoodLoss",
                                           "FoodGain", "UnAssimLoss","ActiveRespLoss",
                                           "DetritalGain","FishingGain", "MzeroLoss",
                                           "FishingLoss", "DetritalLoss")])
  )
    row.names(deriv.table) <-  step_run$params$spname
    
  for(sp in target_species){
    output_list[[sp]][as.character(years[i]), ] <- deriv.table[sp,]
    
  }
}
for(sp in target_species){
  colnames(output_list[[sp]]) <- colnames(deriv.table)
  output_list[[sp]][, 2:13] <- output_list[[sp]][,2:13]/output_list[[sp]]$Biomass
}
return(output_list)
}

#' species derivates plots 
#' 
#' Plots the derivatives (biomass, fishing loss, net feeding rate) for a list of target species over time.
#' 
#' @param output_list a list of data frames generated by `get_species_derivatives()`
#' @param target_species a character vector of species to plot
#' @param scenario_name optional character string to include in plot titles for context
#' 
#' @return NULL. plots will be rendered in the plot window
#' @importFrom graphics par plot
#' @export
#' 

plot_sp_derivates <- function(output_list, target_species, scenario_name = NULL) {
  oldpar <- graphics::par(mfrow = c(length(target_species),3), mar=c(4,4,3,1))
  on.exit(graphics::par(oldpar))
  
  for(sp in target_species){
    out <- output_list[[sp]]
    
    graphics::plot(rownames(out),log(out$Biomass),xlab="year", ylab="log Biomass", 
                   main=paste(sp, scenario_name, sep=" "))
    
    graphics::plot(rownames(out), out$FishingLoss, xlab="year", ylab="F rate", 
                   main=paste(sp, scenario_name, sep=" "))
    
    graphics::plot(rownames(out), out$FoodGain - out$UnAssimLoss - out$ActiveRespLoss, xlab="year", 
         ylab="Net feeding rate", main=paste(sp, scenario_name, sep=" "))
    
    graphics::plot(rownames(out), out$ActiveRespLoss, xlab="year", ylab="F rate",
                   ylab ="resp loss", main = paste(sp, scenario_name, sep=" "))
    
  }
}

#' changes in derivatives heatmap plot 
#' 
#' Percent change to calculate the derivatives. results are shown in heatmap form
#' 
#' @param species_data a data.frame of the derivatives for a single species 
#' @param species_name a character string of the species names for the plot titles
#' @param color_cap a numeric value to cap the percent change to +-100% we will use color
#' scale squish with the default of 1 (100%)
#' 
#' @return A ggplot object.
#' @importFrom dplyr mutate across lag arrange filter %>% 
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal theme element_text labs
#' @importFrom scales percent squish
#' @export
#' 

plot_deriv_change_heatmap <- function(species_data, species_name){
  df <- species_data
  df$Year <- as.numeric(rownames(df))
  
  #year-to-year rel change
  long_df <- df %>% 
    #dplyr::arrange(Year) %>% 
    #dplyr::mutate(dplyr::across(-Year, ~(.x -dplyr::lag(.x))/abs(dplyr::lag(.x)))) %>% 
    #tidyr::pivot_longer(cols=-Year, names_to="Derivative", values_to = "Year_change") %>% 
    #dplyr::filter(!is.na(Year_change))
    
    dplyr::arrange(Year) %>% 
    dplyr::mutate(dplyr::across(-Year, ~{
      Year_change <- .x-dplyr::lag(.x)
      max_swing <- max(abs(Year_change), na.rm=TRUE)
                       if(max_swing ==0){ 
                         return(0)
                       } else {
                           return(Year_change/max_swing)
                         }
    })) %>% 
      tidyr::pivot_longer(cols=-Year, names_to="Derivative", values_to = "Scaled_change") %>% 
      dplyr::filter(!is.na(Scaled_change))
    
  y_levels <- rev(sort(unique(long_df$Year)))
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x=Derivative,y=factor(Year, levels=y_levels), fill = Scaled_change))+
    ggplot2::geom_tile(color = "gray95", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low= "#ca0020",
      mid= "white",
      high = "#0571b0",
      midpoint = 0,
      limits=c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      labels = c(" -100", "-50", "0", "50", "100")#oob= scales::squish,
      #labels=scales::percent
    )+
    ggplot2::theme_minimal()+
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle=45, hjust=1,size=7),
      axis.text.y = ggplot2::element_text(size=5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank())+
    ggplot2::labs(
      title = paste("Year to year change:", species_name),
      x= " ",
      y= " ",
      fill= "% change"
      
      
    )
  print(p)
}


