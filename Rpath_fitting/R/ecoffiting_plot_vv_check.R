# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 26 January 2026
#
# Vulnerability chack function
# 
# This function compares the vulnerability (VV) values between the original scene and the fitted results.
# usage:
# my_vv_table <- get_vv_comparison(scene_bioen, fit_results_63par$Bioen, scene_name = "63par$Bioen")
# View(my_vv_table)
# my_vv_table %>% 
#  mutate(Multiplier = New_VV / Old_VV) %>% 
#  arrange(desc(Multiplier)) %>% 
#  head(10)
#
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #



get_vv_comparison <- function(original_scene, fit_results, plot = TRUE, scene_name=NULL) {
  library(ggplot2)
  library(dplyr)
  library(scales)

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
    mutate(
      Rel_Change = (New_VV / Old_VV) - 1,
      Log_Ratio  = log10(New_VV / Old_VV)
    )

# optim bounds analysis. Max multiplier ceiling is exp(4) ~ 54.6x increase, and floor is exp(-4) ~ 0.018x decrease.
  # if max multiplier is close to 54.6, then some parameters are hitting the upper bound. If min multiplier is close to 0.018, then some parameters are hitting the lower bound. If many parameters are at the bounds, it may indicate that the optimization is struggling to find a good fit within the specified parameter space, and you may want to consider adjusting the bounds or investigating potential issues with the model or data.
  # Based on your optim bounds: lower=-4, upper=4 (log space)

  boundary_report <- vv_comp %>%
    summarise(
      Total_Links    = n(),
      Decreased      = sum(Rel_Change < -1e-5),
      Increased      = sum(Rel_Change > 1e-5),
      At_Floor_1.0   = sum(New_VV <= 1.0001),
      Max_Multiplier = max(New_VV / Old_VV, na.rm = TRUE),
      Min_Multiplier = min(New_VV / Old_VV, na.rm = TRUE)
    )
  
  message("--- Vulnerability Fit Summary ---")
  print(boundary_report)
  
  if (plot) {
    plot_data <- vv_comp %>% filter(abs(Rel_Change) > 1e-6)
    
    if(nrow(plot_data) > 0) {
      max_abs_change <- max(abs(plot_data$Rel_Change), na.rm = TRUE)
      max_log <- max(abs(vv_comp$Log_Ratio), na.rm = TRUE)
      
      sub_title <- "Relative change (New/Old - 1)"
      if(!is.null(scene_name)) sub_title <- paste(sub_title, "for Scenario:", scene_name)
      
      p <- ggplot(plot_data, aes(x = Predator, y = Prey, fill = Rel_Change)) +
        geom_tile(color = "gray95", linewidth = 0.05) +
        scale_fill_gradientn(
          colours = c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020"),
          limits = c(-max_log, max_log),
          rescaler = ~ scales::rescale_mid(.x, mid = 0),
          oob= scales::squish,
          labels = scales::label_number(accuracy=0.1)
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              axis.text.y = element_text(size = 7)) +
        labs(
          title = "Vulnerability Changes",
          subtitle = sub_title,
          fill = "Rel Change"
        )
      print(p)
    }
  }
  
  return(vv_comp)
}


