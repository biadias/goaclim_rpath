library(ggplot2)
library(dplyr)
library(readr)


#roms_temp_month <- read_csv("Rpath_fitting/GOA/wgoa_data_rpath_fitting/Long_WGOA_temp_monthly_1000.csv")
roms_temp_annual<- read_csv("Rpath_fitting/GOA/wgoa_data_rpath_fitting/Long_WGOA_temp_annual_1000.csv")

# ---------------------------------------------------------------------------- # 
# 1. Plot for Bottom depthclass ####
# ---------------------------------------------------------------------------- # 
# bottom temp
bottom_roms_temp_annual<- roms_temp_annual%>% filter(depthclass == "Bottom")

# Plot
bottom_plot <- ggplot(bottom_df, aes(x = year, y = annual_avg_temp , color = simulation)) +
  geom_line(alpha = 0.8, size = 1) +
  labs(
    title = "Annual Average Bottom Temperature by Simulation Scenario",
    x = "Year",
    y = "Temperature",
    color = "Simulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray", size = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

bottom_plot


# ---------------------------------------------------------------------------- #
# 2.  Surface depthclass ####
# ---------------------------------------------------------------------------- #

surface_roms_temp_annual<- roms_temp_annual%>% filter(depthclass == "Surface")

surface_plot <- ggplot(surface_df, aes(x = year, y = annual_avg_temp, color = simulation)) +
  geom_line(alpha = 0.8, size = 1) +
  labs(
    title = "Annual Average Surface Temperature by Simulation Scenario",
    x = "Year",
    y = "Temperature",
    color = "Simulation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(linetype = "dashed", color = "gray", size = 0.5),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

surface_plot