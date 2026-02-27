#------------------------------------------------------------------------------#
#AUTHORS: Bia Dias 
#AFFILIATIONS: Alaska Fisheries Science Center / CICOES University of Washington
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#
# This function creates the scenario comparison plots for biomass fitting for all the species.
#------------------------------------------------------------------------------#



library(tidyverse)

goa_results <- readRDS("GOA/GOA_fit_results_59M04par.rds")# This was my best model
goa_results$Bioen

catch_raw <- goa_results$Bioen$final_run$annual_Catch
biomass_raw <- goa_results$Bioen$final_run$annual_Biomass

long_fish_dt <- function(mat, value_name) {
  mat %>%
    as.data.frame() %>%
    rownames_to_column("Year") %>%
    pivot_longer(-Year, names_to = "Species", values_to = value_name) %>%
    mutate(Year = as.numeric(Year))
}

# F calculation
f_rates_data <- long_fish_dt(catch_raw, "Catch") %>%
  left_join(long_fish_dt(biomass_raw, "Biomass"), by = c("Year", "Species")) %>%
  mutate(F_rate = Catch / Biomass)

f_summary <- f_rates_data %>%
  mutate(
    Period = case_when(
      Year >= 2016 & Year <= 2020 ~ "F 2016-2020",
      Year >= 2012 & Year <= 2021 ~ "F 2011-2020",
      Year >= 1991 & Year <= 2021 ~ "F 1991-2020",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period)) %>%
  group_by(Species, Period) %>%
  summarize(Mean_F = mean(F_rate, na.rm = TRUE), .groups = "drop")
f_summary

#write.csv(f_summary, "GOA/f_rate_summary.csv", row.names = FALSE)

plot_data <- f_summary %>%
  filter(
    Species %in% c(
      "arrowtooth_flounder_adult",
      "pacific_cod_adult",
      "sablefish_adult",
      "walleye_pollock_adult"
    )
  )

plot_data$Period <- factor(plot_data$Period,
                           levels = c("F 1991-2020", "F 2011-2020", "F 2016-2020"))

species_order <- colnames(catch_raw)
plot_data$Species <- factor(plot_data$Species, levels = species_order)

#write.csv(plot_data, "GOA/f_rate_summary_gf.csv", row.names = FALSE)

# Plot
f_rate_plot <- ggplot(plot_data, aes(x = Species, y = Mean_F, fill = Period)) +
  geom_histogram(stat = "identity", position = position_dodge(), alpha=0.7) +
  scale_fill_manual(values = c(
    "F 1991-2020" = "#e69f00",
    "F 2011-2020" = "#56b4e9",
    "F 2016-2020" = "#009e73"
  )) +
  labs(title = "Mean F Rates by Species and Period", x = "Species", y = "Mean F Rate") +
  theme_minimal() +
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    size = 12
  ))
