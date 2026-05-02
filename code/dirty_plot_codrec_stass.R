library(tidyverse)
library(tidyverse)

# Subset for 1991-2020 (hindcast)
data_comp <- tibble(
  year = 1991:2020,
  rec_191b = c(0.45, 0.41, 0.29, 0.33, 0.44, 0.29, 0.30, 0.24, 0.33, 0.38, 
               0.27, 0.18, 0.22, 0.26, 0.39, 0.58, 0.45, 0.57, 0.43, 0.42, 
               0.54, 1.05, 0.69, 0.27, 0.26, 0.26, 0.20, 0.16, 0.09, 0.15),
  
  #scene_full_pcod$forcing$ForcedRecs[,"pacific_cod_adult"]
  pcod_rec = c(0.838451, 0.985377, 0.944768, 0.952599, 0.905512, 0.986808, 
               0.984284, 0.931014, 0.944777, 0.983264, 0.946236, 0.992783, 
               0.824655, 0.999290, 0.965411, 0.998714, 0.882177, 0.942311, 
               0.929522, 0.979133, 0.973260, 0.860452, 1.000017, 0.935743, 
               0.763282, 0.698916, 0.973151, 0.924720, 0.723440, 0.999535)
)

# center model 19.1b around 1
# Note: The mean used here (0.42) is based on the full (1977-2023)-2 time series
full_series_mean <- 0.42

data_comp <- data_comp %>%
  mutate(
    rec_191b_centered = rec_191b / full_series_mean
  )

# 3. Visualization
ggplot(data_comp, aes(x = year)) +
  # Model 19.1b Centered Index
  geom_line(aes(y = rec_191b_centered, color = "Model 19.1b (Centered)"), linewidth = 1.2,
            linetype = 4) +
  geom_point(aes(y = rec_191b_centered, color = "Model 19.1b (Centered)")) +
  # Pcod Recruitment Index
  geom_line(aes(y = pcod_rec, color = "Pcod Recruitment"), linewidth = 1.2) +
  geom_point(aes(y = pcod_rec, color = "Pcod Recruitment")) +
  # Baseline average
  geom_hline(yintercept = 1, color = "grey", linetype = "dotted", linewidth= 1.5) +
  scale_color_manual(values = c("Model 19.1b (Centered)" = "#2c7fb8", 
                                "Pcod Recruitment" = "#31a354")) +
  labs(
    title = "Comparison of Centered Recruitment Indices",
    subtitle = "Both indices centered around 1.0 (Long-term Average)",
    y = "Relative Index (Value / Mean)",
    x = "Year",
    color = " "
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
