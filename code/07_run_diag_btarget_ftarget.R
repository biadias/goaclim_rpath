# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 08 June 2026
#
# code/07_run_diag_btarget_ftarget.R
# Purpose: Diagnostics of reference points (B0, Btarget, Blim).
# Compare Ftarget results from simultaneous (optim) vs 
# sequential (optimize) approaches, generating diagnostic plots.
# ---------------------------------------------------------------------------- #

library(dplyr)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------------------------- #
# 1. Configuration & Data Loading
# ---------------------------------------------------------------------------- #

bftarget_dir <- "data/bftarget"
plot_dir <- "data/bftarget/diagnostics"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Load the simultaneous optimization results (optim)
fopt_file <- file.path(bftarget_dir, "Fopt_target_WGOA_GFDL_persist.csv")
if (!file.exists(fopt_file)) stop("Simultaneous optim results not found.")
fopt_data <- read.csv(fopt_file, row.names = 1) %>%
  rename(Fopt_Simultaneous = persist) %>%
  tibble::rownames_to_column(var = "Species")

# Load the sequential optimization results (optimize)
ftier3_file <- file.path(bftarget_dir, "F_Tier3_WGOA_GFDL_persist.csv")
if (!file.exists(ftier3_file)) stop("Sequential optimize results not found.")
ftier3_data <- read.csv(ftier3_file, row.names = 1) %>%
  tibble::rownames_to_column(var = "Species")

# ---------------------------------------------------------------------------- #
# 2. Data Merging & Deviation Calculations
# ---------------------------------------------------------------------------- #
# We compare Fopt_Simultaneous to F40_ABC since both targeted B40/B40_SSB.

compare_df <- left_join(fopt_data, ftier3_data, by = "Species") %>%
  mutate(
    # Absolute difference: Positive means Simultaneous estimated a higher F
    Abs_Diff = Fopt_Simultaneous - F40_ABC,
    
    # Relative difference (percentage)
    # Using a small offset to prevent division by zero if F40_ABC is exactly 0
    Rel_Diff_Pct = (Abs_Diff / (F40_ABC + 1e-9)) * 100
  )

# Save the raw comparison table
write.csv(compare_df, file.path(plot_dir, "F_optimization_comparison_table.csv"), row.names = FALSE)
message("Comparison table saved to: ", plot_dir)

# Print a quick summary to console
cat("\n=== Optimization Comparison Summary ===\n")
print(summary(compare_df[, c("Fopt_Simultaneous", "F40_ABC", "Abs_Diff")]))

# ---------------------------------------------------------------------------- #
# 3. Diagnostic Plots
# ---------------------------------------------------------------------------- #

## Plot 1: 1-to-1 Scatter Plot 
# Shows how closely the two methods align overall. Deviations from the dashed line 
# indicate ecological trade-offs captured by the simultaneous method.
p1 <- ggplot(compare_df, aes(x = F40_ABC, y = Fopt_Simultaneous)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(size = 3, alpha = 0.7, color = "dodgerblue4") +
  ggrepel::geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15) +
  labs(
    title = "Comparison of Target F Methodologies",
    subtitle = "Simultaneous (optim) vs Sequential (optimize)",
    x = "Sequential 1D F40 (optimize)",
    y = "Simultaneous Multispecies F40 (optim)"
  ) +
  theme_minimal() +
  coord_fixed(ratio = 1, xlim = c(0, max(compare_df$F40_ABC, compare_df$Fopt_Simultaneous)), 
              ylim = c(0, max(compare_df$F40_ABC, compare_df$Fopt_Simultaneous)))

ggsave(file.path(plot_dir, "Diag1_1to1_Scatter.png"), plot = p1, width = 8, height = 8, bg = "white")


## Plot 2: Absolute Difference Bar Chart
# Highlights which specific species are most sensitive to multi-species interactions
p2 <- ggplot(compare_df, aes(x = reorder(Species, Abs_Diff), y = Abs_Diff, fill = Abs_Diff > 0)) +
  geom_col(color = "black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "darkorange"), 
                    labels = c("TRUE" = "Simultaneous > Sequential", 
                               "FALSE" = "Sequential > Simultaneous")) +
  labs(
    title = "Absolute Difference in F Estimates",
    subtitle = "Fopt (Simultaneous) minus F40 (Sequential)",
    x = "Species",
    y = "Absolute Difference in F",
    fill = "Direction"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "Diag2_Abs_Difference_Bar.png"), plot = p2, width = 8, height = 6, bg = "white")


## Plot 3: Side-by-Side Dumbbell/Lollipop Plot
# Good for seeing the absolute magnitude of F for both methods side-by-side
compare_long <- compare_df %>%
  select(Species, Fopt_Simultaneous, F40_ABC) %>%
  pivot_longer(cols = c("Fopt_Simultaneous", "F40_ABC"), 
               names_to = "Method", values_to = "F_Estimate")

p3 <- ggplot(compare_long, aes(x = F_Estimate, y = reorder(Species, F_Estimate), color = Method)) +
  geom_line(aes(group = Species), color = "gray70", linewidth = 1) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("F40_ABC" = "darkorange", "Fopt_Simultaneous" = "steelblue"),
                     labels = c("F40_ABC" = "Sequential (1D)", "Fopt_Simultaneous" = "Simultaneous (N-D)")) +
  labs(
    title = "Target F Estimates by Species",
    subtitle = "Connecting lines show the gap between methodologies",
    x = "Fishing Mortality (F)",
    y = "Species",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(file.path(plot_dir, "Diag3_Dumbbell_Plot.png"), plot = p3, width = 8, height = 8, bg = "white")

message("Diagnostic plots successfully generated in: ", plot_dir)