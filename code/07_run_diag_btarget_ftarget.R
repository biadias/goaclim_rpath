# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 28 July 2026
#
# code/07_run_diag_btarget_ftarget.R
# Purpose: Diagnostics of reference points (B0, Btarget, Blim).
# Compare Ftarget results from simultaneous (optim) vs
# sequential (optimize) approaches, generating diagnostic plots.
# ---------------------------------------------------------------------------- #

library(tidyverse)

# ---------------------------------------------------------------------------- #
# 1. Load data ####
# ---------------------------------------------------------------------------- #

bftarget_dir <- "data/bftarget"
plot_dir <- "data/bftarget/diagnostics"
if (!dir.exists(plot_dir))
  dir.create(plot_dir, recursive = TRUE)

# Load the simultaneous optimization results (optim)
fopt_file <- file.path(bftarget_dir, "Fopt_target_WGOA_persist.csv")
if (!file.exists(fopt_file))
  stop("Simultaneous optim results not found.")
fopt_data <- read.csv(fopt_file, row.names = 1) %>%
  rename(Fopt_Simultaneous = persist) %>%
  tibble::rownames_to_column(var = "Species")

# Load the sequential optimization results (optimize)
ftier3_file <- file.path(bftarget_dir, "F_Tier3_WGOA_persist_v2.csv")
if (!file.exists(ftier3_file))
  stop("Sequential optimize results not found.")
ftier3_data <- read.csv(ftier3_file, row.names = 1) %>%
  tibble::rownames_to_column(var = "Species")

# Load the Fmsy approach
ftier3_grid_file <- file.path(bftarget_dir, "F_Tier3_WGOA_persist_GridSearch.csv")
if (!file.exists(ftier3_grid_file))
  stop("Sequential optimize results not found.")
ftier3_grid_data <- read.csv(ftier3_grid_file, row.names = 1) %>%
  tibble::rownames_to_column(var = "Species")

# ---------------------------------------------------------------------------- #
# 2. Combining the outputs ####
# ---------------------------------------------------------------------------- #
# We compare Fopt_Simultaneous to F40_ABC since both targeted B40/B40_SSB.

compare_df <- fopt_data |> 
  left_join( ftier3_data,  by = "Species") |> 
  left_join(ftier3_grid_data, by = "Species") |> 
  mutate(
    # Absolute difference: Positive means Simultaneous estimated a higher F
    Abs_diff_fseq = Fopt_Simultaneous - F40_ABC,
    Abs_diff_fgrid = Fopt_Simultaneous - F40_ABC_Grid,
    
    # Relative difference (%)
    # Using a small offset to prevent division by zero if F40_ABC is exactly 0
    Rel_Diff_Pct_seq = (Abs_diff_fseq / (F40_ABC + 1e-9)) * 100,
    Rel_Diff_Pct_grid = (Abs_diff_fgrid / (F40_ABC_Grid + 1e-9)) * 100
  )


write.csv(
  compare_df,
  file.path(plot_dir, "F_optimization_comparison_table.csv"),
  row.names = FALSE
)
message("Comparison table saved to: ", plot_dir)

cat("\n=== Optimization Comparison Summary ===\n")
print(summary(compare_df[, c("Fopt_Simultaneous", "F40_ABC","F40_ABC_Grid", "Abs_diff_fseq")]))

# ---------------------------------------------------------------------------- #
# 3. Diagnostic Plots ####
# ---------------------------------------------------------------------------- #

## Plot 1: 1-to-1 Scatter Plot
# Shows how closely the two methods align overall.
# Deviations from the dashed line indicate ecological trade-offs captured by the simultaneous method.
p1a <- ggplot(compare_df, aes(x = F40_ABC, y = F40_ABC_Grid)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "darkred",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_point(size = 3,
             alpha = 0.7,
             color = "#0072B2") +
  ggrepel::geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15) +
  labs(
    title = "Comparison of Target F Methodologies",
    subtitle = "Sequential (optimize) vs FMSY grid",
    x = "Sequential 1D F40 (optimize)",
    y = "FMSY_F40 (grid)"
  ) +
  theme_minimal() +
  coord_fixed(
    ratio = 1,
    xlim = c(0, 1.5),
    #max(compare_df$F40_ABC, compare_df$Fopt_Simultaneous)),
    ylim = c(0, 1.5)
  ) #max(compare_df$F40_ABC, compare_df$Fopt_Simultaneous)))

ggsave(
  file.path(plot_dir, "Diag1_1to1_Scatter_1a.png"),
  plot = p1a,
  width = 8,
  height = 8,
  bg = "white"
)

p1b <- ggplot(compare_df, aes(x = F40_ABC, y = Fopt_Simultaneous)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "darkred",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_point(size = 3,
             alpha = 0.7,
             color = "#0072B2") +
  ggrepel::geom_text_repel(aes(label = Species), size = 3, max.overlaps = 15) +
  labs(
    title = "Comparison of Target F Methodologies",
    subtitle = "Sequential (optimize) vs Simultaneous (optim)",
    x = "Sequential 1D F40 (optimize)",
    y = "simultaneous F40 (optim)"
  ) +
  theme_minimal() +
  coord_fixed(
    ratio = 1,
    xlim = c(0, 1.5),
    #max(compare_df$F40_ABC, compare_df$Fopt_Simultaneous)),
    ylim = c(0, 1.5)
  ) #max(compare_df$F40_ABC, compare_df$Fopt_Simultaneous)))

ggsave(
  file.path(plot_dir, "Diag1_1to1_Scatter_1b.png"),
  plot = p1b,
  width = 8,
  height = 8,
  bg = "white"
)


## Plot 2: Absolute Difference Bar Chart (I don't really like this one!)
# Highlights which specific species are most sensitive to multi-species interactions
p2 <- ggplot(compare_df, aes(
  x = reorder(Species, Abs_Diff),
  y = Abs_Diff,
  fill = Abs_Diff > 0
)) +
  geom_col(color = "black", alpha = 0.8) +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#0072B2", "FALSE" = "#E69F00"),
    labels = c("TRUE" = "Simultaneous > Sequential", "FALSE" = "Sequential > Simultaneous")
  ) +
  labs(title = "Fopt (Simultaneous) minus F40 (Sequential)",
       x = "Species",
       y = "Absolute Difference in F",
       fill = "Direction") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(
  file.path(plot_dir, "Diag2_Abs_Difference_Bar.png"),
  plot = p2,
  width = 8,
  height = 6,
  bg = "white"
)


## Plot 3: Side-by-Side Dumbbell/Lollipop Plot (good one!)
# Good for seeing the absolute magnitude of F for both methods side-by-side
compare_long <- compare_df %>%
  select(Species, Fopt_Simultaneous, F40_ABC, F40_ABC_Grid) %>%
  pivot_longer(
    cols = c("Fopt_Simultaneous", "F40_ABC", "F40_ABC_Grid"),
    names_to = "Method",
    values_to = "F_Estimate"
  )

p3 <- ggplot(compare_long, aes(
  x = F_Estimate,
  y = reorder(Species, F_Estimate),
  color = Method
)) +
  geom_line(aes(group = Species),
            color = "gray70",
            linewidth = 1) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(
    values = c(
      "F40_ABC" = "#E69F00",
      "Fopt_Simultaneous" = "#0072B2",
      "F40_ABC_Grid" = "#D95F02"
    ),
    labels = c("F40_ABC" = "Sequential (1D)", "Fopt_Simultaneous" = "Simultaneous (optim)", "F40_ABC_Grid" = "Grid (1D)")
  ) +
  labs(
    title = "Target F Estimates by Species",
    x = "Fishing Mortality (F)",
    y = "Species",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(
  file.path(plot_dir, "Diag3_Dumbbell_Plot.png"),
  plot = p3,
  width = 8,
  height = 8,
  bg = "white"
)

message("Diagnostic plots successfully generated in: ", plot_dir)

# ---------------------------------------------------------------------------- #
# 4. Diag for simulated and target ####
# ---------------------------------------------------------------------------- #

## 4.1. Define the species and the capped F values ####
# print(round(Fcomp_final, 4)) #check the test F using the F40_Opt column for the FG you wanna test
# note to self: please make a for loop for this

sp <- "pacific_ocean_perch_adult"
test_F <- 0.0370 # with penalty for stock depletion look at the sumsq_btarg_single() from 06_run script
## 4.2. Apply the test F to the projection rows ####
scene_persist_fopt$fishing$ForcedFRate[proj_rows, sp] <- test_F

## 4.3. Run the simulation ####
run_test <- rsim.run(scene_persist_fopt, method = "AB", years = all_years)

## 4.4. Extract the simulated end-of-century biomass (checking SSB if it's an SSB stock) ####
is_ssb <- sp %in% ssb_stocks
simulated_bio <- if (is_ssb) {
  end_cent_SSB(run_test)[sp, ]
} else {
  end_cent_biomass(run_test)[sp, ]
}

## 4.5. Extract the B40 target for comparison ####
target_bio <- if (is_ssb) {
  Btarg_all[["persist"]][sp, "B40_SSB"]
} else {
  Btarg_all[["persist"]][sp, "B40"]
}

## 4.6. Print the results side-by-side ####
message("Diagnostic Check for: ", sp)
message("Simulated Biomass at F=", test_F, " : ", round(simulated_bio, 2))
message("Target B40                 : ", round(target_bio, 2))
message("Ratio (Simulated / Target) : ",
        round(simulated_bio / target_bio, 3))
#If the Ratio is exactly 1.0: The optimizer perfectly hit the target right at the boundary limit
#If the Ratio is > 1.0 (e.g., 1.5): This means the simulated biomass is still higher than B40
#If the Ratio is < 1.0: The optimizer overshot the target


# Reset POP back to status quo just to be safe
scene_persist_fopt$fishing$ForcedFRate[proj_rows, sp] <- F_meanlast[sp]

# ---------------------------------------------------------------------------- #
# 5. Test for sablefish specific hump on the SSB ####
# ---------------------------------------------------------------------------- #

sp <- "sablefish_adult"

## 5.1 Build the baseline scene (same as scene_persist_f40 in Step 9 - 06_run_btarget script) ####
scene_sab <- F_clim_sim_scene(
  scene=scene_bioen_best, ssp="persist", cons=TRUE, resp=TRUE, buf=FALSE,
  pcod_rec=TRUE, pcod_rec_method="cauchy", bioen_sp=bioen_sp,
  tdc_hind=tdc_hind, tdr_hind=tdr_hind, managed_sp=managed_sp_list,
  f_equil=F_equil, f_zero=F_zero, f_scenario="mean", f_ref_yrs=2016:2020,
  climate_dir="data/climate/", hind_yrs=hind_years, proj_yrs=2021:2099,
  hind_data_start_yr=1991, climate_data_start_yr=1980, verbose=FALSE
)
frate_cols   <- colnames(scene_sab$fishing$ForcedFRate)
equil_cols   <- intersect(frate_cols, names(F_equil))
managed_cols <- intersect(frate_cols, managed_sp_list)
scene_sab$fishing$ForcedFRate[proj_rows, equil_cols]   <-
  matrix(F_equil[equil_cols],   nrow=length(proj_rows), ncol=length(equil_cols),   byrow=TRUE)
scene_sab$fishing$ForcedFRate[proj_rows, managed_cols] <-
  matrix(F_meanlast[managed_cols], nrow=length(proj_rows), ncol=length(managed_cols), byrow=TRUE)

## 5.2 Dense F scan: 0 to 1.0 ####
f_scan <- c(seq(0, 0.05, by=0.005), seq(0.05, 0.5, by=0.02), seq(0.5, 1.0, by=0.05))
f_scan <- sort(unique(f_scan))

ssb_vals <- numeric(length(f_scan))
bio_vals <- numeric(length(f_scan))

# Check what juvenile stanza F is set to in the projection
juve_col <- grep("sablefish_juvenile", frate_cols, value=TRUE)
message("Sablefish juvenile in ForcedFRate?:", length(juve_col) > 0, "\n")

if (length(juve_col) > 0)
  message("Juvenile F in proj (first value):", scene_sab$fishing$ForcedFRate[proj_rows[1], juve_col], "\n")


ptm <- proc.time()
for (idx in seq_along(f_scan)) {
  scene_sab$fishing$ForcedFRate[proj_rows, sp] <- f_scan[idx]
  run_tmp <- rsim.run(scene_sab, method="AB", years=all_years)
  ssb_vals[idx] <- end_cent_SSB(run_tmp)[sp, ]
  bio_vals[idx] <- end_cent_biomass(run_tmp)[sp, ]
}
message("Scan elapsed:", round((proc.time()-ptm)[3]/60, 2), "min\n")

sab_scan <- data.frame(f=f_scan, SSB=ssb_vals, Biomass=bio_vals)
message("\nB40_SSB target:", round(Btarg_all[["persist"]][sp, "B40_SSB"], 4), "\n")
message("Min SSB in scan:", round(min(ssb_vals), 4), "at F =", f_scan[which.min(ssb_vals)], "\n")
message("Max SSB in scan:", round(max(ssb_vals), 4), "at F =", f_scan[which.max(ssb_vals)], "\n")
message("SSB at F=0     :", round(ssb_vals[1], 4), "\n")
message("Biomass at F=0 :", round(bio_vals[1], 4), "\n")

library(ggplot2)

b40_ssb <- Btarg_all[["persist"]][sp, "B40_SSB"]
b35_ssb <- Btarg_all[["persist"]][sp, "B35_SSB"]
b40_bio <- Btarg_all[["persist"]][sp, "B40"]
b35_bio <- Btarg_all[["persist"]][sp, "B35"]

## Where does SSB cross B40_SSB? ####
above <- sab_scan$SSB > b40_ssb
cross_idx <- which(diff(above) != 0)  # indices where sign flips
message("B40_SSB crossings at F ≈:", round(sab_scan$f[cross_idx + 1], 4), "\n")

cross_bio <- sab_scan$Biomass > b40_bio
message("B40_bio crossings at F ≈:", round(sab_scan$f[which(diff(cross_bio) != 0) + 1], 4), "\n")

# Plot both curves with reference lines
df_long <- data.frame(
  f       = rep(sab_scan$f, 2),
  value   = c(sab_scan$SSB, sab_scan$Biomass),
  metric  = rep(c("SSB", "Total Biomass"), each=nrow(sab_scan))
)

ggplot(df_long, aes(x=f, y=value)) +
  geom_line(color="#0072B2", linewidth=1) +
  geom_point(color="#0072B2", size=1.5) +
  geom_hline(aes(yintercept=ifelse(metric=="SSB", b40_ssb, b40_bio),
                 linetype="B40"), color="#009E73", linewidth=0.8) +
  geom_hline(aes(yintercept=ifelse(metric=="SSB", b35_ssb, b35_bio),
                 linetype="B35"), color="#D95F02", linewidth=0.8) +
  geom_vline(xintercept=F_meanlast[sp], linetype="dotted", color="gray40") +
  scale_linetype_manual(name="Target", values=c("B40"="dashed","B35"="dotdash")) +
  facet_wrap(~metric, scales="free_y") +
  labs(title="Sablefish adult: end-of-century metric vs. F",
       subtitle=paste0("Dotted vertical = F_meanlast (", round(F_meanlast[sp],4), ")"),
       x="Fishing Mortality Rate (F)", y="End-of-century value (t/km²)") +
  theme_bw()

