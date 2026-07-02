
#Checks

run_Fmean    <- rsim.run(scene_persist_pcod_f0c,    method = "AB", years = all_years)
 run_F0    <- rsim.run(scene_persist_pcod_f0g,    method = "AB", years = all_years)
 check_sp <- c("pacific_ocean_perch_adult", "sablefish_adult", "flathead_sole_adult",
                +               "pacific_herring_adult", "large_copepods", "euphausiids")
 
   for (sp in check_sp) {
         diff_pct <- 100 * (run_F0$annual_Biomass[, sp] - run_Fmean$annual_Biomass[, sp]) / 
               run_Fmean$annual_Biomass[, sp]
         cat(sprintf("%-30s max abs diff: %.2f%%\n", sp, max(abs(diff_pct))))
     }
#pacific_ocean_perch_adult      max abs diff: 10.07%
#sablefish_adult                max abs diff: 120.09%
#flathead_sole_adult            max abs diff: 54.81%
#pacific_herring_adult          max abs diff: 63.49%
#large_copepods                 max abs diff: 8.51%
#euphausiids                    max abs diff: 1.27%
 for (sp in c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")) {
       end_F0    <- mean(tail(run_F0$annual_Biomass[, sp], 20))
       end_Fmean <- mean(tail(run_Fmean$annual_Biomass[, sp], 20))
       cat(sprintf("%-30s F0/Fmean end-of-century: %.2f\n", sp, end_F0/end_Fmean))
   }
#walleye_pollock_adult          F0/Fmean end-of-century: 3.22
#pacific_cod_adult              F0/Fmean end-of-century: 3.41
#arrowtooth_flounder_adult      F0/Fmean end-of-century: 0.39
 for (sp in c("walleye_pollock_juvenile", "pacific_cod_juvenile", "arrowtooth_flounder_juvenile")) {
       end_F0    <- mean(tail(run_F0$annual_Biomass[, sp], 20))
       end_Fmean <- mean(tail(run_Fmean$annual_Biomass[, sp], 20))
       cat(sprintf("%-30s F0/Fmean end-of-century: %.2f\n", sp, end_F0/end_Fmean))}

 #PCOD testing####

source("code/Function_F_clim_sim_scene_v2.R")

 # Make sure the patched pcod_rec block is in F_clim_sim_scene.R first.
 # Then source the function and run:
 
 scene_persist_F0_norec <- F_clim_sim_scene(
   scene = scene_bioen_best,
   ssp = "persist",
   cons = TRUE, resp = TRUE, buf = FALSE,
   pcod_rec = FALSE,
   bioen_sp = bioen_sp,
   tdc_hind = tdc_hind, tdr_hind = tdr_hind,
   managed_sp = managed_sp_list,
   f_equil = F_equil, f_zero = F_zero,
   f_scenario = "zero",
   f_ref_yrs = 2016:2020,
   climate_dir = "data/climate/",
   hind_yrs = hind_years,
   proj_yrs = 2021:2099,
   hind_data_start_yr = 1991,
   climate_data_start_yr = 1980,
   verbose = TRUE)
 
 scene_persist_F0_rec <- F_clim_sim_scene(
   scene = scene_bioen_best,
   ssp = "persist",
   cons = TRUE, resp = TRUE, buf = FALSE,
   pcod_rec = TRUE,
   pcod_rec_method = "cauchy",
   bioen_sp = bioen_sp,
   tdc_hind = tdc_hind, tdr_hind = tdr_hind,
   managed_sp = managed_sp_list,
   f_equil = F_equil, f_zero = F_zero,
   f_scenario = "zero",
   f_ref_yrs = 2016:2020,
   climate_dir = "data/climate/",
   hind_yrs = hind_years,
   proj_yrs = 2021:2099,
   hind_data_start_yr = 1991,
   climate_data_start_yr = 1980,
   verbose = TRUE)
 
 run_norec <- rsim.run(scene_persist_F0_norec, method = "AB", years = all_years)
 run_rec   <- rsim.run(scene_persist_F0_rec,   method = "AB", years = all_years)

 scene_126_F0_norec <- F_clim_sim_scene(
   scene = scene_bioen_best,
   ssp = "126",
   cons = TRUE, resp = TRUE, buf = FALSE,
   pcod_rec = FALSE,
   bioen_sp = bioen_sp,
   tdc_hind = tdc_hind, tdr_hind = tdr_hind,
   managed_sp = managed_sp_list,
   f_equil = F_equil, f_zero = F_zero,
   f_scenario = "zero",
   f_ref_yrs = 2016:2020,
   climate_dir = "data/climate/",
   hind_yrs = hind_years,
   proj_yrs = 2021:2099,
   hind_data_start_yr = 1991,
   climate_data_start_yr = 1980,
   verbose = TRUE)
 
 scene_126_F0_rec <- F_clim_sim_scene(
   scene = scene_bioen_best,
   ssp = "126",
   cons = TRUE, resp = TRUE, buf = FALSE,
   pcod_rec = TRUE,
   pcod_rec_method = "cauchy",
   bioen_sp = bioen_sp,
   tdc_hind = tdc_hind, tdr_hind = tdr_hind,
   managed_sp = managed_sp_list,
   f_equil = F_equil, f_zero = F_zero,
   f_scenario = "zero",
   f_ref_yrs = 2016:2020,
   climate_dir = "data/climate/",
   hind_yrs = hind_years,
   proj_yrs = 2021:2099,
   hind_data_start_yr = 1991,
   climate_data_start_yr = 1980,
   verbose = TRUE)
 
 run_126norec <- rsim.run(scene_126_F0_norec, method = "AB", years = all_years)
 run_126rec   <- rsim.run(scene_126_F0_rec,   method = "AB", years = all_years)

 # --- Check 1: verbose output from the function call -----------------------
 # When you run scene_persist_F0_rec, the verbose=TRUE will print:
 #   "Hindcast Feb-Apr btemp climatology: X.XXX C (proxy = Y.YYYY)"
 #   "Projection multiplier range: [A.AAA, B.BBB], mean: C.CCC"
 #   "End-of-century (last 20 yrs) mean multiplier: D.DDD"
 # 
 # Expected: under persist, all three projection lines should show 1.000
 
 # --- Check 2: ForcedRecs matrix inspection --------------------------------
 fr <- scene_persist_F0_rec$forcing$ForcedRecs[, "pacific_cod_adult"]
 cat(sprintf("ForcedRecs cod range: [%.4f, %.4f]\n", min(fr), max(fr)))
 cat(sprintf("Hindcast (first %d months) mean: %.4f\n",
             length(hind_years)*12, mean(fr[1:(length(hind_years)*12)])))
 cat(sprintf("Projection (last %d months) mean: %.4f\n",
             length(2021:2099)*12, mean(fr[(length(hind_years)*12+1):length(fr)])))
 # Expected: hindcast = 1.0000, projection = 1.0000 (persist matches climatology)
 
 # --- Check 3: cod biomass with vs without recruitment forcing -------------
 end_norec <- mean(tail(run_norec$annual_Biomass[, "pacific_cod_adult"], 20))
 end_rec   <- mean(tail(run_rec$annual_Biomass[, "pacific_cod_adult"],   20))
 cat(sprintf("\nPersist cod end-of-century biomass:\n"))
 cat(sprintf("  rec_off: %.4f\n", end_norec))
 cat(sprintf("  rec_on:  %.4f\n", end_rec))
 cat(sprintf("  diff:    %.3f%%\n", 100 * (end_rec - end_norec) / end_norec))
 # Expected: <1% diff (essentially identical)
 
 # --- Check 4: hindcast must be byte-identical -----------------------------
 hind_rows <- 1:length(hind_years)
 identical(run_rec$annual_Biomass[hind_rows, ], run_norec$annual_Biomass[hind_rows, ])
 # Expected: TRUE
 
 # --- Check 5: full trajectory (optional, for visual confirmation) ---------
 plot(run_norec$annual_Biomass[, "pacific_cod_adult"], type = "l",
      xlab = "year index", ylab = "cod adult biomass",
      main = "585 Fmean: cod with vs without recruitment forcing",
      col = "black", lwd = 2)
 
 lines(run_rec$annual_Biomass[, "pacific_cod_adult"], col = "red", lwd = 2, lty = 2)
 lines(run_126rec$annual_Biomass[, "pacific_cod_adult"], col = "blue", lwd = 2, lty = 2)
 lines(run_126norec$annual_Biomass[, "pacific_cod_adult"], col = "blue", lwd = 2, lty = 1)
 legend("topleft", legend = c("rec_off_persist", "rec_on_persist (cauchy)"),
        col = c("black", "red"), lty = c(1, 2), lwd = 2, bty = "n")
 # Expected: the two lines should be visually identical

 ssps <- c("126", "245", "585", "persist")
 rec_options <- c(FALSE, TRUE)
 
 # 2. Nested loop: SSPs by Recruitment status
 for (s in ssps) {
   for (r in rec_options) {
     
     # Define a suffix for naming (e.g., "rec" or "norec")
     rec_suffix <- ifelse(r, "rec", "norec")
     
     # Create the scene object name (e.g., scene_126_F0_rec)
     scene_name <- paste0("scene_", s, "_Fmean_", rec_suffix)
     
     # Generate the scene using assign() to save it to the environment
     assign(scene_name, F_clim_sim_scene(
       scene = scene_bioen_best,
       ssp = s,
       cons = TRUE, resp = TRUE, buf = FALSE,
       pcod_rec = r,
       pcod_rec_method = if(r) "cauchy" else NULL, # Only use cauchy if rec is TRUE
       bioen_sp = bioen_sp,
       tdc_hind = tdc_hind, tdr_hind = tdr_hind,
       managed_sp = managed_sp_list,
       f_equil = F_equil, f_zero = F_zero,
       f_scenario = "zero",
       f_ref_yrs = 2016:2020,
       climate_dir = "data/climate/",
       hind_yrs = hind_years,
       proj_yrs = 2021:2099,
       hind_data_start_yr = 1991,
       climate_data_start_yr = 1980,
       verbose = FALSE # Set to FALSE to keep the console clean
     ))
     
     # Create the run object name (e.g., run_126rec)
     # Note: Using your naming convention where 'persist' doesn't have the SSP number
     run_name <- if(s == "persist") paste0("run_", rec_suffix) else paste0("run_", s, rec_suffix)
     
     # Execute the simulation and save the result
     # get(scene_name) retrieves the object we just created above
     assign(run_name, rsim.run(get(scene_name), method = "AB", years = all_years))
     
     message("Finished: ", run_name)
   }
 }
  
 # --------------------------------------------------------------------------- #
 #PLOT ####
 # --------------------------------------------------------------------------- #


 species_list <- c("arrowtooth_flounder_adult", "pacific_cod_adult", "walleye_pollock_adult",
                   "arrowtooth_flounder_juvenile", "pacific_cod_juvenile", "walleye_pollock_juvenile")
 ssps <- c("persist", "126", "245", "585")
 rec_types <- c("norec", "rec")
 
base_cb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 
 cols <- c("persist" = "black", "126" = "#009E73", "245" = "#56B4E9", "585" = "#E69F00")
 ltys <- c("norec" = 2, "rec" = 1)
 
 
 # Optional: If you want all plots to appear in a grid in the same window, uncomment below:
  par(mfrow = c(2, 3)) # Adjust dimensions (rows, columns) based on how many species you have
 
 # Outer Loop: Iterate through each species
 for (sp in species_list) {
   
   # --- NEW: Find the global maximum biomass for this species across ALL runs ---
   max_bio <- 0
   for (s in ssps) {
     for (r in rec_types) {
       obj_name <- if(s == "persist") paste0("run_", r) else paste0("run_", s, r)
       if (exists(obj_name)) {
         temp_run <- get(obj_name)
         run_max <- max(temp_run$annual_Biomass[, sp], na.rm = TRUE)
         if (run_max > max_bio) {
           max_bio <- run_max
         }
       }
     }
   }
   
   # 3. Initialize the plot
   first_name <- if(ssps[1] == "persist") paste0("run_", rec_types[1]) else paste0("run_", ssps[1], rec_types[1])
   initial_obj <- get(first_name)
   
   # Note the dynamic y-axis and dynamic titles using 'sp'
   plot(initial_obj$annual_Biomass[, sp], 
        type = "n", 
        ylim = c(0, max_bio * 1.1), # 10% buffer above the absolute highest value
        xlab = "Year Index", 
        ylab = paste(sp, "Biomass"),
        main = paste(sp, "Trajectories by SSP & Rec\nF0"))
   
   # 4. Loop through and add lines
   for (s in ssps) {
     for (r in rec_types) {
       
       obj_name <- if(s == "persist") paste0("run_", r) else paste0("run_", s, r)
       
       if (exists(obj_name)) {
         current_run <- get(obj_name)
         lines(current_run$annual_Biomass[, sp], 
               col = cols[s], 
               lty = ltys[r], 
               lwd = 2)
       }
     }
   }
   
   # 5. Create a Dynamic Legend
   legend_labels <- expand.grid(rec_types, ssps)
   legend_text <- paste(legend_labels$Var2, legend_labels$Var1)
   legend_cols <- cols[as.character(legend_labels$Var2)]
   legend_ltys <- ltys[as.character(legend_labels$Var1)]
   
   legend("bottomleft", 
          legend = legend_text,
          col = legend_cols, 
          lty = legend_ltys, 
          lwd = 2, 
          bty = "n",
          cex = 1,
          ncol = 2)
 }
 




# ---------------------------------------------------------------------------- #
# WORK IN PROGRESS ####
# ---------------------------------------------------------------------------- #
 


# Run status quo scenario (Mean F from 2016-2020 applied to projection)
scene_status_quo <- F_clim_sim_scene(
  scene           = scene_full, 
  ssp             = "126", 
  managed_sp      = managed_sp, 
  F_equil         = F_equil, 
  F_zero          = F_zero,
  zero_fishing_sp = NULL,        # No one is zeroed out
  f_ref_yrs       = 2016:2020,
  bioen_sp        = bioen_sp_noceph,
  tdc_hind_bt     = tdc_hind_bt,
  tdr_hind_bt     = tdr_hind_bt,
  verbose         = FALSE
)

run_sq <- rsim.run(scene_status_quo, method = "AB", years = all_years)

# Extract Current Biomass (Mean of 2016-2020 from the hindcast)
# Or use tail(run_sq$annual_Biomass, 1) for the very end of the run
B_current <- apply(run_sq$annual_Biomass[which(all_years %in% 2016:2020), managed_sp], 2, mean)

# 1. Prepare results table
bratio_results <- data.frame(
  Species = managed_sp,
  B_current = B_current[managed_sp],
  B0_2099 = NA,
  B_ratio = NA
)

# 2. Loop to solve for B0 and Ratio
for (i in 1:nrow(bratio_results)) {
  sp <- bratio_results$Species[i]
  
  # Create the B0 scene (F=0 for target species)
  scene_b0 <- F_clim_sim_scene(
    scene           = scene_full, 
    ssp             = "126", 
    managed_sp      = managed_sp, 
    F_equil         = F_equil, 
    F_zero          = F_zero,
    zero_fishing_sp = sp, 
    bioen_sp        = bioen_sp_noceph,
    tdc_hind_bt     = tdc_hind_bt,
    tdr_hind_bt     = tdr_hind_bt,
    verbose         = FALSE
  )
  
  run_b0 <- rsim.run(scene_b0, method = "AB", years = all_years)
  
  # Calculate B0 (Mean of last 10 years)
  b0_val <- mean(tail(run_b0$annual_Biomass[, sp], 10))
  
  # Update table
  bratio_results$B0_2099[i] <- b0_val
  bratio_results$B_ratio[i] <- bratio_results$B_current[i] / b0_val
  
  message(sprintf("Species: %s | Ratio: %.3f", sp, bratio_results$B_ratio[i]))
}

# 3. Final Output
print(bratio_results)
write.csv(bratio_results, "WGOA_Bratio_SSP126.csv", row.names = FALSE)

#What these Ratios mean:$B_{ratio} > 1.0$: The species is currently more abundant 
#than its long-term unexploited equilibrium (rare, usually indicates a massive 
#recruitment pulse or release from predation).
#$0.4 < B_{ratio} < 1.0$: The species is in a "healthy" managed state 
#(many management targets aim for $B_{40\%}$ or $0.4$).
#$B_{ratio} < 0.2$: The species might be considered depleted or overfished relative 
#to its unexploited state.


# 1. Setup
ssp_list   <- c("126", "585")
target_spp <- managed_sp # Or a subset of interest
final_results <- data.frame()

# 2. Outer loop for Climate Scenarios
for (s in ssp_list) {
  message(sprintf("--- Starting SSP %s ---", s))
  
  # A. Run Status Quo for this SSP to get B_current (last 5 years of hindcast)
  scene_sq <- F_clim_sim_scene(
    scene = scene_full, ssp = s, managed_sp = managed_sp, 
    F_equil = F_equil, F_zero = F_zero, f_ref_yrs = 2016:2020,
    bioen_sp = bioen_sp_noceph, tdc_hind_bt = tdc_hind_bt, tdr_hind_bt = tdr_hind_bt,
    verbose = FALSE
  )
  run_sq <- rsim.run(scene_sq, method = "AB", years = all_years)
  b_curr_vec <- apply(run_sq$annual_Biomass[which(all_years %in% 2016:2020), target_spp], 2, mean)
  
  # B. Inner loop for Species B0
  for (sp in target_spp) {
    scene_b0 <- F_clim_sim_scene(
      scene = scene_full, ssp = s, managed_sp = managed_sp, 
      F_equil = F_equil, F_zero = F_zero, zero_fishing_sp = sp, 
      bioen_sp = bioen_sp_noceph, tdc_hind_bt = tdc_hind_bt, tdr_hind_bt = tdr_hind_bt,
      verbose = FALSE
    )
    run_b0 <- rsim.run(scene_b0, method = "AB", years = all_years)
    
    # Calculate values
    b0_val    <- mean(tail(run_b0$annual_Biomass[, sp], 10))
    b_curr    <- b_curr_vec[sp]
    b_ratio   <- b_curr / b0_val
    
    # Store results
    final_results <- rbind(final_results, data.frame(
      SSP = paste0("SSP", s),
      Species = sp,
      B_current = b_curr,
      B0 = b0_val,
      B_ratio = b_ratio
    ))
  }
}

library(ggplot2)

ggplot(final_results, aes(x = Species, y = B_ratio, fill = SSP)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  # Add a horizontal line for the B40% target
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", size = 1) +
  # Add a horizontal line for B100% (Unfished equilibrium)
  geom_hline(yintercept = 1.0, linetype = "solid", color = "darkgrey") +
  scale_fill_manual(values = c("SSP126" = "#2c7bb6", "SSP585" = "#d7191c")) +
  labs(
    title = "Comparison of B-Ratios across Climate Scenarios",
    subtitle = "B_ratio = Current Biomass / Unfished Biomass (B0)",
    y = expression(B[Current] / B[0]),
    x = ""
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 1, y = 0.45, label = "B40% Target", color = "red", hjust = 0)


# POLLOCK ####

# Test a range of HandleSelf values for pollock
pollock_idx <- which(scene_bioen_best$params$spname == "walleye_pollock_adult")

hs_values <- c(0, 0.5, 1, 2)

hs_runs <- purrr::map(hs_values, \(hs) {
  s <- scenes[["126_pcod_zero_all"]]  # start from the F0 scene
  s$params$HandleSelf[pollock_idx] <- hs
  rsim.run(s, method = "AB", years = all_years)
})
names(hs_runs) <- paste0("HandleSelf_", hs_values)

# Plot pollock biomass across HandleSelf values
purrr::map_dfr(names(hs_runs), \(nm) {
  tibble(
    year     = 1991:2099,
    biomass  = hs_runs[[nm]]$annual_Biomass[, "walleye_pollock_adult"],
    scenario = nm
  )
}) |>
  filter(year >= 2021) |>
  ggplot(aes(year, biomass, colour = scenario)) +
  geom_line() +
  geom_hline(yintercept = bal$Biomass["walleye_pollock_adult"],
             linetype = "dashed", alpha = 0.5) +
  labs(title = "SSP126 F0: pollock sensitivity to HandleSelf",
       y = "Biomass (t/km²)", x = NULL)
