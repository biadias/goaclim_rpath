# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 08 June 2026
#
# code/06_run_btarget_ftarget.R
# Purpose: Estimate biomass reference points (B0, Btarget, Blim) and
# climate-informed target fishing mortality (Ftarget) under persist scenario
# for the WGOA fit Rpath model.
# ---------------------------------------------------------------------------- #

#pak::pkg_install("noaa-edab/Rpath@ssb_output_diag")
#library(Rpath)
source("code/01_load_best_model.R")
source("code/Function_F_clim_sim_scene.R")

# Run a quick test
#test_run <- rsim.run(scene_bioen_pcod, method = "AB", years = all_years)

# Check the exact row names the model outputs
#rownames(end_cent_biomass(test_run))
#rownames(end_cent_SSB(test_run))

# Compare them to what you are passing in
#print(managed_sp_list)



# ---------------------------------------------------------------------------- #
# 0. Configuration ####
# ---------------------------------------------------------------------------- #

managed_sp_list <- c(
  "walleye_pollock_adult",
  "pacific_cod_adult",
  "arrowtooth_flounder_adult",
  "flathead_sole_adult",
  "octopus",
  "deep_water_flatfish",
  "sablefish_adult",
  "shallow_water_flatfish",
  "rex_sole_adult",
  "pacific_ocean_perch_adult",
  "slope_rockfish",
  "demersal_shelf_rockfish",
  "pelagic_shelf_rockfish",
  "atka_mackerel",
  "big_skate",
  "longnose_skate",
  "other_skates",
  "pacific_halibut_adult",
  "salmon_shark",
  "pacific_sleeper_shark"
)

tier3_sp_list <- c("walleye_pollock_adult","pacific_cod_adult", "sablefish_adult",
                   "atka_mackerel","slope_rockfish", 
                   #"demersal_shelf_rockfish", # check with Bridget about this one
                   "pelagic_shelf_rockfish",
                   "arrowtooth_flounder_adult","shallow_water_flatfish",
                   "rex_sole_adult","flathead_sole_adult", "deep_water_flatfish" )


ssb_stocks_list <- c("pacific_herring_adult","pacific_halibut_adult",
                "arrowtooth_flounder_adult", "rex_sole_adult",
                "flathead_sole_adult","pacific_ocean_perch_adult",
                "sablefish_adult", "pacific_cod_adult",
                "walleye_pollock_adult")


ssps <- c("persist", "126", "245", "585")

# Annual row indices (rows in annual_Biomass / annual_Catch / ForcedFRate)
hind_rows       <- seq_along(hind_years)                             # 1:30  (1991-2020)
proj_rows       <- (length(hind_years) + 1):length(all_years)       # 31:109 (2021-2099)
ref5_hind_rows   <- (length(hind_years) - 4):length(hind_years)     # 26:30  (2016-2020)

# Output directory
bftarget_dir <- "data/bftarget"
if (!dir.exists(bftarget_dir))
  dir.create(bftarget_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
# 1. Helper functions ####
# ---------------------------------------------------------------------------- #

# Mean end-of-century biomass: mean of the last 60 months (5 years) of
# out_Biomass (monthly output). Column 1 of out_Biomass is a time index;
# species start at column 2. Returns a single-column named matrix (species
# as row names), consistent with EBS indexing: result[species_name, ].
end_cent_biomass <- function(rsim) {
  mean_bio_groups <- matrix(data = NA,
                            nrow = length(2:dim(rsim$out_Biomass)[2]),
                            ncol = 1)
  for (i in 1:(length(2:dim(rsim$out_Biomass)[2]))) {
    mean_bio_groups[i, ] <- mean(rsim$out_Biomass[(dim(rsim$out_Biomass)[1] - 59):dim(rsim$out_Biomass)[1], i + 1])
  }
  row.names(mean_bio_groups) <-
    colnames(rsim$out_Biomass[, 2:dim(rsim$out_Biomass)[2]])
  return(mean_bio_groups)
}

# Mean end-of-century SSB: same last-60-month window as end_cent_biomass but
# reads out_SSB. NOTE: unlike out_Biomass, out_SSB has NO leading time/dummy
# column — columns are directly the stanza "Oldest" species. Do NOT offset by +1.
# Returns a single-column named matrix (species as row names).
end_cent_SSB <- function(rsim) {
  n_spp <- ncol(rsim$out_SSB)
  n_row <- nrow(rsim$out_SSB)
  mean_ssb_groups <- matrix(data = NA, nrow = n_spp, ncol = 1)
  for (i in 1:n_spp) {
    mean_ssb_groups[i, ] <- mean(rsim$out_SSB[(n_row - 59):n_row, i])
  }
  row.names(mean_ssb_groups) <- colnames(rsim$out_SSB)
  return(mean_ssb_groups)
}

# Objective function for cif_wgoa: normalized sum of squared deviations
# between simulated end-of-century biomass/SSB and target biomass/SSB.
#
# ssb_stocks: character vector of species whose reference point is SSB
#             (derived from scene$stanzas$Oldest %in% managed_sp).
#             For those species, end_cent_SSB is used; all others use
#             end_cent_biomass. Pass NULL or character(0) to use biomass only.
sumsq_btarg <- function(F_pars,
                        scene,
                        target_bio,
                        managed_sp,
                        ssb_stocks,
                        all_years,
                        proj_rows) {
  F_pars <- pmax(F_pars, 0)
  scene$fishing$ForcedFRate[proj_rows, managed_sp] <-
    matrix(
      F_pars,
      nrow = length(proj_rows),
      ncol = length(managed_sp),
      byrow = TRUE
    )
  run_tmp <- rsim.run(scene, method = "AB", years = all_years)
  
  # Start with total biomass for all managed species
  sim_bio <- end_cent_biomass(run_tmp)[managed_sp, ]
  
  # Replace SSB-species values with their spawning stock biomass
  if (length(ssb_stocks) > 0) {
    ssb_in_managed <- intersect(ssb_stocks, managed_sp)
    if (length(ssb_in_managed) > 0) {
      sim_bio[ssb_in_managed] <- end_cent_SSB(run_tmp)[ssb_in_managed, ]
    }
  }
  
  sum(((sim_bio - target_bio) / pmax(target_bio, 1e-9))^2, na.rm = TRUE)
}



# Climate-informed F (cif): L-BFGS-B optimization across all managed species
# simultaneously. Starting values are F_meanlast.
#
# ssb_stocks: character vector of stanza "Oldest" species (from
#             intersect(scene$stanzas$Oldest, managed_sp)). The caller should
#             set target_biomass[ssb_stocks] to the SSB-based reference point
#             (e.g., B40_SSB) before calling this function.
cif_wgoa <- function(scene,
                     target_biomass,
                     managed_sp,
                     ssb_stocks,
                     F_meanlast,
                     all_years,
                     proj_rows) {
  pars <- F_meanlast[managed_sp]
  
  optim(
    par        = pars,
    fn         = sumsq_btarg,
    scene      = scene,
    target_bio = target_biomass[managed_sp],
    managed_sp = managed_sp,
    ssb_stocks = ssb_stocks,
    all_years  = all_years,
    proj_rows  = proj_rows,
    method     = "L-BFGS-B",
    lower      = rep(0, length(pars)),
    upper      = rep(max(pars) * 5 + 0.01, length(pars)),
    control    = list(maxit = 200, factr = 1e9)
  )
}

# ---------------------------------------------------------------------------- #
# 2. Hindcast run and F_meanlast ####
# ---------------------------------------------------------------------------- #
# Build a persistence/hindcast-only scene using the fitted best model,
# run through all_years, then derive mean F from hindcast years
# 2016-2020 as catch / biomass.

scene_hind <- F_clim_sim_scene(
  scene      = scene_bioen_best,
  ssp = "persist",
  cons = TRUE,
  resp = TRUE,
  buf = FALSE,
  pcod_rec = TRUE,
  pcod_rec_method = "cauchy",
  bioen_sp = bioen_sp,
  tdc_hind = tdc_hind,
  tdr_hind = tdr_hind,
  managed_sp  = managed_sp_list,
  f_equil     = F_equil,
  f_zero = F_zero,
  f_scenario = "mean",
  f_ref_yrs   = 2016:2020,
  climate_dir = "data/climate/",
  hind_yrs = hind_years,
  proj_yrs = 2021:2099,
  hind_data_start_yr = 1991,
  climate_data_start_yr = 1980,
  verbose = FALSE
)

run_hind <- rsim.run(scene_hind, method = "AB", years = all_years)

# Stanza "Oldest" species that are also in managed_sp_list.
# These species use SSB (not total biomass) as their reference point.
ssb_stocks <- intersect(scene_hind$stanzas$Oldest, managed_sp_list)
message("SSB stocks (stanza Oldest intersect managed_sp_list): ",
        paste(ssb_stocks, collapse = ", "))

# F = catch / biomass, averaged over 2016-2020
catch_ref5 <- run_hind$annual_Catch[ref5_hind_rows, managed_sp_list]
bio_ref5   <- run_hind$annual_Biomass[ref5_hind_rows, managed_sp_list]
F_ref5     <- catch_ref5 / bio_ref5
F_meanlast <- colMeans(F_ref5, na.rm = TRUE)   # kept as F_meanlast for downstream compatibility
F_meanlast[is.nan(F_meanlast) | is.infinite(F_meanlast)] <- 0

# ---------------------------------------------------------------------------- #
# 3. B0 – persistence scenario only ####
# ---------------------------------------------------------------------------- #
# Approach (one species at a time):
#   - Background (non-focal) species: F_meanlast (mean 2016-2020)
#   - Focal managed species:          F = 0
#   - Run to 2100 under stable climate (persist)
#   - Record both end_cent_biomass and end_cent_SSB
#
# This gives B0 for each managed species under climate-stable conditions.

message("\n============================================================")
message("B0 calculation | SSP: persist")
message("============================================================")

# Build persistence climate scene
scene_persist <- F_clim_sim_scene(
  scene      = scene_bioen_best,
  ssp = "persist",
  cons = TRUE,
  resp = TRUE,
  buf = FALSE,
  pcod_rec = TRUE,
  pcod_rec_method = "cauchy",
  bioen_sp = bioen_sp,
  tdc_hind = tdc_hind,
  tdr_hind = tdr_hind,
  managed_sp  = managed_sp_list,
  f_equil     = F_equil,
  f_zero = F_zero,
  f_scenario = "mean",
  f_ref_yrs   = 2016:2020,
  climate_dir = "data/climate/",
  hind_yrs = hind_years,
  proj_yrs = 2021:2099,
  hind_data_start_yr = 1991,
  climate_data_start_yr = 1980,
  verbose = FALSE
)

# Projection period baseline:
#   - All groups:    F_equil  (background low-level F)
#   - Managed groups: F_meanlast (mean 2016-2020), overwritten below per focal species
frate_cols   <- colnames(scene_persist$fishing$ForcedFRate)
equil_cols   <- intersect(frate_cols, names(F_equil))
managed_cols <- intersect(frate_cols, managed_sp_list)

scene_persist$fishing$ForcedFRate[proj_rows, equil_cols] <-
  matrix(
    F_equil[equil_cols],
    nrow = length(proj_rows),
    ncol = length(equil_cols),
    byrow = TRUE
  )
scene_persist$fishing$ForcedFRate[proj_rows, managed_cols] <-
  matrix(
    F_meanlast[managed_cols],
    nrow = length(proj_rows),
    ncol = length(managed_cols),
    byrow = TRUE
  )

# Containers for B0 (biomass and SSB)
B0_biomass_persist <- setNames(numeric(length(managed_sp_list)), managed_sp_list)
B0_SSB_persist     <- setNames(rep(NA, length(managed_sp_list)), managed_sp_list)

ptm <- proc.time()
for (sp in managed_sp_list) {
  message("  B0 persist | F=0 for: ", sp)
  
  # Set focal species F = 0
  scene_persist$fishing$ForcedFRate[proj_rows, sp] <- 0
  
  run_b0 <- rsim.run(scene_persist, method = "AB", years = all_years)
  
  # Record total end-of-century biomass (EVERY species has this)
  B0_biomass_persist[sp] <- end_cent_biomass(run_b0)[sp, ]
  
  # Record end-of-century SSB (ONLY if the species is a multi-stanza stock)
  if (sp %in% ssb_stocks) {
    B0_SSB_persist[sp] <- end_cent_SSB(run_b0)[sp, ]
  }
  
  # Reset focal species F back to status quo 
  # Note: Ensure this variable name matches what you defined earlier! (F_meanlast or F_meanlast10)
  scene_persist$fishing$ForcedFRate[proj_rows, sp] <- F_meanlast[sp]
}
message("persist B0 elapsed: ", round((proc.time() - ptm)[3] / 60, 2), " min")

# Collect into a named list for consistency with downstream reference point code
B0_all <- list(
  persist = list(
    biomass = B0_biomass_persist,
    SSB     = B0_SSB_persist
  )
)

# ---------------------------------------------------------------------------- #
# 4. Biomass reference points – persistence scenario ####
# ---------------------------------------------------------------------------- #
# B40 is the standard target for WGOA groundfish (no crab species).
# Blim rules follow NPFMC FMP conventions -- verify before use in assessments:
#   - Most groundfish:  Blim = 0.05 * B40
#   - Pacific cod:      Blim = 0.50 * B40 (precautionary buffer)
#
# Reference points are derived from B0_biomass (total biomass) and also
# tabulated from B0_SSB (spawning stock biomass) for stocks where SSB is used.

# end of July for the reference points and August the code for HCRs

Btarg_all <- vector("list", 1L)
names(Btarg_all) <- "persist"

ssp <- "persist"

B0     <- B0_all[[ssp]]$biomass
B0_SSB <- B0_all[[ssp]]$SSB

B40 <- 0.40 * B0
B35 <- 0.35 * B0

# Btarget: B40 for all managed groundfish in the WGOA
Btarg <- B40

# Blim (verify against NPFMC FMP guidance before use)
Blim                          <- 0.05 * B40   # 0.02*B0 / 0.05 of B40, 0.02 of B0
Blim["pacific_cod_adult"]     <- 0.50 * B40["pacific_cod_adult"]  # REVIEW FOR COD – Ask Andy ####

target_bio <- cbind(
  B20        = 0.20 * B0,
  B25        = 0.25 * B0,
  B30        = 0.30 * B0,
  B35        = B35, # this is the closes to Fmsy
  B40        = B40,
  B50        = 0.50 * B0,
  B60        = 0.60 * B0,
  B70        = 0.70 * B0,
  B0         = B0,
  B0_SSB     = B0_SSB,
  B40_SSB    = 0.40 * B0_SSB,
  B35_SSB    = 0.35 * B0_SSB,
  Blim       = Blim,
  Btarget_SQ = Btarg #SQ = Status quo target*
)

outfile <- file.path(bftarget_dir, paste0("B_SSB_target_WGOA_GFDL_", ssp, ".csv"))
write.csv(target_bio, file = outfile, row.names = TRUE)
message("Saved: ", outfile)

Btarg_all[[ssp]] <- target_bio

# ---------------------------------------------------------------------------- #
# 5. Climate-informed Ftarget – persistence scenario Fopt_persist ####
# ---------------------------------------------------------------------------- #
# Find the vector of F (one per managed species) that minimizes the sum of
# squared deviations between simulated end-of-century biomass/SSB and Btarget.
#
# Target reference points:
#   - ssb_stocks:       B40_SSB  (spawning stock biomass × 0.40)
#   - all other managed: B40     (total biomass × 0.40)
#
# Starting values = F_meanlast (mean 2016-2020). Optimiser: L-BFGS-B.
#
# Explanation: 
# Fopt_persist:  target fishing mortalities (F) for all 20 managed species. 
# Because optim() was solving for all species simultaneously, Fopt_persist$par 
# is the vector of 20 optimized (managed_species) F values (one for each species) 
# that collectively minimized the overall difference between the simulated 
# biomasses and the B40 targets.

message("\n============================================================")
message("Ftarget optimisation | SSP: persist")
message("============================================================")

scene_persist_fopt <- F_clim_sim_scene(
  scene      = scene_bioen_best,
  ssp = "persist",
  cons = TRUE,
  resp = TRUE,
  buf = FALSE,
  pcod_rec = TRUE,
  pcod_rec_method = "cauchy",
  bioen_sp = bioen_sp,
  tdc_hind = tdc_hind,
  tdr_hind = tdr_hind,
  managed_sp  = managed_sp_list,
  f_equil     = F_equil,
  f_zero = F_zero,
  f_scenario = "mean",
  f_ref_yrs   = 2016:2020,
  climate_dir = "data/climate/",
  hind_yrs = hind_years,
  proj_yrs = 2021:2099,
  hind_data_start_yr = 1991,
  climate_data_start_yr = 1980,
  verbose = FALSE
)

# Projection period baseline: F_equil for all groups, F_meanlast for managed
frate_cols   <- colnames(scene_persist_fopt$fishing$ForcedFRate)
equil_cols   <- intersect(frate_cols, names(F_equil))
managed_cols <- intersect(frate_cols, managed_sp_list)

scene_persist_fopt$fishing$ForcedFRate[proj_rows, equil_cols] <-
  matrix(
    F_equil[equil_cols],
    nrow = length(proj_rows),
    ncol = length(equil_cols),
    byrow = TRUE
  )
scene_persist_fopt$fishing$ForcedFRate[proj_rows, managed_cols] <-
  matrix(
    F_meanlast[managed_cols],
    nrow = length(proj_rows),
    ncol = length(managed_cols),
    byrow = TRUE
  )

# Build target biomass vector:
#   - ssb_stocks get B40_SSB (SSB-based reference point)
#   - all other managed species get B40 (total biomass reference point)
target_biomass_fopt                <- Btarg_all[["persist"]][, "B40"]
target_biomass_fopt[ssb_stocks]    <- Btarg_all[["persist"]][ssb_stocks, "B40_SSB"]

ptm <- proc.time()
Fopt_persist <- cif_wgoa(
  scene          = scene_persist_fopt,
  target_biomass = target_biomass_fopt,
  managed_sp     = managed_sp_list,
  ssb_stocks     = ssb_stocks,
  F_meanlast     = F_meanlast,
  all_years      = all_years,
  proj_rows      = proj_rows
)
message("persist Ftarget elapsed: ", round((proc.time() - ptm)[3] / 60, 2), " min")
message("Convergence: ", Fopt_persist$convergence, "  (0 = success)")

# ---------------------------------------------------------------------------- #
# 6. Save Ftarget results simultaneously ####
# ---------------------------------------------------------------------------- #

Ftarg_matrix <- matrix(
  Fopt_persist$par,
  ncol = 1,
  dimnames = list(managed_sp_list, "persist")
)

write.csv(
  Ftarg_matrix,
  file = file.path(bftarget_dir, "Fopt_target_WGOA_GFDL_persist.csv"),
  row.names = TRUE
)
message("Saved: ", file.path(bftarget_dir, "Fopt_target_WGOA_GFDL_persist.csv"))

# Quick check: compare Ftarget to F_meanlast
Fcomp <- data.frame(
  F_meanlast    = F_meanlast[managed_sp_list],
  Ftarg_persist   = Ftarg_matrix[, "persist"]
)
print(round(Fcomp, 4))


#NOTES ####
# Observation error from Andy code! Ask Andy about it.
# Single run, but I also will have the observation error/variance that Andy developed on it.
# One B0 under the persist scenario (Climate stable)

# ---------------------------------------------------------------------------- #
# 7.Climate-informed Ftarget– persistence scenario (Individual sp Loop)#####
# ---------------------------------------------------------------------------- #
# Find the F for EACH managed species individually that drives its end-of-century 
# biomass/SSB to its specific B40 and B35 target, while all other species are held 
# at status quo (F_meanlast) or background level (F_equil).

message("\n============================================================")
message("Ftarget individual optimization | SSP: persist")
message("============================================================")

# Helper function for single-species optimization
sumsq_btarg_single <- function(F_val,
                               sp,
                               scene,
                               target_bio_val,
                               is_ssb,
                               all_years,
                               proj_rows) {
  # Non-negative F, just in case
  F_val <- max(F_val, 0)
  
  # Set F only for the focal species
  scene$fishing$ForcedFRate[proj_rows, sp] <- F_val
  
  # Run the simulation
  run_tmp <- rsim.run(scene, method = "AB", years = all_years)
  
  # Extract simulated end-of-century metric for the focal species
  if (is_ssb) {
    sim_val <- end_cent_SSB(run_tmp)[sp, ]
  } else {
    sim_val <- end_cent_biomass(run_tmp)[sp, ]
  }
  
  # Return squared relative difference
  return(((sim_val - target_bio_val) / max(target_bio_val, 1e-9))^2)
}


# Initialize the baseline scenario
scene_persist_f40 <- F_clim_sim_scene(
  scene                 = scene_bioen_best,
  ssp                   = "persist",
  cons                  = TRUE,
  resp                  = TRUE,
  buf                   = FALSE,
  pcod_rec              = TRUE,
  pcod_rec_method       = "cauchy",
  bioen_sp              = bioen_sp,
  tdc_hind              = tdc_hind,
  tdr_hind              = tdr_hind,
  managed_sp            = managed_sp_list,
  f_equil               = F_equil,
  f_zero                = F_zero,
  f_scenario            = "mean",
  f_ref_yrs             = 2016:2020,
  climate_dir           = "data/climate/",
  hind_yrs              = hind_years,
  proj_yrs              = 2021:2099,
  hind_data_start_yr    = 1991,
  climate_data_start_yr = 1980,
  verbose               = FALSE
)



# Set baseline projection rates: background F for non-managed, status quo F for managed
frate_cols   <- colnames(scene_persist_f40$fishing$ForcedFRate)
equil_cols   <- intersect(frate_cols, names(F_equil))
managed_cols <- intersect(frate_cols, managed_sp_list)

scene_persist_f40$fishing$ForcedFRate[proj_rows, equil_cols] <-
  matrix(F_equil[equil_cols], nrow = length(proj_rows), ncol = length(equil_cols), byrow = TRUE)

scene_persist_f40$fishing$ForcedFRate[proj_rows, managed_cols] <-
  matrix(F_meanlast10[managed_cols], nrow = length(proj_rows), ncol = length(managed_cols), byrow = TRUE)


## 7.1 optimize F40 (target Fishing Mortality) Run the individual optimization loop ####
F40_persist <- setNames(numeric(length(managed_sp_list)), managed_sp_list)

ptm <- proc.time()
for (sp in managed_sp_list) {
  message("  Optimizing F40 for: ", sp)
  
  # Determine if this species uses SSB or total biomass
  is_ssb <- sp %in% ssb_stocks
  
  # Get target value (B40_SSB or B40)
  target_val <- if (is_ssb) {
    Btarg_all[["persist"]][sp, "B40_SSB"]
  } else {
    Btarg_all[["persist"]][sp, "B40"]
  }
  
  # Establish a reasonable upper search bound for F (5x status quo F, with a floor of 2.0)
  upper_bound <- F_meanlast[sp] * 5 + 0.1
  if (upper_bound < 0.2) {
    upper_bound <- 2.0  # Safe window if status quo F is near zero
  }
  
  # Run 1D optimization
  opt_res <- optimize( # 1D optimization with optimize() instead of optim() 
    f              = sumsq_btarg_single,
    interval       = c(0, upper_bound),
    sp             = sp,
    scene          = scene_persist_f40,
    target_bio_val = target_val,
    is_ssb         = is_ssb,
    all_years      = all_years,
    proj_rows      = proj_rows,
    tol            = 1e-4
  )
  
  # Record results
  F40_persist[sp] <- opt_res$minimum
  
  # RESET focal species F back to status quo (F_meanlast) before moving to next species
  scene_persist_f40$fishing$ForcedFRate[proj_rows, sp] <- F_meanlast[sp]
}

message("persist F40 individual optimization elapsed: ", round((proc.time() - ptm)[3] / 60, 2), " min")
  
## 7.2 optimize F35 = OFL target Fishing Mortality. Run the individual optimization loop ####

F35_persist <- setNames(numeric(length(managed_sp_list)), managed_sp_list)

ptm <- proc.time()
for (sp in managed_sp_list) {
  message("  Optimizing F35 for: ", sp)
  is_ssb <- sp %in% ssb_stocks
  
  target_val_35 <- if (is_ssb) {
    Btarg_all[["persist"]][sp, "B35_SSB"]
  } else {
    Btarg_all[["persist"]][sp, "B35"]
  }
  
  upper_bound <- F_meanlast[sp] * 5 + 0.1
  if (upper_bound < 0.2) {
    upper_bound <- 2.0  # Safe window if status quo F is near zero
  }
  
  opt_res_35 <- optimize(
    f              = sumsq_btarg_single,
    interval       = c(0, upper_bound),
    sp             = sp,
    scene          = scene_persist_f40,
    target_bio_val = target_val_35,
    is_ssb         = is_ssb,
    all_years      = all_years,
    proj_rows      = proj_rows,
    tol            = 1e-4
  )
  
  F35_persist[sp] <- opt_res_35$minimum
  scene_persist_f40$fishing$ForcedFRate[proj_rows, sp] <- F_meanlast[sp]
}
message("persist F35 individual optimization elapsed: ", round((proc.time() - ptm)[3] / 60, 2), " min")

# ---------------------------------------------------------------------------- #
# 8. Save Ftarget and Tier 3 F results ####
# ---------------------------------------------------------------------------- #

Tier3_F_matrix <- data.frame(
  F40_ABC = F40_persist,
  F35_OFL = F35_persist
)

write.csv(
  Tier3_F_matrix,
  file = file.path(bftarget_dir, "F_Tier3_WGOA_GFDL_persist.csv"),
  row.names = TRUE
)
message("Saved: ", file.path(bftarget_dir, "F_Tier3_WGOA_GFDL_persist.csv"))

# Quick check comparison
Fcomp <- data.frame(
  F_meanlast  = F_meanlast[managed_sp_list],
  F40_ABC       = F40_persist,
  F35_OFL       = F35_persist
)
print(round(Fcomp, 4))

