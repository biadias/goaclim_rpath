# ---------------------------------------------------------------------------- #
# AUTHORS: Bia Dias
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
# DATE: 08 June 2026
#
# code/06_run_btarget_ftarget.R
# Purpose: Estimate biomass reference points (B0, Btarget, Blim) and
# climate-informed target fishing mortality (Ftarget) under GFDL SSP 126,
# SSP 245, and SSP 585, plus a persistence scenario for the WGOA fit Rpath model.
# ---------------------------------------------------------------------------- #

source("code/01_load_best_model.R")
source("code/Function_F_clim_sim_scene.R")

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

ssps <- c("persist", "126", "245", "585")

# Annual row indices (rows in annual_Biomass / annual_Catch / ForcedFRate)
hind_rows       <- seq_along(hind_years)                             # 1:30  (1991-2020)
proj_rows       <- (length(hind_years) + 1):length(all_years)       # 31:109 (2021-2099)
last10_hind_rows <- (length(hind_years) - 9):length(hind_years)     # 21:30  (2011-2020)

# Output directory
bftarget_dir <- "data/bftarget"
if (!dir.exists(bftarget_dir)) dir.create(bftarget_dir, recursive = TRUE)

# ---------------------------------------------------------------------------- #
# 1. Helper functions ####
# ---------------------------------------------------------------------------- #

# Mean end-of-century biomass: mean of the last 60 months (5 years) of
# out_Biomass (monthly output). Column 1 of out_Biomass is a time index;
# species start at column 2. Returns a single-column named matrix (species
# as row names), consistent with EBS indexing: result[species_name, ].
end_cent_biomass <- function(rsim) {
  mean_bio_groups <- matrix(
    data = NA,
    nrow = length(2:dim(rsim$out_Biomass)[2]),
    ncol = 1
  )
  for (i in 1:(length(2:dim(rsim$out_Biomass)[2]))) {
    mean_bio_groups[i, ] <- mean(
      rsim$out_Biomass[(dim(rsim$out_Biomass)[1] - 59):dim(rsim$out_Biomass)[1], i + 1]
    )
  }
  row.names(mean_bio_groups) <-
    colnames(rsim$out_Biomass[, 2:dim(rsim$out_Biomass)[2]])
  return(mean_bio_groups)
}

# Objective function for cif_wgoa: normalized sum of squared deviations
# between simulated end-of-century biomass and target biomass
sumsq_btarg <- function(F_pars, scene, target_bio, managed_sp, all_years, proj_rows) {
  F_pars <- pmax(F_pars, 0)
  scene$fishing$ForcedFRate[proj_rows, managed_sp] <-
    matrix(F_pars, nrow = length(proj_rows), ncol = length(managed_sp), byrow = TRUE)
  run_tmp <- rsim.run(scene, method = "AB", years = all_years)
  sim_bio <- end_cent_biomass(run_tmp)[managed_sp, ]   # named vector from matrix rows
  sum(((sim_bio - target_bio) / pmax(target_bio, 1e-9))^2, na.rm = TRUE)
}

# Climate-informed F (cif): L-BFGS-B optimization across all managed species
# simultaneously. Starting values are F_meanlast10.
cif_wgoa <- function(scene, target_biomass, managed_sp,
                     F_meanlast10, all_years, proj_rows) {
  pars <- F_meanlast10[managed_sp]

  optim(
    par     = pars,
    fn      = sumsq_btarg,
    scene   = scene,
    target_bio = target_biomass[managed_sp],
    managed_sp = managed_sp,
    all_years  = all_years,
    proj_rows  = proj_rows,
    method  = "L-BFGS-B",
    lower   = rep(0,       length(pars)),
    upper   = rep(max(pars) * 5 + 0.01, length(pars)),
    control = list(maxit = 200, factr = 1e9)
  )
}

# ---------------------------------------------------------------------------- #
# 2. Hindcast run and F_meanlast10 ####
# ---------------------------------------------------------------------------- #
# Build a persistence/hindcast-only scene using the fitted best model,
# run through all_years, then derive mean F from the last 10 hindcast years
# (2011-2020) as catch / biomass.

scene_hind <- F_clim_sim_scene(
  scene      = scene_bioen_best, ssp = "persist",
  cons = TRUE, resp = TRUE, buf = FALSE,
  pcod_rec = TRUE, pcod_rec_method = "cauchy",
  bioen_sp = bioen_sp, tdc_hind = tdc_hind, tdr_hind = tdr_hind,
  managed_sp  = managed_sp_list,
  f_equil     = F_equil, f_zero = F_zero, f_scenario = "mean",
  f_ref_yrs   = 2016:2020, #change that!!!
  climate_dir = "data/climate/",
  hind_yrs = hind_years, proj_yrs = 2021:2099,
  hind_data_start_yr = 1991, climate_data_start_yr = 1980,
  verbose = FALSE
)

run_hind <- rsim.run(scene_hind, method = "AB", years = all_years)

# F = catch / biomass, averaged over 2011-2020
catch_last10 <- run_hind$annual_Catch[last10_hind_rows,   managed_sp_list]
bio_last10   <- run_hind$annual_Biomass[last10_hind_rows, managed_sp_list]
F_last10     <- catch_last10 / bio_last10
F_meanlast10 <- colMeans(F_last10, na.rm = TRUE)
F_meanlast10[is.nan(F_meanlast10) | is.infinite(F_meanlast10)] <- 0

# ---------------------------------------------------------------------------- #
# 3. B0 per SSP #### change this! 
# ---------------------------------------------------------------------------- #
# For each SSP, loop over managed species:
#   - all other species remain at F_meanlast10 in the projection period
#   - set F = 0 for the focal species, run the sim, record end-century biomass
#   - that biomass is B0 for the focal species under that climate scenario 

#NOTES: background F (low level F, fixed F for everything that was not focused)
#F=0 for the managed species and Fmean2016-2020 for the other species. 
#We want the SSBout 
#run to 2100
#F=0 one at the time with Fmean for the other species. 
#F40 is the same approach than described above. 

B0_all <- vector("list", length(ssps))
names(B0_all) <- ssps

for (ssp in ssps) {
  message("\n============================================================")
  message("B0 calculation | SSP: ", ssp)
  message("============================================================")

  # Build climate scene for this SSP
  scene_ssp <- F_clim_sim_scene(
    scene      = scene_bioen_best, ssp = ssp,
    cons = TRUE, resp = TRUE, buf = FALSE,
    pcod_rec = TRUE, pcod_rec_method = "cauchy",
    bioen_sp = bioen_sp, tdc_hind = tdc_hind, tdr_hind = tdr_hind,
    managed_sp  = managed_sp_list,
    f_equil     = F_equil, f_zero = F_zero, f_scenario = "mean",
    f_ref_yrs   = 2016:2020,
    climate_dir = "data/climate/",
    hind_yrs = hind_years, proj_yrs = 2021:2099,
    hind_data_start_yr = 1991, climate_data_start_yr = 1980,
    verbose = FALSE
  )

  # Projection period: F_equil for all groups, then F_meanlast10 for managed
  frate_cols   <- colnames(scene_ssp$fishing$ForcedFRate)
  equil_cols   <- intersect(frate_cols, names(F_equil))
  managed_cols <- intersect(frate_cols, managed_sp_list)

  scene_ssp$fishing$ForcedFRate[proj_rows, equil_cols] <-
    matrix(F_equil[equil_cols],
           nrow = length(proj_rows), ncol = length(equil_cols), byrow = TRUE)
  scene_ssp$fishing$ForcedFRate[proj_rows, managed_cols] <-
    matrix(F_meanlast10[managed_cols],
           nrow = length(proj_rows), ncol = length(managed_cols), byrow = TRUE)

  B0_vec <- setNames(numeric(length(managed_sp_list)), managed_sp_list)

  ptm <- proc.time()
  for (sp in managed_sp_list) {
    message("  B0 for: ", sp)
    scene_ssp$fishing$ForcedFRate[proj_rows, sp] <- 0          # set focal F = 0
    run_b0     <- rsim.run(scene_ssp, method = "AB", years = all_years)
    B0_vec[sp] <- end_cent_biomass(run_b0)[sp, ]   # [species, ] from named matrix
    scene_ssp$fishing$ForcedFRate[proj_rows, sp] <- F_meanlast10[sp]  # reset
  }
  message("SSP ", ssp, " B0 elapsed: ",
          round((proc.time() - ptm)[3] / 60, 2), " min")

  B0_all[[ssp]] <- B0_vec
}

# ---------------------------------------------------------------------------- #
# 4. Biomass reference points per SSP ####
# ---------------------------------------------------------------------------- #
# B40 is the standard target for WGOA groundfish (no crab species).
# Blim rules follow NPFMC FMP conventions -- verify before use in assessments:
#   - Most groundfish:  Blim = 0.05 * B40
#   - Pacific cod:      Blim = 0.50 * B40 (precautionary buffer)

#B0 stable climate only. 
# end of July for the reference points and August the code for HCRs

Btarg_all <- vector("list", length(ssps))
names(Btarg_all) <- ssps

for (ssp in ssps) {
  B0 <- B0_all[[ssp]]

  B40 <- 0.40 * B0
  B35 <- 0.35 * B0

  # Btarget: B40 for all managed groundfish in the WGOA
  Btarg <- B40

  # Blim (verify against NPFMC FMP guidance before use)
  Blim                          <- 0.05 * B40 #0.02*B0/ 0.05 of B40, 0.02 of B0
  Blim["pacific_cod_adult"]     <- 0.50 * B40["pacific_cod_adult"] # REVIEW THIS FOR COD - Ask Andy####

  target_bio <- cbind(
    B20           = 0.20 * B0,
    B25           = 0.25 * B0,
    B30           = 0.30 * B0,
    B35           = B35,
    B40           = B40,
    B50           = 0.50 * B0,
    B60           = 0.60 * B0,
    B70           = 0.70 * B0,
    B0            = B0,
    Blim          = Blim,
    Btarget_SQ    = Btarg
  )

  outfile <- file.path(bftarget_dir, paste0("B_target_WGOA_GFDL_", ssp, ".csv"))
  write.csv(target_bio, file = outfile, row.names = TRUE)
  message("Saved: ", outfile)

  Btarg_all[[ssp]] <- target_bio
}

# ---------------------------------------------------------------------------- #
# 5. Climate-informed Ftarget per SSP ####
# ---------------------------------------------------------------------------- #
# For each SSP, find the vector of fishing mortalities (one per managed species)
# that minimizes the sum of squared deviations between simulated end-of-century
# biomass and Btarget (= B40 for all species here).
# Starting values = F_meanlast10. Uses L-BFGS-B via cif_wgoa().

#review this: have one species at time. 


Fopt_all <- vector("list", length(ssps))
names(Fopt_all) <- ssps

for (ssp in ssps) {
  message("\n============================================================")
  message("Ftarget optimisation | SSP: ", ssp)
  message("============================================================")

  scene_ssp <- F_clim_sim_scene(
    scene      = scene_bioen_best, ssp = ssp,
    cons = TRUE, resp = TRUE, buf = FALSE,
    pcod_rec = TRUE, pcod_rec_method = "cauchy",
    bioen_sp = bioen_sp, tdc_hind = tdc_hind, tdr_hind = tdr_hind,
    managed_sp  = managed_sp_list,
    f_equil     = F_equil, f_zero = F_zero, f_scenario = "mean",
    f_ref_yrs   = 2016:2020,
    climate_dir = "data/climate/",
    hind_yrs = hind_years, proj_yrs = 2021:2099,
    hind_data_start_yr = 1991, climate_data_start_yr = 1980,
    verbose = FALSE
  )

  # Projection period baseline: F_equil for all, F_meanlast10 for managed
  frate_cols   <- colnames(scene_ssp$fishing$ForcedFRate)
  equil_cols   <- intersect(frate_cols, names(F_equil))
  managed_cols <- intersect(frate_cols, managed_sp_list)

  scene_ssp$fishing$ForcedFRate[proj_rows, equil_cols] <-
    matrix(F_equil[equil_cols],
           nrow = length(proj_rows), ncol = length(equil_cols), byrow = TRUE)
  scene_ssp$fishing$ForcedFRate[proj_rows, managed_cols] <-
    matrix(F_meanlast10[managed_cols],
           nrow = length(proj_rows), ncol = length(managed_cols), byrow = TRUE)

  target_biomass <- Btarg_all[[ssp]][, "Btarget_SQ"]

  ptm <- proc.time()
  Fopt_all[[ssp]] <- cif_wgoa(
    scene          = scene_ssp,
    target_biomass = target_biomass,
    managed_sp     = managed_sp_list,
    F_meanlast10   = F_meanlast10,
    all_years      = all_years,
    proj_rows      = proj_rows
  )
  message("SSP ", ssp, " Ftarget elapsed: ",
          round((proc.time() - ptm)[3] / 60, 2), " min")
  message("Convergence: ", Fopt_all[[ssp]]$convergence,
          "  (0 = success)")
}

# ---------------------------------------------------------------------------- #
# 6. Save Ftarget results ####
# ---------------------------------------------------------------------------- #

Ftarg_matrix <- do.call(cbind, lapply(ssps, function(ssp) Fopt_all[[ssp]]$par))
colnames(Ftarg_matrix) <- ssps
rownames(Ftarg_matrix) <- managed_sp_list

write.csv(Ftarg_matrix,
          file = file.path(bftarget_dir, "Ftarget_WGOA_GFDL_allSSPs.csv"),
          row.names = TRUE)
message("Saved: ", file.path(bftarget_dir, "Ftarget_WGOA_GFDL_allSSPs.csv"))

# Quick check: compare Ftarget to F_meanlast10
Fcomp <- data.frame(
  F_meanlast10 = F_meanlast10[managed_sp_list],
  Ftarg_matrix
)
print(round(Fcomp, 4))


#NOTES ####
# Observation error from Andy code! Ask Andy about it. 
# Single run, but I also will have the observation error/variance that Andy developed on it.
# One B0 under the persist scenario (Climate stable) 
# 
