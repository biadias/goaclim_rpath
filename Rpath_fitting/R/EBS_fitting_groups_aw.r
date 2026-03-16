# ---------------------------------------------------------------------------- #
# AUTHORS: George (Andy) Whitehouse
# AFFILIATIONS: CICOES University of Washington
# E-MAIL OF CORRESPONDENCE AUTHOR: gaw@uw.edu
# UPDATED: 1 August 2025
#
# Purpose: Script for initial fitting exploration of EBS Rpath for ACLIM 3.0. Run
# Rpath_fitting/R/EBS_fitting_setup.r to set up base scenario object.
# ---------------------------------------------------------------------------- #

# General instructions:
# Fit vulnerabilities from low to high trophic level.
# First, fit each group separately for vulnerabilities (get realistic results?). 
# Record the change in NLL and qualitatively describe change in fit of the group 
# being fit.
# second, work up from lower to upper trophic levels fitting and retaining 
# vulnerabilities. Again, note change in NLL.
# ---------------------------------------------------------------------------- #

# ebs abbreviated name lookup
ebs_name_lookup <- read.csv("C:/Users/andy.whitehouse/Work/src/aclim3_rpath/a3_group_3letter_lookup.csv")
source("R/EBS_fitting_setup.r")
# base model without fitting any parameters ---------------------------------- #
# we want our outputs to plot all the species - make a list 
plot.species <- c(rpath.living(bal),rpath.detrital(bal))  

# Run and calculate without changing any fitting values  
run.base <- rsim.fit.run(NA, NA, NA, scene=scene_fit, run_method="AB", 
                         run_years=hind_years, verbose=T)
# test to show output negative log likelihood
run.base.nll <- rsim.fit.run(NA, NA, NA, scene=scene_fit, run_method="AB", 
                             run_years=hind_years, verbose=F); run.base.nll # default verbosity

# Save likelihood table                            
base.like <- rsim.fit.table(scene_fit, run.base)

# Plot panels    
rsim.runplot(scene_fit, run.base, plot.species)
rsim.runplot(scene_fit, run.base, "euphausiids")
rsim.runplot(scene_fit, run.base, "copepods")
# ---------------------------------------------------------------------------- #  


# fitting groups one at a time and storing their VVs and NLLs ---------------- #
source("R/vv_fit_loop.R")
ptm <- proc.time()
# test_vv_loop <- vv_fit_loop(bal, scene_fit, hind_years, run.base.nll)
vv_loop_7605 <- vv_fit_loop(bal, scene_fit, hind_years, run.base.nll)
proc.time() - ptm
# this took about 8 minutes
round(vv_loop_7605[[1]][order(vv_loop_7605[[1]][,5]),1:7],3)
write.csv(vv_loop_7605[[2]], "C:/Users/andy.whitehouse/Work/Andy/REEM/ACLIM/ACLIM_3/fitting/vv_loop_7605_change_nlls.csv")

# exponentiate the VVs
vv_loop_7605[[1]]$Vprey <- 1 + exp(vv_loop_7605[[1]]$preyvul)
vv_loop_7605[[1]]$Vpred <- 1 + exp(vv_loop_7605[[1]]$predvul)
write.csv(vv_loop_7605[[1]], "C:/Users/andy.whitehouse/Work/Andy/REEM/ACLIM/ACLIM_3/fitting/vv_loop_7605_initial_notes.csv")
round(vv_loop_7605[[1]][order(vv_loop_7605[[1]][,5]),c(1:7,9:10)],3)

# apply all the vulnerabilities to a new scenario object --------------------- #
# all bio VVs
# all_bio_vv <- c(test_vv_loop$preyvul, test_vv_loop$predvul)
# all_bio_sp <- rep(row.names(test_vv_loop),2)
# all_bio_vartype <- c(rep("preyvul", length(row.names(test_vv_loop))), 
#                      rep("predvul",length(row.names(test_vv_loop)))) 
# scene.all_bio <- scene_fit
# scene.all_bio$params <- rsim.fit.apply(all_bio_vv, all_bio_sp, 
#                                        all_bio_vartype, scene_fit$params)
# # Test that we get the correct likelihood in the new scene
# all_bio_nll <- rsim.fit.run(NA, NA, NA, scene = scene.all_bio, run_method="AB", 
#                             run_years = hind_years, verbose=F); all_bio_nll  
# all_bio_nll_out <- rsim.fit.run(NA, NA, NA, scene = scene.all_bio, run_method="AB", 
#                                 run_years = hind_years, verbose=T)  
# # Now let's look
# rsim.runplot(scene.all_bio, all_bio_nll_out, plot.species)
# 
# all_bio.like <- rsim.fit.table(scene.all_bio, all_bio_nll_out)  
# # what changed?            
# # View(all_bio.like - base.like)
# like.diff.all_bio <- all_bio.like - base.like
# like.diff.all_bio[order(like.diff.all_bio$Biomass),]
# all_bio_nll - run.base.nll
# ---------------------------------------------------------------------------- #
# apply the VVs one-at-a-time and run

for(i in 1:dim(vv_loop_7605[[1]])[1]) {
  # Apply the parameters to a new scene ("permanently").  Note: this loses the
  # distinction between preyvul and predvul, unfortunately.
  scene.group        <- scene_fit
  group_vartype      <- c("preyvul", "predvul") 
  group.result       <- as.numeric(vv_loop_7605[[1]][i, 2:3])
  func_group         <- row.names(vv_loop_7605[[1]])[i]
  group.species      <- rep(func_group,2)
  scene.group$params <- rsim.fit.apply(group.result, group.species, 
                                       group_vartype, scene.group$params)
  # Test that we get the correct likelihood in the new scene
  group_nll <- rsim.fit.run(NA, NA, NA, scene = scene.group, run_method = "AB", 
                            run_years = hind_years, verbose=F)  
  print(paste0(row.names(vv_loop_7605)[i], " ", group_nll))
  # Save the actual full run (same as above with verbose=T)
  run.group <- rsim.fit.run(group.result, group.species, 
                              group_vartype, scene = scene_fit, 
                              run_method = "AB", run_years = hind_years, 
                              verbose=T)
  # save and rename scenario object
  assign(paste0("scene_", ebs_name_lookup$abbreviation[i]), scene.group)
  assign(paste0("run_", ebs_name_lookup$abbreviation[i]), run.group)
  
}
# Now let's look
# gray_whales
rsim.runplot(scene_grw, run_grw, "gray_whales")
rsim.runplot(scene_fit, run.base, "gray_whales")
#nf_seal_adu
rsim.runplot(scene_nfa, run_nfa, "nf_seal_adu")
rsim.runplot(scene_fit, run.base, "nf_seal_adu")
# oth_birds
rsim.runplot(scene_obd, run_obd, "oth_birds")
rsim.runplot(scene_fit, run.base, "oth_birds")
# murres_puffins
rsim.runplot(scene_mpf, run_mpf, "murres_puffins")
rsim.runplot(scene_fit, run.base, "murres_puffins")
# kittiwakes
rsim.runplot(scene_ktw, run_ktw, "kittiwakes")
rsim.runplot(scene_fit, run.base, "kittiwakes")
# auklets
rsim.runplot(scene_auk, run_auk, "auklets")
rsim.runplot(scene_fit, run.base, "auklets")
# fulmars
rsim.runplot(scene_ful, run_ful, "fulmars")
rsim.runplot(scene_fit, run.base, "fulmars")
# sharks
rsim.runplot(scene_srk, run_srk, "sharks")
rsim.runplot(scene_fit, run.base, "sharks")
# pollock_juv
rsim.runplot(scene_plj, run_plj, "pollock_juv")
rsim.runplot(scene_fit, run.base, "pollock_juv")
# pollock_adu
rsim.runplot(scene_pol, run_pol, "pollock_adu")
rsim.runplot(scene_fit, run.base, "pollock_adu")
# pcod_juv
rsim.runplot(scene_coj, run_coj, "pcod_juv")
rsim.runplot(scene_fit, run.base, "pcod_juv")
# pcod_adu
rsim.runplot(scene_cod, run_cod, "pcod_adu")
rsim.runplot(scene_fit, run.base, "pcod_adu")
# herring
rsim.runplot(scene_her, run_her, "herring")
rsim.runplot(scene_fit, run.base, "herring")
# arrowtooth_juv
rsim.runplot(scene_atj, run_atj, "arrowtooth_juv")
rsim.runplot(scene_fit, run.base, "arrowtooth_juv")
# arrowtooth_adu
rsim.runplot(scene_atf, run_atf, "arrowtooth_adu")
rsim.runplot(scene_fit, run.base, "arrowtooth_adu")
# kamchatka
rsim.runplot(scene_kam, run_kam, "kamchatka")
rsim.runplot(scene_fit, run.base, "kamchatka")
# gr_turbot_juv
rsim.runplot(scene_gtj, run_gtj, "gr_turbot_juv")
rsim.runplot(scene_fit, run.base, "gr_turbot_juv")
# gr_turbot_adu
rsim.runplot(scene_gta, run_gta, "gr_turbot_adu")
rsim.runplot(scene_fit, run.base, "gr_turbot_adu")
# halibut_juv
rsim.runplot(scene_haj, run_haj, "halibut_juv")
rsim.runplot(scene_fit, run.base, "halibut_juv")
# halibut_adu
rsim.runplot(scene_hal, run_hal, "halibut_adu")
rsim.runplot(scene_fit, run.base, "halibut_adu")
# yf_sole
rsim.runplot(scene_yfs, run_yfs, "yf_sole")
rsim.runplot(scene_fit, run.base, "yf_sole")
# fh_sole
rsim.runplot(scene_fhs, run_fhs, "fh_sole")
rsim.runplot(scene_fit, run.base, "fh_sole")
# nr_sole
rsim.runplot(scene_nrs, run_nrs, "nr_sole")
rsim.runplot(scene_fit, run.base, "nr_sole")
# ak_plaice
rsim.runplot(scene_akp, run_akp, "ak_plaice")
rsim.runplot(scene_fit, run.base, "ak_plaice")
# oth_flatfish
rsim.runplot(scene_off, run_off, "oth_flatfish")
rsim.runplot(scene_fit, run.base, "oth_flatfish")
# ak_skate
rsim.runplot(scene_aks, run_aks, "ak_skate")
rsim.runplot(scene_fit, run.base, "ak_skate")
# oth_skate
rsim.runplot(scene_osk, run_osk, "oth_skate")
rsim.runplot(scene_fit, run.base, "oth_skate")
# sablefish_juv
rsim.runplot(scene_sbj, run_sbj, "sablefish_juv")
rsim.runplot(scene_fit, run.base, "sablefish_juv")
# sablefish_adu
rsim.runplot(scene_sba, run_sba, "sablefish_adu")
rsim.runplot(scene_fit, run.base, "sablefish_adu")
# eelpouts
rsim.runplot(scene_eel, run_eel, "eelpouts")
rsim.runplot(scene_fit, run.base, "eelpouts")
# deep_demersals
rsim.runplot(scene_ddm, run_ddm, "deep_demersals")
rsim.runplot(scene_fit, run.base, "deep_demersals")
# pac_ocean_perch
rsim.runplot(scene_pop, run_pop, "pac_ocean_perch")
rsim.runplot(scene_fit, run.base, "pac_ocean_perch")
# oth_rockfish
rsim.runplot(scene_orf, run_orf, "oth_rockfish")
rsim.runplot(scene_fit, run.base, "oth_rockfish")
# north_rockfish
rsim.runplot(scene_nrf, run_nrf, "north_rockfish")
rsim.runplot(scene_fit, run.base, "north_rockfish")
# shortraker_rock
rsim.runplot(scene_shr, run_shr, "shortraker_rock")
rsim.runplot(scene_fit, run.base, "shortraker_rock")
# rougheye_rock
rsim.runplot(scene_rer, run_rer, "rougheye_rock")
rsim.runplot(scene_fit, run.base, "rougheye_rock")
# atka
rsim.runplot(scene_atk, run_atk, "atka")
rsim.runplot(scene_fit, run.base, "atka")
# shallow_demersals
rsim.runplot(scene_sdm, run_sdm, "shallow_demersals")
rsim.runplot(scene_fit, run.base, "shallow_demersals")
# lg_sculpins
rsim.runplot(scene_lsc, run_lsc, "lg_sculpins")
rsim.runplot(scene_fit, run.base, "lg_sculpins")
# octopus
rsim.runplot(scene_oct, run_oct, "octopus")
rsim.runplot(scene_fit, run.base, "octopus")
# mycto_bathy
rsim.runplot(scene_myb, run_myb, "mycto_bathy")
rsim.runplot(scene_fit, run.base, "mycto_bathy")
# capelin
rsim.runplot(scene_cap, run_cap, "capelin")
rsim.runplot(scene_fit, run.base, "capelin")
# sandlance
rsim.runplot(scene_snd, run_snd, "sandlance")
rsim.runplot(scene_fit, run.base, "sandlance")
# oth_forage
rsim.runplot(scene_ofg, run_ofg, "oth_forage")
rsim.runplot(scene_fit, run.base, "oth_forage")
# bairdi
rsim.runplot(scene_bdi, run_bdi, "bairdi")
rsim.runplot(scene_fit, run.base, "bairdi")
# king_crabs
rsim.runplot(scene_kng, run_kng, "king_crabs")
rsim.runplot(scene_fit, run.base, "king_crabs")
# opilio
rsim.runplot(scene_opi, run_opi, plot.species)
rsim.runplot(scene_opi, run_opi, "opilio")
rsim.runplot(scene_fit, run.base, "opilio")
# pandalid
rsim.runplot(scene_pnd, run_pnd, "pandalid")
rsim.runplot(scene_fit, run.base, "pandalid")
# snails_sea_stars
rsim.runplot(scene_mot, run_mot, "snails_sea_stars")
rsim.runplot(scene_fit, run.base, "snails_sea_stars")
# herb_echinoderms
rsim.runplot(scene_hec, run_hec, "herb_echinoderms")
rsim.runplot(scene_fit, run.base, "herb_echinoderms")
# misc_crabs
rsim.runplot(scene_mcr, run_mcr, "misc_crabs")
rsim.runplot(scene_fit, run.base, "misc_crabs")
# structural_epifauna
rsim.runplot(scene_sep, run_sep, "structural_epifauna")
rsim.runplot(scene_fit, run.base, "structural_epifauna")
# jellyfish
rsim.runplot(scene_jyf, run_jyf, "jellyfish")
rsim.runplot(scene_fit, run.base, "jellyfish")
# euphausiids
rsim.runplot(scene_eup, run_eup, "euphausiids")
rsim.runplot(scene_fit, run.base, "euphausiids")
# copepods
rsim.runplot(scene_cop, run_cop, "copepods")
rsim.runplot(scene_fit, run.base, "copepods")

# ---------------------------------------------------------------------------- #
# fitting vulnerabilities only for groups with biomass time series
ts_groups <- c("gray_whales", "nf_seal_adu", "oth_birds", "murres_puffins",
               "kittiwakes", "auklets", "fulmars", "sharks", "pollock_juv",
               "pollock_adu", "pcod_juv", "pcod_adu", "herring", 
               "arrowtooth_juv", "arrowtooth_adu", "kamchatka", "gr_turbot_juv",
               "gr_turbot_adu", "halibut_juv", "halibut_adu", "yf_sole", 
               "fh_sole", "nr_sole", "ak_plaice", "oth_flatfish", "ak_skate",
               "oth_skate", "sablefish_adu", "eelpouts", "deep_demersals",
               "pac_ocean_perch", "oth_rockfish", "north_rockfish", 
               "shortraker_rock", "atka", "shallow_demersals", "lg_sculpins",
               "octopus", "mycto_bathy", "sandlance", "oth_forage", "bairdi",
               "king_crabs", "opilio", "pandalid", "snails_sea_stars",
               "herb_echinoderms", "misc_crabs", "structural_epifauna",
               "jellyfish", "euphausiids", "copepods")
ts_groups_vartype <- c(rep("preyvul", length(ts_groups)), 
                     rep("predvul",length(ts_groups))) 
ts_groups_species <- c(ts_groups, ts_groups)
ts_groups_start   <- rep(0,length(ts_groups_species))

ptm <- proc.time()
ts_groups_opt <- optim(ts_groups_start, rsim.fit.run, method = "L-BFGS-B", 
                     lower = -4.5, upper = 4.5, hessian = FALSE, 
                     control = list(trace = 4, factr=1e13),
                     species=ts_groups_species, vartype=ts_groups_vartype,
                     scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
ts_groups_opt
#View(data.frame(ts_groups_species, ts_groups_vartype, ts_groups_start, ts_groups_opt$par))
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.ts_groups <- scene_fit
scene.ts_groups$params <- rsim.fit.apply(ts_groups_opt$par, ts_groups_species, 
                                       ts_groups_vartype, scene_fit$params)
# View(data.frame(scene.ts_groups$params$spname[scene.ts_groups$params$spname[scene.ts_groups$params$PreyFrom+1]],
#                 scene.ts_groups$params$spname[scene.ts_groups$params$spname[scene.ts_groups$params$PreyTo+1]],
#                 scene.ts_groups$params$VV))

# Test that we get the correct likelihood in the new scene
ts_groups_nll <- rsim.fit.run(NA, NA, NA, scene=scene.ts_groups, run_method="AB", 
                            run_years=hind_years, verbose=F); ts_groups_nll  
# But can also preserve parameters and run same thing with the original.    
ts_groups.result <- ts_groups_opt$par
rsim.fit.run(ts_groups.result, ts_groups_species, ts_groups_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.ts_groups <- rsim.fit.run(ts_groups.result, 
                              ts_groups_species, 
                              ts_groups_vartype, 
                              scene=scene_fit, 
                              run_method="AB", 
                              run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.ts_groups, plot.species)

ts_groups.like <- rsim.fit.table(scene_fit, run.ts_groups)  
# what changed?            
# View(ts_groups.like - base.like)
like.diff.ts_groups <- ts_groups.like - base.like
like.diff.ts_groups[order(like.diff.ts_groups$Biomass),]
ts_groups_nll - run.base.nll

cbind(ts_groups, 
      round(1 + exp(ts_groups_opt$par[1:52]),4), 
      round(1 + exp(ts_groups_opt$par[51:104]),4))

# ---------------------------------------------------------------------------- #
# fitting vulnerabilities for groups with biomass time series and lower TL
# groups without time series
ts_groups_ltl <- c("gray_whales", "nf_seal_adu", "oth_birds", "murres_puffins",
               "kittiwakes", "auklets", "fulmars", "sharks", "pollock_juv",
               "pollock_adu", "pcod_juv", "pcod_adu", "herring", 
               "arrowtooth_juv", "arrowtooth_adu", "kamchatka", "gr_turbot_juv",
               "gr_turbot_adu", "halibut_juv", "halibut_adu", "yf_sole", 
               "fh_sole", "nr_sole", "ak_plaice", "oth_flatfish", "ak_skate",
               "oth_skate", "sablefish_adu", "eelpouts", "deep_demersals",
               "pac_ocean_perch", "oth_rockfish", "north_rockfish", 
               "shortraker_rock", "atka", "shallow_demersals", "lg_sculpins",
               "octopus", "mycto_bathy", "sandlance", "oth_forage", "bairdi",
               "king_crabs", "opilio", "pandalid", "snails_sea_stars",
               "herb_echinoderms", "misc_crabs", "structural_epifauna",
               "jellyfish", "euphausiids", "copepods",
               # adding lower TL groups without time series here:
               "lg_phytoplankton", "sm_phytoplankton",
               "pelagic_microbes", "benthic_microbes",
               "pel_zooplankton", "ben_zooplankton",
               "infauna")
ts_groups_ltl_vartype <- c(rep("preyvul", length(ts_groups_ltl)), 
                           rep("predvul",length(ts_groups_ltl))) 
ts_groups_ltl_species <- c(ts_groups_ltl, ts_groups_ltl)
ts_groups_ltl_start   <- rep(0,length(ts_groups_ltl_species))

ptm <- proc.time()
ts_groups_ltl_opt <- optim(ts_groups_ltl_start, rsim.fit.run, method = "L-BFGS-B", 
                           lower = -4.5, upper = 4.5, hessian = FALSE, 
                           control = list(trace = 4, factr=1e13), 
                           species=ts_groups_ltl_species, 
                           vartype=ts_groups_ltl_vartype, scene=scene_fit, 
                           run_method="AB", run_years=hind_years)
proc.time() - ptm
ts_groups_ltl_opt
#View(data.frame(ts_groups_ltl_species, ts_groups_ltl_vartype, ts_groups_ltl_start, ts_groups_ltl_opt$par))
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.ts_groups_ltl <- scene_fit
scene.ts_groups_ltl$params <- rsim.fit.apply(ts_groups_ltl_opt$par, ts_groups_ltl_species, 
                                         ts_groups_ltl_vartype, scene_fit$params)
# View(data.frame(scene.ts_groups_ltl$params$spname[scene.ts_groups_ltl$params$spname[scene.ts_groups_ltl$params$PreyFrom+1]],
#                 scene.ts_groups_ltl$params$spname[scene.ts_groups_ltl$params$spname[scene.ts_groups_ltl$params$PreyTo+1]],
#                 scene.ts_groups_ltl$params$VV))

# Test that we get the correct likelihood in the new scene
ts_groups_ltl_nll <- rsim.fit.run(NA, NA, NA, scene=scene.ts_groups_ltl, run_method="AB", 
                              run_years=hind_years, verbose=F); ts_groups_ltl_nll  
# But can also preserve parameters and run same thing with the original.    
ts_groups_ltl.result <- ts_groups_ltl_opt$par
rsim.fit.run(ts_groups_ltl.result, ts_groups_ltl_species, ts_groups_ltl_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.ts_groups_ltl <- rsim.fit.run(ts_groups_ltl.result, ts_groups_ltl_species, ts_groups_ltl_vartype, 
                              scene=scene_fit, run_method="AB", 
                              run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.ts_groups_ltl, plot.species)

ts_groups_ltl.like <- rsim.fit.table(scene_fit, run.ts_groups_ltl)  
# what changed?            
# View(ts_groups_ltl.like - base.like)
like.diff.ts_groups_ltl <- ts_groups_ltl.like - base.like
like.diff.ts_groups_ltl[order(like.diff.ts_groups_ltl$Biomass),]
ts_groups_ltl_nll - run.base.nll

cbind(ts_groups_ltl, 
      round(1 + exp(ts_groups_ltl_opt$par[1:52]),4), 
      round(1 + exp(ts_groups_ltl_opt$par[51:104]),4))





################################################################################
# To do them manually one at a time and look at plots, etc. ------------------ #
# Create object to store the vulnerabilities
# ebs abbreviated name lookup
ebs_name_lookup <- read.csv("C:/Users/andy.whitehouse/Work/src/aclim3_rpath/a3_group_3letter_lookup.csv")
ebs_vv <- matrix(NA, ncol = 7, nrow = (bal$NUM_LIVING + bal$NUM_DEAD))
colnames(ebs_vv)  <- c("spnum","short_name","preyvul","predvul","nll","dnll","convrg")
row.names(ebs_vv) <- ebs_name_lookup$Group
ebs_vv[,1] <- scene_fit$params$spnum[2:77]
ebs_vv[,2] <- ebs_name_lookup$short_name
ebs_vv <- as.data.frame(ebs_vv)

# ---------------------------------------------------------------------------- #
ben_det <- c("benthic_detritus")
ben_det_vartype <- c(rep("preyvul", length(ben_det)), 
                     rep("predvul",length(ben_det))) 
ben_det_species <- c(ben_det, ben_det)
ben_det_start   <- rep(0,length(ben_det_species))

ben_det_opt <- optim(ben_det_start, rsim.fit.run, method = "L-BFGS-B", 
                       lower = -4.5, upper = 4.5, hessian = FALSE, 
                       control = list(trace = 5, factr=1e11),
                       species=ben_det_species, vartype=ben_det_vartype,
                       scene=scene_fit, run_method="AB", run_years=hind_years)
ben_det_opt
#View(data.frame(ben_det_species, ben_det_vartype, ben_det_start, ben_det_opt$par))
# store vulnerabilities
ebs_vv["benthic_detritus", 3:4] <- ben_det_opt$par
# store convergence
ebs_vv["benthic_detritus", 7] <- ben_det_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.ben_det <- scene_fit
scene.ben_det$params <- rsim.fit.apply(ben_det_opt$par, ben_det_species, 
                                       ben_det_vartype, scene_fit$params)
# View(data.frame(scene.ben_det$params$spname[scene.ben_det$params$spname[scene.ben_det$params$PreyFrom+1]],
#                 scene.ben_det$params$spname[scene.ben_det$params$spname[scene.ben_det$params$PreyTo+1]],
#                 scene.ben_det$params$VV))

# Test that we get the correct likelihood in the new scene
ben_det_nll <- rsim.fit.run(NA, NA, NA, scene=scene.ben_det, run_method="AB", 
                            run_years=hind_years, verbose=F); ben_det_nll  
ebs_vv["benthic_detritus", 5] <- ben_det_nll
# But can also preserve parameters and run same thing with the original.    
ben_det.result <- ben_det_opt$par
rsim.fit.run(ben_det.result, ben_det_species, ben_det_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.ben_det <- rsim.fit.run(ben_det.result, ben_det_species, ben_det_vartype, 
                            scene=scene_fit, run_method="AB", 
                            run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.ben_det, plot.species)

ben_det.like <- rsim.fit.table(scene_fit, run.ben_det)  
# what changed?            
# View(ben_det.like - base.like)
like.diff.ben_det <- ben_det.like - base.like
like.diff.ben_det[order(like.diff.ben_det$Biomass),]
ben_det_nll - run.base.nll
ebs_vv["benthic_detritus", 6] <- ben_det_nll - run.base.nll

# ---------------------------------------------------------------------------- #
pel_det <- c("pelagic_detritus")
pel_det_vartype <- c(rep("preyvul", length(pel_det)), 
                     rep("predvul",length(pel_det))) 
pel_det_species <- c(pel_det, pel_det)
pel_det_start   <- rep(0,length(pel_det_species))

pel_det_opt <- optim(pel_det_start, rsim.fit.run, method = "L-BFGS-B", 
                     lower = -4.5, upper = 4.5, hessian = FALSE, 
                     control = list(trace = 5, factr=1e11),
                     species=pel_det_species, vartype=pel_det_vartype,
                     scene=scene_fit, run_method="AB", run_years=hind_years)
pel_det_opt
#View(data.frame(pel_det_species, pel_det_vartype, pel_det_start, pel_det_opt$par))
# store vulnerabilities
ebs_vv["pelagic_detritus", 3:4] <- pel_det_opt$par
# store convergence
ebs_vv["pelagic_detritus", 7] <- pel_det_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.pel_det <- scene_fit
scene.pel_det$params <- rsim.fit.apply(pel_det_opt$par, pel_det_species, 
                                       pel_det_vartype, scene_fit$params)
# View(data.frame(scene.pel_det$params$spname[scene.pel_det$params$spname[scene.pel_det$params$PreyFrom+1]],
#                 scene.pel_det$params$spname[scene.pel_det$params$spname[scene.pel_det$params$PreyTo+1]],
#                 scene.pel_det$params$VV))

# Test that we get the correct likelihood in the new scene
pel_det_nll <- rsim.fit.run(NA, NA, NA, scene=scene.pel_det, run_method="AB", 
                            run_years=hind_years, verbose=F); pel_det_nll  
ebs_vv["pelagic_detritus", 5] <- pel_det_nll
# But can also preserve parameters and run same thing with the original.    
pel_det.result <- pel_det_opt$par
rsim.fit.run(pel_det.result, pel_det_species, pel_det_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.pel_det <- rsim.fit.run(pel_det.result, pel_det_species, pel_det_vartype, 
                            scene=scene_fit, run_method="AB", 
                            run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.pel_det, plot.species)

pel_det.like <- rsim.fit.table(scene_fit, run.pel_det)  
# what changed?            
# View(pel_det.like - base.like)
like.diff.pel_det <- pel_det.like - base.like
like.diff.pel_det[order(like.diff.pel_det$Biomass),]
pel_det_nll - run.base.nll
ebs_vv["pelagic_detritus", 6] <- pel_det_nll - run.base.nll


# ---------------------------------------------------------------------------- #
disc_off <- c("discards_offal")
disc_off_vartype <- c(rep("preyvul", length(disc_off)), 
                     rep("predvul",length(disc_off))) 
disc_off_species <- c(disc_off, disc_off)
disc_off_start   <- rep(0,length(disc_off_species))

disc_off_opt <- optim(disc_off_start, rsim.fit.run, method = "L-BFGS-B", 
                     lower = -4.5, upper = 4.5, hessian = FALSE, 
                     control = list(trace = 5, factr=1e11),
                     species=disc_off_species, vartype=disc_off_vartype,
                     scene=scene_fit, run_method="AB", run_years=hind_years)
disc_off_opt
#View(data.frame(disc_off_species, disc_off_vartype, disc_off_start, disc_off_opt$par))
# store vulnerabilities
ebs_vv["discards_offal", 3:4] <- disc_off_opt$par
# store convergence
ebs_vv["discards_offal", 7] <- disc_off_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.disc_off <- scene_fit
scene.disc_off$params <- rsim.fit.apply(disc_off_opt$par, disc_off_species, 
                                       disc_off_vartype, scene_fit$params)
# View(data.frame(scene.disc_off$params$spname[scene.disc_off$params$spname[scene.disc_off$params$PreyFrom+1]],
#                 scene.disc_off$params$spname[scene.disc_off$params$spname[scene.disc_off$params$PreyTo+1]],
#                 scene.disc_off$params$VV))

# Test that we get the correct likelihood in the new scene
disc_off_nll <- rsim.fit.run(NA, NA, NA, scene=scene.disc_off, run_method="AB", 
                            run_years=hind_years, verbose=F); disc_off_nll  
ebs_vv["discards_offal", 5] <- disc_off_nll
# But can also preserve parameters and run same thing with the original.    
disc_off.result <- disc_off_opt$par
rsim.fit.run(disc_off.result, disc_off_species, disc_off_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.disc_off <- rsim.fit.run(disc_off.result, disc_off_species, disc_off_vartype, 
                            scene=scene_fit, run_method="AB", 
                            run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.disc_off, plot.species)

disc_off.like <- rsim.fit.table(scene_fit, run.disc_off)  
# what changed?            
# View(disc_off.like - base.like)
like.diff.disc_off <- disc_off.like - base.like
like.diff.disc_off[order(like.diff.disc_off$Biomass),]
disc_off_nll - run.base.nll
ebs_vv["discards_offal", 6] <- disc_off_nll - run.base.nll


# ---------------------------------------------------------------------------- #
sm_pp <- c("sm_phytoplankton")
sm_pp_vartype <- c(rep("preyvul", length(sm_pp)), 
                      rep("predvul",length(sm_pp))) 
sm_pp_species <- c(sm_pp, sm_pp)
sm_pp_start   <- rep(0,length(sm_pp_species))

sm_pp_opt <- optim(sm_pp_start, rsim.fit.run, method = "L-BFGS-B", 
                      lower = -4.5, upper = 4.5, hessian = FALSE, 
                      control = list(trace = 5, factr=1e11),
                      species=sm_pp_species, vartype=sm_pp_vartype,
                      scene=scene_fit, run_method="AB", run_years=hind_years)
sm_pp_opt
#View(data.frame(sm_pp_species, sm_pp_vartype, sm_pp_start, sm_pp_opt$par))
# store vulnerabilities
ebs_vv["sm_phytoplankton", 3:4] <- sm_pp_opt$par
# store convergence
ebs_vv["sm_phytoplankton", 7] <- sm_pp_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.sm_pp <- scene_fit
scene.sm_pp$params <- rsim.fit.apply(sm_pp_opt$par, sm_pp_species, 
                                        sm_pp_vartype, scene_fit$params)
# View(data.frame(scene.sm_pp$params$spname[scene.sm_pp$params$spname[scene.sm_pp$params$PreyFrom+1]],
#                 scene.sm_pp$params$spname[scene.sm_pp$params$spname[scene.sm_pp$params$PreyTo+1]],
#                 scene.sm_pp$params$VV))

# Test that we get the correct likelihood in the new scene
sm_pp_nll <- rsim.fit.run(NA, NA, NA, scene=scene.sm_pp, run_method="AB", 
                             run_years=hind_years, verbose=F); sm_pp_nll  
ebs_vv["sm_phytoplankton", 5] <- sm_pp_nll
# But can also preserve parameters and run same thing with the original.    
sm_pp.result <- sm_pp_opt$par
rsim.fit.run(sm_pp.result, sm_pp_species, sm_pp_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.sm_pp <- rsim.fit.run(sm_pp.result, sm_pp_species, sm_pp_vartype, 
                             scene=scene_fit, run_method="AB", 
                             run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.sm_pp, plot.species)

sm_pp.like <- rsim.fit.table(scene_fit, run.sm_pp)  
# what changed?            
# View(sm_pp.like - base.like)
like.diff.sm_pp <- sm_pp.like - base.like
like.diff.sm_pp[order(like.diff.sm_pp$Biomass),]
sm_pp_nll - run.base.nll
ebs_vv["sm_phytoplankton", 6] <- sm_pp_nll - run.base.nll


# ---------------------------------------------------------------------------- #
lg_pp <- c("lg_phytoplankton")
lg_pp_vartype <- c(rep("preyvul", length(lg_pp)), 
                   rep("predvul",length(lg_pp))) 
lg_pp_species <- c(lg_pp, lg_pp)
lg_pp_start   <- rep(0,length(lg_pp_species))

lg_pp_opt <- optim(lg_pp_start, rsim.fit.run, method = "L-BFGS-B", 
                   lower = -4.5, upper = 4.5, hessian = FALSE, 
                   control = list(trace = 5, factr=1e11),
                   species=lg_pp_species, vartype=lg_pp_vartype,
                   scene=scene_fit, run_method="AB", run_years=hind_years)
lg_pp_opt
#View(data.frame(lg_pp_species, lg_pp_vartype, lg_pp_start, lg_pp_opt$par))
# store vulnerabilities
ebs_vv["lg_phytoplankton", 3:4] <- lg_pp_opt$par
# store convergence
ebs_vv["lg_phytoplankton", 7] <- lg_pp_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.lg_pp <- scene_fit
scene.lg_pp$params <- rsim.fit.apply(lg_pp_opt$par, lg_pp_species, 
                                     lg_pp_vartype, scene_fit$params)
# View(data.frame(scene.lg_pp$params$spname[scene.lg_pp$params$spname[scene.lg_pp$params$PreyFrom+1]],
#                 scene.lg_pp$params$spname[scene.lg_pp$params$spname[scene.lg_pp$params$PreyTo+1]],
#                 scene.lg_pp$params$VV))

# Test that we get the correct likelihood in the new scene
lg_pp_nll <- rsim.fit.run(NA, NA, NA, scene=scene.lg_pp, run_method="AB", 
                          run_years=hind_years, verbose=F); lg_pp_nll  
ebs_vv["lg_phytoplankton", 5] <- lg_pp_nll
# But can also preserve parameters and run same thing with the original.    
lg_pp.result <- lg_pp_opt$par
rsim.fit.run(lg_pp.result, lg_pp_species, lg_pp_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.lg_pp <- rsim.fit.run(lg_pp.result, lg_pp_species, lg_pp_vartype, 
                          scene=scene_fit, run_method="AB", 
                          run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.lg_pp, plot.species)

lg_pp.like <- rsim.fit.table(scene_fit, run.lg_pp)  
# what changed?            
# View(lg_pp.like - base.like)
like.diff.lg_pp <- lg_pp.like - base.like
like.diff.lg_pp[order(like.diff.lg_pp$Biomass),]
lg_pp_nll - run.base.nll
ebs_vv["lg_phytoplankton", 6] <- lg_pp_nll - run.base.nll


# ---------------------------------------------------------------------------- #
ben_micr <- c("benthic_microbes")
ben_micr_vartype <- c(rep("preyvul", length(ben_micr)), 
                   rep("predvul",length(ben_micr))) 
ben_micr_species <- c(ben_micr, ben_micr)
ben_micr_start   <- rep(0,length(ben_micr_species))

ben_micr_opt <- optim(ben_micr_start, rsim.fit.run, method = "L-BFGS-B", 
                   lower = -4.5, upper = 4.5, hessian = FALSE, 
                   control = list(trace = 5, factr=1e11),
                   species=ben_micr_species, vartype=ben_micr_vartype,
                   scene=scene_fit, run_method="AB", run_years=hind_years)
ben_micr_opt
#View(data.frame(ben_micr_species, ben_micr_vartype, ben_micr_start, ben_micr_opt$par))
# store vulnerabilities
ebs_vv["benthic_microbes", 3:4] <- ben_micr_opt$par
# store convergence
ebs_vv["benthic_microbes", 7] <- ben_micr_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.ben_micr <- scene_fit
scene.ben_micr$params <- rsim.fit.apply(ben_micr_opt$par, ben_micr_species, 
                                     ben_micr_vartype, scene_fit$params)
# View(data.frame(scene.ben_micr$params$spname[scene.ben_micr$params$spname[scene.ben_micr$params$PreyFrom+1]],
#                 scene.ben_micr$params$spname[scene.ben_micr$params$spname[scene.ben_micr$params$PreyTo+1]],
#                 scene.ben_micr$params$VV))

# Test that we get the correct likelihood in the new scene
ben_micr_nll <- rsim.fit.run(NA, NA, NA, scene=scene.ben_micr, run_method="AB", 
                          run_years=hind_years, verbose=F); ben_micr_nll  
ebs_vv["benthic_microbes", 5] <- ben_micr_nll
# But can also preserve parameters and run same thing with the original.    
ben_micr.result <- ben_micr_opt$par
rsim.fit.run(ben_micr.result, ben_micr_species, ben_micr_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.ben_micr <- rsim.fit.run(ben_micr.result, ben_micr_species, ben_micr_vartype, 
                          scene=scene_fit, run_method="AB", 
                          run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.ben_micr, plot.species)

ben_micr.like <- rsim.fit.table(scene_fit, run.ben_micr)  
# what changed?            
# View(ben_micr.like - base.like)
like.diff.ben_micr <- ben_micr.like - base.like
like.diff.ben_micr[order(like.diff.ben_micr$Biomass),]
ben_micr_nll - run.base.nll
ebs_vv["benthic_microbes", 6] <- ben_micr_nll - run.base.nll


# ---------------------------------------------------------------------------- #
pel_micr <- c("pelagic_microbes")
pel_micr_vartype <- c(rep("preyvul", length(pel_micr)), 
                   rep("predvul",length(pel_micr))) 
pel_micr_species <- c(pel_micr, pel_micr)
pel_micr_start   <- rep(0,length(pel_micr_species))

pel_micr_opt <- optim(pel_micr_start, rsim.fit.run, method = "L-BFGS-B", 
                   lower = -4.5, upper = 4.5, hessian = FALSE, 
                   control = list(trace = 5, factr=1e11),
                   species=pel_micr_species, vartype=pel_micr_vartype,
                   scene=scene_fit, run_method="AB", run_years=hind_years)
pel_micr_opt
#View(data.frame(pel_micr_species, pel_micr_vartype, pel_micr_start, pel_micr_opt$par))
# store vulnerabilities
ebs_vv["pelagic_microbes", 3:4] <- pel_micr_opt$par
# store convergence
ebs_vv["pelagic_microbes", 7] <- pel_micr_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.pel_micr <- scene_fit
scene.pel_micr$params <- rsim.fit.apply(pel_micr_opt$par, pel_micr_species, 
                                     pel_micr_vartype, scene_fit$params)
# View(data.frame(scene.pel_micr$params$spname[scene.pel_micr$params$spname[scene.pel_micr$params$PreyFrom+1]],
#                 scene.pel_micr$params$spname[scene.pel_micr$params$spname[scene.pel_micr$params$PreyTo+1]],
#                 scene.pel_micr$params$VV))

# Test that we get the correct likelihood in the new scene
pel_micr_nll <- rsim.fit.run(NA, NA, NA, scene=scene.pel_micr, run_method="AB", 
                          run_years=hind_years, verbose=F); pel_micr_nll  
ebs_vv["pelagic_microbes", 5] <- pel_micr_nll
# But can also preserve parameters and run same thing with the original.    
pel_micr.result <- pel_micr_opt$par
rsim.fit.run(pel_micr.result, pel_micr_species, pel_micr_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.pel_micr <- rsim.fit.run(pel_micr.result, pel_micr_species, pel_micr_vartype, 
                          scene=scene_fit, run_method="AB", 
                          run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.pel_micr, plot.species)

pel_micr.like <- rsim.fit.table(scene_fit, run.pel_micr)  
# what changed?            
# View(pel_micr.like - base.like)
like.diff.pel_micr <- pel_micr.like - base.like
like.diff.pel_micr[order(like.diff.pel_micr$Biomass),]
pel_micr_nll - run.base.nll
ebs_vv["pelagic_microbes", 6] <- pel_micr_nll - run.base.nll


# ---------------------------------------------------------------------------- #
cop <- c("copepods")
cop_vartype <- c(rep("preyvul", length(cop)), 
                      rep("predvul",length(cop))) 
cop_species <- c(cop, cop)
cop_start   <- rep(0,length(cop_species))

cop_opt <- optim(cop_start, rsim.fit.run, method = "L-BFGS-B", 
                      lower = -4.5, upper = 4.5, hessian = FALSE, 
                      control = list(trace = 5, factr=1e11),
                      species=cop_species, vartype=cop_vartype,
                      scene=scene_fit, run_method="AB", run_years=hind_years)
cop_opt
#View(data.frame(cop_species, cop_vartype, cop_start, cop_opt$par))
# store vulnerabilities
ebs_vv["copepods", 3:4] <- cop_opt$par
# store convergence
ebs_vv["copepods", 7] <- cop_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.cop <- scene_fit
scene.cop$params <- rsim.fit.apply(cop_opt$par, cop_species, 
                                        cop_vartype, scene_fit$params)
# View(data.frame(scene.cop$params$spname[scene.cop$params$spname[scene.cop$params$PreyFrom+1]],
#                 scene.cop$params$spname[scene.cop$params$spname[scene.cop$params$PreyTo+1]],
#                 scene.cop$params$VV))

# Test that we get the correct likelihood in the new scene
cop_nll <- rsim.fit.run(NA, NA, NA, scene=scene.cop, run_method="AB", 
                             run_years=hind_years, verbose=F); cop_nll  
ebs_vv["copepods", 5] <- cop_nll
# But can also preserve parameters and run same thing with the original.    
cop.result <- cop_opt$par
rsim.fit.run(cop.result, cop_species, cop_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.cop <- rsim.fit.run(cop.result, cop_species, cop_vartype, 
                             scene=scene_fit, run_method="AB", 
                             run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.cop, plot.species)

cop.like <- rsim.fit.table(scene_fit, run.cop)  
# what changed?            
# View(cop.like - base.like)
like.diff.cop <- cop.like - base.like
like.diff.cop[order(like.diff.cop$Biomass),]
cop_nll - run.base.nll
ebs_vv["copepods", 6] <- cop_nll - run.base.nll


# ---------------------------------------------------------------------------- #
eup <- c("euphausiids")
eup_vartype <- c(rep("preyvul", length(eup)), 
                 rep("predvul",length(eup))) 
eup_species <- c(eup, eup)
eup_start   <- rep(0,length(eup_species))

eup_opt <- optim(eup_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=eup_species, vartype=eup_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
eup_opt
#View(data.frame(eup_species, eup_vartype, eup_start, eup_opt$par))
# store vulnerabilities
ebs_vv["euphausiids", 3:4] <- eup_opt$par
# store convergence
ebs_vv["euphausiids", 7] <- eup_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.eup <- scene_fit
scene.eup$params <- rsim.fit.apply(eup_opt$par, eup_species, 
                                   eup_vartype, scene_fit$params)
# View(data.frame(scene.eup$params$spname[scene.eup$params$spname[scene.eup$params$PreyFrom+1]],
#                 scene.eup$params$spname[scene.eup$params$spname[scene.eup$params$PreyTo+1]],
#                 scene.eup$params$VV))

# Test that we get the correct likelihood in the new scene
eup_nll <- rsim.fit.run(NA, NA, NA, scene=scene.eup, run_method="AB", 
                        run_years=hind_years, verbose=F); eup_nll  
ebs_vv["euphausiids", 5] <- eup_nll
# But can also preserve parameters and run same thing with the original.    
eup.result <- eup_opt$par
rsim.fit.run(eup.result, eup_species, eup_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.eup <- rsim.fit.run(eup.result, eup_species, eup_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.eup, plot.species)

eup.like <- rsim.fit.table(scene_fit, run.eup)  
# what changed?            
# View(eup.like - base.like)
like.diff.eup <- eup.like - base.like
like.diff.eup[order(like.diff.eup$Biomass),]
eup_nll - run.base.nll
ebs_vv["euphausiids", 6] <- eup_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.eup$Biomass < -9.499)
row.names(like.diff.eup)[down_nll]
print(cbind(row.names(like.diff.eup)[down_nll], like.diff.eup[down_nll,1]))

up_nll <- which(like.diff.eup$Biomass > 9.499)
row.names(like.diff.eup)[up_nll]
print(cbind(row.names(like.diff.eup)[up_nll], like.diff.eup[up_nll,1]))


# ---------------------------------------------------------------------------- #
inf <- c("infauna")
inf_vartype <- c(rep("preyvul", length(inf)), 
                 rep("predvul",length(inf))) 
inf_species <- c(inf, inf)
inf_start   <- rep(0,length(inf_species))

inf_opt <- optim(inf_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=inf_species, vartype=inf_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
inf_opt
#View(data.frame(inf_species, inf_vartype, inf_start, inf_opt$par))
# store vulnerabilities
ebs_vv["infauna", 3:4] <- inf_opt$par
# store convergence
ebs_vv["infauna", 7] <- inf_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.inf <- scene_fit
scene.inf$params <- rsim.fit.apply(inf_opt$par, inf_species, 
                                   inf_vartype, scene_fit$params)
# View(data.frame(scene.inf$params$spname[scene.inf$params$spname[scene.inf$params$PreyFrom+1]],
#                 scene.inf$params$spname[scene.inf$params$spname[scene.inf$params$PreyTo+1]],
#                 scene.inf$params$VV))

# Test that we get the correct likelihood in the new scene
inf_nll <- rsim.fit.run(NA, NA, NA, scene=scene.inf, run_method="AB", 
                        run_years=hind_years, verbose=F); inf_nll  
ebs_vv["infauna", 5] <- inf_nll
# But can also preserve parameters and run same thing with the original.    
inf.result <- inf_opt$par
rsim.fit.run(inf.result, inf_species, inf_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.inf <- rsim.fit.run(inf.result, inf_species, inf_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.inf, plot.species)

inf.like <- rsim.fit.table(scene_fit, run.inf)  
# what changed?            
# View(inf.like - base.like)
like.diff.inf <- inf.like - base.like
like.diff.inf[order(like.diff.inf$Biomass),]
inf_nll - run.base.nll
ebs_vv["infauna", 6] <- inf_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.inf$Biomass < -9.499)
row.names(like.diff.inf)[down_nll]
print(cbind(row.names(like.diff.inf)[down_nll], like.diff.inf[down_nll,1]))

up_nll <- which(like.diff.inf$Biomass > 9.499)
row.names(like.diff.inf)[up_nll]
print(cbind(row.names(like.diff.inf)[up_nll], like.diff.inf[up_nll,1]))


# ---------------------------------------------------------------------------- #
bdi <- c("bairdi")
bdi_vartype <- c(rep("preyvul", length(bdi)), 
                 rep("predvul",length(bdi))) 
bdi_species <- c(bdi, bdi)
bdi_start   <- rep(0,length(bdi_species))

bdi_opt <- optim(bdi_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=bdi_species, vartype=bdi_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
bdi_opt
#View(data.frame(bdi_species, bdi_vartype, bdi_start, bdi_opt$par))
# store vulnerabilities
ebs_vv["bairdi", 3:4] <- bdi_opt$par
# store convergence
ebs_vv["bairdi", 7] <- bdi_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.bdi <- scene_fit
scene.bdi$params <- rsim.fit.apply(bdi_opt$par, bdi_species, 
                                   bdi_vartype, scene_fit$params)
View(data.frame(scene.bdi$params$spname[scene.bdi$params$spname[scene.bdi$params$PreyFrom+1]],
                scene.bdi$params$spname[scene.bdi$params$spname[scene.bdi$params$PreyTo+1]],
                scene.bdi$params$VV))

# Test that we get the correct likelihood in the new scene
bdi_nll <- rsim.fit.run(NA, NA, NA, scene=scene.bdi, run_method="AB", 
                        run_years=hind_years, verbose=F); bdi_nll  
ebs_vv["bairdi", 5] <- bdi_nll
# But can also preserve parameters and run same thing with the original.    
bdi.result <- bdi_opt$par
rsim.fit.run(bdi.result, bdi_species, bdi_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.bdi <- rsim.fit.run(bdi.result, bdi_species, bdi_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.bdi, plot.species)

bdi.like <- rsim.fit.table(scene_fit, run.bdi)  
# what changed?            
# View(bdi.like - base.like)
like.diff.bdi <- bdi.like - base.like
like.diff.bdi[order(like.diff.bdi$Biomass),]
bdi_nll - run.base.nll
ebs_vv["bairdi", 6] <- bdi_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.bdi$Biomass < -9.499)
row.names(like.diff.bdi)[down_nll]
print(cbind(row.names(like.diff.bdi)[down_nll], like.diff.bdi[down_nll,1]))

up_nll <- which(like.diff.bdi$Biomass > 9.499)
row.names(like.diff.bdi)[up_nll]
print(cbind(row.names(like.diff.bdi)[up_nll], like.diff.bdi[up_nll,1]))


# ---------------------------------------------------------------------------- #
ofg <- c("oth_forage")
ofg_vartype <- c(rep("preyvul", length(ofg)), 
                 rep("predvul",length(ofg))) 
ofg_species <- c(ofg, ofg)
ofg_start   <- rep(0,length(ofg_species))

ofg_opt <- optim(ofg_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=ofg_species, vartype=ofg_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
ofg_opt
#View(data.frame(ofg_species, ofg_vartype, ofg_start, ofg_opt$par))
# store vulnerabilities
ebs_vv["oth_forage", 3:4] <- ofg_opt$par
# store convergence
ebs_vv["oth_forage", 7] <- ofg_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.ofg <- scene_fit
scene.ofg$params <- rsim.fit.apply(ofg_opt$par, ofg_species, 
                                   ofg_vartype, scene_fit$params)
# View(data.frame(scene.ofg$params$spname[scene.ofg$params$spname[scene.ofg$params$PreyFrom+1]],
#                 scene.ofg$params$spname[scene.ofg$params$spname[scene.ofg$params$PreyTo+1]],
#                 scene.ofg$params$VV))

# Test that we get the correct likelihood in the new scene
ofg_nll <- rsim.fit.run(NA, NA, NA, scene=scene.ofg, run_method="AB", 
                        run_years=hind_years, verbose=F); ofg_nll  
ebs_vv["oth_forage", 5] <- ofg_nll
# But can also preserve parameters and run same thing with the original.    
ofg.result <- ofg_opt$par
rsim.fit.run(ofg.result, ofg_species, ofg_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.ofg <- rsim.fit.run(ofg.result, ofg_species, ofg_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.ofg, plot.species)

ofg.like <- rsim.fit.table(scene_fit, run.ofg)  
# what changed?            
# View(ofg.like - base.like)
like.diff.ofg <- ofg.like - base.like
like.diff.ofg[order(like.diff.ofg$Biomass),]
ofg_nll - run.base.nll
ebs_vv["oth_forage", 6] <- ofg_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.ofg$Biomass < -9.499)
row.names(like.diff.ofg)[down_nll]
print(cbind(row.names(like.diff.ofg)[down_nll], like.diff.ofg[down_nll,1]))

up_nll <- which(like.diff.ofg$Biomass > 9.499)
row.names(like.diff.ofg)[up_nll]
print(cbind(row.names(like.diff.ofg)[up_nll], like.diff.ofg[up_nll,1]))


# ---------------------------------------------------------------------------- #
bzp <- c("ben_zooplankton")
bzp_vartype <- c(rep("preyvul", length(bzp)), 
                 rep("predvul",length(bzp))) 
bzp_species <- c(bzp, bzp)
bzp_start   <- rep(0,length(bzp_species))

bzp_opt <- optim(bzp_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=bzp_species, vartype=bzp_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
bzp_opt
#View(data.frame(bzp_species, bzp_vartype, bzp_start, bzp_opt$par))
# store vulnerabilities
ebs_vv["ben_zooplankton", 3:4] <- bzp_opt$par
# store convergence
ebs_vv["ben_zooplankton", 7] <- bzp_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.bzp <- scene_fit
scene.bzp$params <- rsim.fit.apply(bzp_opt$par, bzp_species, 
                                   bzp_vartype, scene_fit$params)
# View(data.frame(scene.bzp$params$spname[scene.bzp$params$spname[scene.bzp$params$PreyFrom+1]],
#                 scene.bzp$params$spname[scene.bzp$params$spname[scene.bzp$params$PreyTo+1]],
#                 scene.bzp$params$VV))

# Test that we get the correct likelihood in the new scene
bzp_nll <- rsim.fit.run(NA, NA, NA, scene=scene.bzp, run_method="AB", 
                        run_years=hind_years, verbose=F); bzp_nll  
ebs_vv["ben_zooplankton", 5] <- bzp_nll
# But can also preserve parameters and run same thing with the original.    
bzp.result <- bzp_opt$par
rsim.fit.run(bzp.result, bzp_species, bzp_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.bzp <- rsim.fit.run(bzp.result, bzp_species, bzp_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.bzp, plot.species)

bzp.like <- rsim.fit.table(scene_fit, run.bzp)  
# what changed?            
# View(bzp.like - base.like)
like.diff.bzp <- bzp.like - base.like
like.diff.bzp[order(like.diff.bzp$Biomass),]
bzp_nll - run.base.nll
ebs_vv["ben_zooplankton", 6] <- bzp_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.bzp$Biomass < -9.499)
row.names(like.diff.bzp)[down_nll]
print(cbind(row.names(like.diff.bzp)[down_nll], like.diff.bzp[down_nll,1]))

up_nll <- which(like.diff.bzp$Biomass > 9.499)
row.names(like.diff.bzp)[up_nll]
print(cbind(row.names(like.diff.bzp)[up_nll], like.diff.bzp[up_nll,1]))


# ---------------------------------------------------------------------------- #
nrs <- c("nr_sole")
nrs_vartype <- c(rep("preyvul", length(nrs)), 
                 rep("predvul",length(nrs))) 
nrs_species <- c(nrs, nrs)
nrs_start   <- rep(0,length(nrs_species))

nrs_opt <- optim(nrs_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=nrs_species, vartype=nrs_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
nrs_opt
#View(data.frame(nrs_species, nrs_vartype, nrs_start, nrs_opt$par))
# store vulnerabilities
ebs_vv["nr_sole", 3:4] <- nrs_opt$par
# store convergence
ebs_vv["nr_sole", 7] <- nrs_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.nrs <- scene_fit
scene.nrs$params <- rsim.fit.apply(nrs_opt$par, nrs_species, 
                                   nrs_vartype, scene_fit$params)
# View(data.frame(scene.nrs$params$spname[scene.nrs$params$spname[scene.nrs$params$PreyFrom+1]],
#                 scene.nrs$params$spname[scene.nrs$params$spname[scene.nrs$params$PreyTo+1]],
#                 scene.nrs$params$VV))

# Test that we get the correct likelihood in the new scene
nrs_nll <- rsim.fit.run(NA, NA, NA, scene=scene.nrs, run_method="AB", 
                        run_years=hind_years, verbose=F); nrs_nll  
ebs_vv["nr_sole", 5] <- nrs_nll
# But can also preserve parameters and run same thing with the original.    
nrs.result <- nrs_opt$par
rsim.fit.run(nrs.result, nrs_species, nrs_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.nrs <- rsim.fit.run(nrs.result, nrs_species, nrs_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.nrs, plot.species)

nrs.like <- rsim.fit.table(scene_fit, run.nrs)  
# what changed?            
# View(nrs.like - base.like)
like.diff.nrs <- nrs.like - base.like
like.diff.nrs[order(like.diff.nrs$Biomass),]
nrs_nll - run.base.nll
ebs_vv["nr_sole", 6] <- nrs_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.nrs$Biomass < -9.499)
row.names(like.diff.nrs)[down_nll]
print(cbind(row.names(like.diff.nrs)[down_nll], like.diff.nrs[down_nll,1]))

up_nll <- which(like.diff.nrs$Biomass > 9.499)
row.names(like.diff.nrs)[up_nll]
print(cbind(row.names(like.diff.nrs)[up_nll], like.diff.nrs[up_nll,1]))


# ---------------------------------------------------------------------------- #
kam <- c("kamchatka")
kam_vartype <- c(rep("preyvul", length(kam)), 
                 rep("predvul",length(kam))) 
kam_species <- c(kam, kam)
kam_start   <- rep(0,length(kam_species))

kam_opt <- optim(kam_start, rsim.fit.run, method = "L-BFGS-B", 
                 lower = -4.5, upper = 4.5, hessian = FALSE, 
                 control = list(trace = 5, factr=1e11),
                 species=kam_species, vartype=kam_vartype,
                 scene=scene_fit, run_method="AB", run_years=hind_years)
kam_opt
#View(data.frame(kam_species, kam_vartype, kam_start, kam_opt$par))
# store vulnerabilities
ebs_vv["kamchatka", 3:4] <- kam_opt$par
# store convergence
ebs_vv["kamchatka", 7] <- kam_opt$convergence == 0
# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.kam <- scene_fit
scene.kam$params <- rsim.fit.apply(kam_opt$par, kam_species, 
                                   kam_vartype, scene_fit$params)
# View(data.frame(scene.kam$params$spname[scene.kam$params$spname[scene.kam$params$PreyFrom+1]],
#                 scene.kam$params$spname[scene.kam$params$spname[scene.kam$params$PreyTo+1]],
#                 scene.kam$params$VV))

# Test that we get the correct likelihood in the new scene
kam_nll <- rsim.fit.run(NA, NA, NA, scene=scene.kam, run_method="AB", 
                        run_years=hind_years, verbose=F); kam_nll  
ebs_vv["kamchatka", 5] <- kam_nll
# But can also preserve parameters and run same thing with the original.    
kam.result <- kam_opt$par
rsim.fit.run(kam.result, kam_species, kam_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 
# Save the actual full run (same as above with verbose=T)
run.kam <- rsim.fit.run(kam.result, kam_species, kam_vartype, 
                        scene=scene_fit, run_method="AB", 
                        run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.kam, plot.species)

kam.like <- rsim.fit.table(scene_fit, run.kam)  
# what changed?            
# View(kam.like - base.like)
like.diff.kam <- kam.like - base.like
like.diff.kam[order(like.diff.kam$Biomass),]
kam_nll - run.base.nll
ebs_vv["kamchatka", 6] <- kam_nll - run.base.nll

# which groups fit are most affected?
down_nll <- which(like.diff.kam$Biomass < -9.499)
row.names(like.diff.kam)[down_nll]
print(cbind(row.names(like.diff.kam)[down_nll], like.diff.kam[down_nll,1]))

up_nll <- which(like.diff.kam$Biomass > 9.499)
row.names(like.diff.kam)[up_nll]
print(cbind(row.names(like.diff.kam)[up_nll], like.diff.kam[up_nll,1]))











##### 2.Fitting vulnerabilities of opilio ####################################
opilio <- c("opilio")
opilio_vartype <- c(rep("preyvul", length(opilio)), rep("predvul",length(opilio))) 
opilio_species <- c(opilio, opilio)
opilio_start   <- rep(0,length(opilio_species))

ptm <- proc.time()
opilio_opt_AB <- optim(opilio_start, rsim.fit.run, method = "L-BFGS-B", 
                       lower = -4.5, upper = 4.5, hessian = FALSE, 
                       control = list(trace = 5, factr=1e13),
                       species=opilio_species, vartype=opilio_vartype,
                       scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
opilio_opt_AB

# Result: new vector for values  
opilio_opt_AB$par
#View(data.frame(opilio_species, opilio_vartype, opilio_start, opilio_opt_AB$par))

# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.update <- scene_fit
scene.update$params <- rsim.fit.apply(opilio_opt_AB$par, opilio_species, opilio_vartype, scene_fit$params)

View(data.frame(scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyFrom+1]],
                scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyTo+1]],
                scene.update$params$VV))

# Test that we get the correct likelihood in the new scene
rsim.fit.run(NA, NA, NA, 
             scene=scene.update, run_method="AB", run_years=hind_years, verbose=F)  
# But can also preserve parameters and run same thing with the original.    
opilio.result <- opilio_opt_AB$par
rsim.fit.run(opilio.result, opilio_species, opilio_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 

# Save the actual full run (same as above with verbose=T)
run.opilio <- rsim.fit.run(opilio.result, opilio_species, opilio_vartype, 
                           scene=scene_fit, run_method="AB", run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.opilio, plot.species)

opilio.like <- rsim.fit.table(scene_fit, run.opilio)  

# what changed?            
View(opilio.like - base.like)

##### 2.Fitting vulnerabilities of opilio ####################################
crabs <- c("bairdi", "king_crabs", "opilio")
crabs_vartype <- c(rep("preyvul", length(crabs)), rep("predvul",length(crabs))) 
crabs_species <- c(crabs, crabs)
crabs_start   <- rep(0,length(crabs_species))

ptm <- proc.time()
crabs_opt_AB <- optim(crabs_start, rsim.fit.run, method = "L-BFGS-B", 
                      lower = -4.5, upper = 4.5, hessian = FALSE, 
                      control = list(trace = 5, factr=1e14),
                      species=crabs_species, vartype=crabs_vartype,
                      scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
crabs_opt_AB

# Result: new vector for values  
crabs_opt_AB$par
#View(data.frame(crabs_species, crabs_vartype, crabs_start, crabs_opt_AB$par))

# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.update <- scene_fit
scene.update$params <- rsim.fit.apply(crabs_opt_AB$par, crabs_species, crabs_vartype, scene_fit$params)

View(data.frame(scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyFrom+1]],
                scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyTo+1]],
                scene.update$params$VV))

# Test that we get the correct likelihood in the new scene
rsim.fit.run(NA, NA, NA, 
             scene=scene.update, run_method="AB", run_years=hind_years, verbose=F)  
# But can also preserve parameters and run same thing with the original.    
crabs.result <- crabs_opt_AB$par
rsim.fit.run(crabs.result, crabs_species, crabs_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 

# Save the actual full run (same as above with verbose=T)
run.crabs <- rsim.fit.run(crabs.result, crabs_species, crabs_vartype, 
                          scene=scene_fit, run_method="AB", run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.crabs, plot.species)

crabs.like <- rsim.fit.table(scene_fit, run.crabs)  

# what changed?            
View(crabs.like - base.like)



##### 2.Fitting vulnerabilities of shlw_dem ####################################
shlw_dem <- c("shallow_demersals")
shlw_dem_vartype <- c(rep("preyvul", length(shlw_dem)), rep("predvul",length(shlw_dem))) 
shlw_dem_species <- c(shlw_dem, shlw_dem)
shlw_dem_start   <- rep(0,length(shlw_dem_species))

ptm <- proc.time()
shlw_dem_opt_AB <- optim(shlw_dem_start, rsim.fit.run, method = "L-BFGS-B", 
                         lower = -4.5, upper = 4.5, hessian = FALSE, 
                         control = list(trace = 5, factr=1e11),
                         species=shlw_dem_species, vartype=shlw_dem_vartype,
                         scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
shlw_dem_opt_AB

# Result: new vector for values  
shlw_dem_opt_AB$par
#View(data.frame(shlw_dem_species, shlw_dem_vartype, shlw_dem_start, shlw_dem_opt_AB$par))

# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.update <- scene_fit
scene.update$params <- rsim.fit.apply(shlw_dem_opt_AB$par, shlw_dem_species, shlw_dem_vartype, scene_fit$params)

View(data.frame(scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyFrom+1]],
                scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyTo+1]],
                scene.update$params$VV))

# Test that we get the correct likelihood in the new scene
rsim.fit.run(NA, NA, NA, 
             scene=scene.update, run_method="AB", run_years=hind_years, verbose=F)  
# But can also preserve parameters and run same thing with the original.    
shlw_dem.result <- shlw_dem_opt_AB$par
rsim.fit.run(shlw_dem.result, shlw_dem_species, shlw_dem_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 

# Save the actual full run (same as above with verbose=T)
run.shlw_dem <- rsim.fit.run(shlw_dem.result, shlw_dem_species, shlw_dem_vartype, 
                             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.shlw_dem, plot.species)

shlw_dem.like <- rsim.fit.table(scene_fit, run.shlw_dem)  

# what changed?            
View(shlw_dem.like - base.like)




##### 2.Fitting vulnerabilities of microbes ####################################
microbes <- c("pelagic_microbes", "benthic_microbes")
microbes_vartype <- c(rep("preyvul", length(microbes)), rep("predvul", length(microbes)))
microbes_species <- c(microbes, microbes)
microbes_start   <- rep(0,length(microbes_species))

ptm <- proc.time()
microbes_opt_AB <- optim(microbes_start, rsim.fit.run, method = "L-BFGS-B", 
                         lower = -4.5, upper = 4.5, hessian = FALSE, 
                         control = list(trace = 6, factr=1e11),
                         species=microbes_species, vartype=microbes_vartype,
                         scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
microbes_opt_AB

# Result: new vector for values  
microbes_opt_AB$par
#View(data.frame(microbes_species, microbes_vartype, microbes_start, microbes_opt_AB$par))

# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.update <- scene_fit
scene.update$params <- rsim.fit.apply(microbes_opt_AB$par, microbes_species, microbes_vartype, scene_fit$params)

View(data.frame(scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyFrom+1]],
                scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyTo+1]],
                scene.update$params$VV))

# Test that we get the correct likelihood in the new scene
rsim.fit.run(NA, NA, NA, 
             scene=scene.update, run_method="AB", run_years=hind_years, verbose=F)  
# But can also preserve parameters and run same thing with the original.    
microbes.result <- microbes_opt_AB$par
rsim.fit.run(microbes.result, microbes_species, microbes_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 

# Save the actual full run (same as above with verbose=T)
run.microbes <- rsim.fit.run(microbes.result, microbes_species, microbes_vartype, 
                             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.microbes, plot.species)

microbes.like <- rsim.fit.table(scene_fit, run.microbes)  

# what changed?            
View(microbes.like - base.like)





# ---------------------------------------------------------------------------- #
##### 2.Fitting vulnerabilities of all living groups #########################
living <- scene_fit$params$spname[2:74]
living_vartype <- c(rep("preyvul", length(living)), rep("predvul",length(living))) 
living_species <- c(living, living)
living_start   <- rep(0,length(living_species))

ptm <- proc.time()
living_opt_AB <- optim(living_start, rsim.fit.run, method = "L-BFGS-B", 
                       lower = -4.5, upper = 4.5, hessian = FALSE, 
                       control = list(trace = 6, factr=1e14),
                       species=living_species, vartype=living_vartype,
                       scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
living_opt_AB

# Result: new vector for values  
living_opt_AB$par
living_vv <- living_opt_AB$par
names(living_vv) <- living_species
#View(data.frame(living_species, living_vartype, living_start, living_opt_AB$par))

# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.update <- scene_fit
scene.update$params <- rsim.fit.apply(living_opt_AB$par, living_species, living_vartype, scene_fit$params)

View(data.frame(scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyFrom+1]],
                scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyTo+1]],
                scene.update$params$VV))

# Test that we get the correct likelihood in the new scene
rsim.fit.run(NA, NA, NA, 
             scene=scene.update, run_method="AB", run_years=hind_years, verbose=F)  
# But can also preserve parameters and run same thing with the original.    
living.result <- living_opt_AB$par
rsim.fit.run(living.result, living_species, living_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 

# Save the actual full run (same as above with verbose=T)
run.living <- rsim.fit.run(living.result, living_species, living_vartype, 
                           scene=scene_fit, run_method="AB", run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.living, plot.species)

living.like <- rsim.fit.table(scene_fit, run.living)  

# what changed?            
View(living.like - base.like)   








# ---------------------------------------------------------------------------- #
# Fitting vulnerabilities of zooplankton
plankton <- c("copepods","euphausiids","pelagic_microbes","benthic_microbes",
              "lg_phytoplankton", "sm_phytoplankton")
plankton_vartype <- c(rep("predvul", length(plankton)), rep("preyvul",length(plankton))) 
plankton_species <- c(plankton, plankton)
plankton_start   <- rep(0,length(plankton_species))

# Just to look at these vectors
#View(data.frame(plankton_species,plankton_vartype,plankton_start))

ptm <- proc.time()
plankton_opt_AB <- optim(plankton_start, rsim.fit.run, method = "L-BFGS-B", 
                         lower = -4.5, upper = 4.5, hessian = FALSE, 
                         control = list(trace = 6, factr=1e13),
                         species=plankton_species, vartype=plankton_vartype,
                         scene=scene_fit, run_method="AB", run_years=hind_years)
proc.time() - ptm
plankton_opt_AB

# Result: new vector for values  
plankton_opt_AB$par
plankton_vv <- plankton_opt_AB$par
names(plankton_vv) <- plankton_species
#View(data.frame(plankton_species, plankton_vartype, plankton_start, plankton_opt_AB$par))

# Apply the parameters to a new scene ("permanently").  Note: this loses the
# distinction between preyvul and predvul, unfortunately.
scene.update <- scene_fit
scene.update$params <- rsim.fit.apply(plankton_opt_AB$par, plankton_species, plankton_vartype, scene_fit$params)

View(data.frame(scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyFrom+1]],
                scene.update$params$spname[scene.update$params$spname[scene.update$params$PreyTo+1]],
                scene.update$params$VV))

# Test that we get the correct likelihood in the new scene
rsim.fit.run(NA, NA, NA, 
             scene=scene.update, run_method="AB", run_years=hind_years, verbose=F)  
# But can also preserve parameters and run same thing with the original.    
plankton.result <- plankton_opt_AB$par
rsim.fit.run(plankton.result, plankton_species, plankton_vartype, 
             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=F) 

# Save the actual full run (same as above with verbose=T)
run.plankton <- rsim.fit.run(plankton.result, plankton_species, plankton_vartype, 
                             scene=scene_fit, run_method="AB", run_years=hind_years, verbose=T) 
# Now let's look
rsim.runplot(scene_fit, run.plankton, plot.species)

plankton.like <- rsim.fit.table(scene_fit, run.plankton)  

# what changed?            
View(plankton.like - base.like)   



################################################################################

# Mzero fitting (not cleaned)      


fit_species <- rpath.living(bal)
fit_start   <- as.numeric(scene_fit$params$MzeroMort[fit_species])
fit_vartype <- rep("mzero",length(fit_species))


optim_out <- optim(fit_start, rsim.fit.run, method="Nelder-Mead", control=list(trace=3),
                   species=fit_species, vartype=fit_vartype,
                   scene=scene_fit, run_method="AB", run_years=hind_years)  

fup()

run_fit0 <- rsim.fit.run(fit_start, fit_species, fit_vartype, 
                         scene_fit, run_method="AB", run_years=hind_years, verbose=T)

run_fit1 <- rsim.fit.run(optim_out$par, fit_species, fit_vartype, 
                         scene_fit, run_method="AB", run_years=hind_years, verbose=T)

rsim.runplot(scene_fit, run_fit1, plot.species)

## FLOW PLOTTING ###############################################################
library(ggplot2)  
library(plotly)
library(tidyverse)

scene <- scene6
run   <- run_fit3
X11()
rsim.plot.full(scene6, run_fit3, "sablefish_adu")
# ATTEMPT 1 - works but ugly
# dlinks <- data.frame(t(run$annual_Qlink))
# colnames(dlinks) <- as.character(rownames(run$annual_Qlink))
# diet_table <- cbind(data.frame(
#                     from = as.character(scene$params$spname[scene$params$PreyFrom+1]),
#                     to = as.character(scene$params$spname[scene$params$PreyTo+1])),
#                     dlinks)
# sp <- "pollock_adu"
# pt <- diet_table[diet_table$to==sp,!(names(diet_table) =="to")]
# prey <- pt$from
# valframe <- pt[,!(names(pt) =="from")]
# flow <- as.vector(as.matrix(valframe))
# year <- as.numeric(names(valframe)[col(valframe)])
# diet_long <- data.frame(year,prey,flow)

# Attempt 2 - slightly less ugly
qt     <- t(run$annual_Qlink)
flow   <- as.vector(qt)
years  <- as.numeric(colnames(qt)[col(qt)])
from   <- as.character(scene$params$spname[scene$params$PreyFrom+1])
to     <- as.character(scene$params$spname[scene$params$PreyTo+1])
qframe <- data.frame(years, from, to, flow)

flowdat <- qframe %>%
  group_by(years,to) %>%
  mutate(tflow=sum(flow)) %>%
  ungroup() %>%
  mutate(dc = flow/tflow)

sp <- "pollock_adu"
#predframe <- qframe[qframe$to==sp,]
predframe <- flowdat %>% 
  filter(to==sp)

#X11()
#ggplot(data=predframe,aes(x = qyears,y=qvals,group=qfrom)) + 
#  geom_line(aes(color=qfrom)) 

fig1 <- plot_ly(predframe, x= ~years, y= ~dc, color= ~from)
fig1 <- fig1 %>% add_lines()

fig2 <- plot_ly(predframe, x= ~years, y= ~flow, color= ~from)
fig2 <- fig2 %>% add_lines()

fig <- subplot(fig1,fig2)
fig


t(diet_table[,as.character(hind_years)])

##### BELOW THIS NOT YET INCORPORATED INTO MAIN FITTING SEQUENCE  
#######################################################################################  


X11()
rsim.plot.fitbio.small(scene_fit, run_fit, "oth_rockfish",NA)

rsim.plot.full(scene_fit, run_fit, "gr_turbot_juv")

#render.show.fit(bal, scene_fit, run_fit, "test1")
################################################################################
# M0 Fitting (species one by one)
#  
# This function returns a list of biomass timeseries, species with timeseries, and combo
bio_series   <- rsim.fit.list.bio.series(scene_fit)
catch_series <- rsim.fit.list.catch.series(scene_fit)

# Set fitting weight of all biomass and catch timeseries to 0 
scene_zero <- scene_fit
scene_zero <- rsim.fit.set.bio.wt(scene_zero, bio_series$groups, bio_series$sources, 0)
scene_zero <- rsim.fit.set.catch.wt(scene_zero, catch_series$groups, 0)


this_survey <- "ebs_race_shelf"  
survey_spp  <-rsim.fit.list.bio.bysurvey(scene_fit)[[this_survey]]

final_vartype <- NULL
final_species <- NULL
final_values  <- NULL

for (sp in survey_spp){
  
  scene_m <- rsim.fit.set.bio.wt(scene_zero, sp, this_survey, 1.0)
  scene_m <- rsim.fit.set.q(scene_m, sp, this_survey, years=1991)
  
  # Vectors of parameters to fit.  Type of variable, species, and starting value  
  fit_vartype <- c("mzero") 
  fit_species <- c(sp)
  fit_start   <- as.numeric(scene_m$params$MzeroMort[sp])
  
  optim_out <- optim(fit_start, rsim.fit.run, method="Brent", lower=-5, upper=5,
                     species=fit_species, vartype=fit_vartype,
                     scene=scene_m, run_method="AB", run_years=hind_years)
  
  cat(sp," ",optim_out$par," ",optim_out$value,"\n"); flush.console()
  run_out <- rsim.fit.run(optim_out$par, fit_species, fit_vartype, scene_m, "AB", hind_years, verbose=T)
  X11(width=8,height=4)
  rsim.plot.full(scene_m, run_out, sp)    
  final_vartype <- c(final_vartype, fit_vartype)
  final_species <- c(final_species, fit_species)
  final_values  <- c(final_values, optim_out$par)
  
}  

###############################################################################
##
## BELOW THIS LINE ARE OUT-OF-ORDER TESTS NOT "FINAL" fitting sequence.
##
###############################################################################

# rsim.fit.run produces a run that applies the above fit_vectors to the
# scenario before running (scenario remains unchanged)
run_0 <- rsim.fit.run(fit_start, fit_species, fit_vartype, scene_f3, "AB", hind_years, verbose=T) 
X11(width=8,height=4)
rsim.plot.full(scene_f3,run_0,"yf_sole")  
#rsim.plot(run_0)

# By default, that biomass time series is an index (q internally calculated)
# Trying that fit.
optim_f3 <- optim(fit_start, rsim.fit.run, method="Brent", lower=-3, upper=3,
                  species=fit_species, vartype=fit_vartype,
                  scene=scene_f3, run_method="AB", run_years=hind_years)
optim_f3$par
optim_f3$value

run_f3 <- rsim.fit.run(optim_f3$par, fit_species, fit_vartype, scene_f3, "AB", hind_years, verbose=T) 
X11(width=8,height=4)
rsim.plot.full(scene_f3,run_f3,"yf_sole") 


# Trying fit with that same time series with a fixed q of 1 (very silly)
scene_f4 <- rsim.set.fit.q(scene_f3, "yf_sole", "ebs_race_shelf", q=1)
optim_f4 <- optim(fit_start, rsim.fit.run, method="Brent", lower=-3, upper=3,
                  species=fit_species, vartype=fit_vartype,
                  scene=scene_f4, run_method="AB", run_years=hind_years)
optim_f4$par
optim_f4$value

run_f4 <- rsim.fit.run(optim_f4$par, fit_species, fit_vartype, scene_f4, "AB", hind_years, verbose=T) 
X11(width=8,height=4)
rsim.plot.full(scene_f4,run_f4,"yf_sole") 

# Trying fit with that same time series with a fixed q using 1991 and Ecopath starting value
scene_f5 <- rsim.set.fit.q(scene_f3, "yf_sole", "ebs_race_shelf", years=1991)
optim_f5 <- optim(fit_start, rsim.fit.run, method="Brent", lower=-3, upper=3,
                  species=fit_species, vartype=fit_vartype,
                  scene=scene_f5, run_method="AB", run_years=hind_years)
optim_f5$par
optim_f5$value
run_f5 <- rsim.fit.run(optim_f5$par, fit_species, fit_vartype, scene_f5, "AB", hind_years, verbose=T) 
X11(width=8,height=4)
rsim.plot.full(scene_f5,run_f5,"yf_sole")   





scene_test <- rsim.set.fit.q(scene3, "sandlance", "ebs_race_shelf", years=1991)



run_f2   <- rsim.run(scene_f2, method='AB', years=hind_years)


rsim.fit.obj(scene_f2,run_f2,verbose=F)
rsim.fit.table(scene_f2,run_f2)

X11(width=8,height=4)
rsim.plot.full(scene_f2,run_f2,"sandlance")








fit.res <- rsim.fit.run(test$par,
                        species=fit_species, vartype=fit_vartype,
                        scene=scene_f2, run_method="AB", run_years=hind_years,
                        verbose=T)

X11(width=8,height=4)
rsim.plot.full(scene_f2,fit.res,"sandlance")
rsim.plot(fit.res)



scene_test <- rsim.set.fit.q(scene3, "sandlance", "ebs_race_shelf", years=1985)
scene_test <- rsim.set.fit.q(scene3, "sandlance", "ebs_race_shelf", years=1991:1993)

#
# #From a2 Repo
#
# # HISTORICAL CLIMATE
# # hind_clim <- read.csv("climate/A2hind.csv", row.names=NULL)
# hind_clim <- read.csv("climate/a2_hind.csv", row.names=NULL)
# hind_clim <- hind_clim[253:624,]
# BaseB        <- scene$params$B_BaseRef
# names(BaseB) <- scene$params$spname
# # Set groups to force   
# scene$forcing$ForcedBio[1:371,"euphausiids"]      <- BaseB["euphausiids"]      * hind_clim[2:372,"eup"]
# scene$forcing$ForcedBio[1:371,"copepods"]         <- BaseB["copepods"]         * hind_clim[2:372,"cop"]
# scene$forcing$ForcedBio[1:371,"pelagic_microbes"] <- BaseB["pelagic_microbes"] * hind_clim[2:372,"mzl"]
# scene$forcing$ForcedBio[1:371,"benthic_microbes"] <- BaseB["benthic_microbes"] * hind_clim[2:372,"mzl"]
# scene$forcing$ForcedBio[1:371,"lg_phytoplankton"] <- BaseB["lg_phytoplankton"] * hind_clim[2:372,"phl"]
# scene$forcing$ForcedBio[1:371,"sm_phytoplankton"] <- BaseB["sm_phytoplankton"] * hind_clim[2:372,"phs"]     


# # Set-up bioenergetic forcing
# bioen_sp_noceph <- bioen_sp[! bioen_sp %in% c("octopus","squids")]
# # print(bioen_sp_noceph)
# # consumption multiplier
# if (cons == TRUE) {
#   if(esm == "none" & ssp == "persist") {
#     # hindcast
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedSearch[1:372,i] <- tdc_hind_bt[253:624,i]
#     }
#     # projection
#     # do nothing? or default?
#     # for(i in bioen_sp_noceph){
#     #   scene$forcing$ForcedSearch[373:1308,i] <- 1
#     # }
#   }
#   
#   if(esm == "c" & ssp == 45){
#     # hindcast
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedSearch[1:372,i] <- tdc_hind_bt[253:624,i]
#     }
#     # projection
#     bioen_proj <- bioen_params(cmip, esm, ssp)
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedSearch[373:1068,i] <- bioen_proj[[1]][,i]
#     }
#     
#   }
#   else {
#     # hindcast
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedSearch[1:372,i] <- tdc_hind_bt[253:624,i]
#     }
#     # projection
#     bioen_proj <- bioen_params(cmip, esm, ssp)
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedSearch[373:1308,i] <- bioen_proj[[1]][,i]
#     }
#     
#   }
# }
# 
# # respiration multiplier ---------- ---------- #
# if (resp == TRUE) {
#   if(esm == "none" & ssp == "persist") {
#     # hindcast
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedActresp[1:372,i] <- tdr_hind_bt[253:624,i]
#     }
#     # projection
#     # do nothing? or default?
#     # for(i in bioen_sp_noceph){
#     #   scene$forcing$ForcedActresp[373:1308,i] <- 1
#     # }
#   }
#   
#   if(esm == "c" & ssp == 45){
#     # hindcast
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedActresp[1:372,i] <- tdr_hind_bt[253:624,i]
#     }
#     # projection
#     bioen_proj <- bioen_params(cmip, esm, ssp)
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedActresp[373:1068,i] <- bioen_proj[[2]][,i]
#     }
#   }
#   else {
#     # hindcast
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedActresp[1:372,i] <- tdr_hind_bt[253:624,i]
#     }
#     # projection
#     bioen_proj <- bioen_params(cmip, esm, ssp)
#     for(i in bioen_sp_noceph){
#       scene$forcing$ForcedActresp[373:1308,i] <- bioen_proj[[2]][,i]
#     }
#     # print(scene$forcing$ForcedActresp[373:1308,"pcod_adu"])
#   }
# }

