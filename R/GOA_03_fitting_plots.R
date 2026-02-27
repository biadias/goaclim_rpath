
source("R/ecofitting_plot_multi.R")
#source("R/ecofitting_plot_multi_v2.R")
my_runs <- list(Base_61M02par, Bioen_61M02par, Base_59M03par, Bioen_63par, Bioen_59M03par)
#my_runs <- list(run_par59M02_bioen, run_par59M03_bioen,run_par59M02_base, run_par59M03_primprod, run_par59M03_full)
#run_names <- list("Bioen59M0_Pcod_Poll", "bioen59M0_Pcod_Poll_Arrow", "Base59M02_Pcod_Poll", "PP59M0_Pcod_Poll_Arrow", "PP&Bioen59M0_Pcod_Poll_Arrow")
run_names <- list("Base61M02_Pcod_Poll", "Bioen61M02_Pcod_Poll", "Base_59M03_Pcod_Poll_Arrow", "Bioen_63par","Bioen59M0_Pcod_Poll_Arrow")

plot.species <- c(rpath.living(bal),rpath.detrital(bal))
#plot.species2 <- c( "arrowtooth_flounder_adult" ,   "arrowtooth_flounder_juvenile" ,
#                   "pacific_cod_adult" ,            "pacific_cod_juvenile" ,
#                   "walleye_pollock_adult" ,        "walleye_pollock_juvenile" )

# Call the function
rsim.runplot.multi(scene_bioen, my_runs, plot.species, run_names=run_names)
rsim.runplot(scene_base, fit_results_63par$Base$final_run, plot.species, "63par")



# ---------------------------------------------------------------------------- #
# Run 1: 63 Parameters Fitting (Vulnerabilities) ####
# ---------------------------------------------------------------------------- #
my_runs <- list(fit_results_63par$Base$final_run, fit_results_63par$Bioen$final_run, 
                fit_results_63par$PrimProd$final_run, fit_results_63par$Full$final_run)
run_names <- list("Base 63 vul parameters", "Bioen 63 vul parameters", "PrimProd 6363 vul parameters", "Full 6363 vul parameters")

scene_63par <- rsim.runplot.multi(scene_bioen, my_runs, plot.species, run_names=run_names)



# ---------------------------------------------------------------------------- #
# Run 2: 61 Parameters Fitting (Vulnerabilities + M0 2) ####
# ---------------------------------------------------------------------------- #

my_runs <- list(results_61M02par$Base$final_run, results_61M02par$Bioen$final_run, 
                results_61M02par$PrimProd$final_run, results_61M02par$Full$final_run)
run_names <- list("Base61M02par (Pollock, Pcod)", "Bioen61M02par (Pollock, Pcod)", 
                  "PrimProd61M02par (Pollock, Pcod)", "Full61M02par (Pollock, Pcod)")

scene_61M02par <- rsim.runplot.multi(scene_base, my_runs, plot.species, run_names=run_names)

# ---------------------------------------------------------------------------- #
# Run 3: 59 Parameters Fitting (Vulnerabilities + M03) ####
# ---------------------------------------------------------------------------- #

my_runs <- list(results_59M03par$Base$final_run, results_59M03par$Bioen$final_run, 
                results_59M03par$PrimProd$final_run, results_59M03par$Full$final_run)
run_names <- list("Base59M03par (Pollock, Pcod, Arrowtooth)", "Bioen59M03par (Pollock, Pcod, Arrowtooth)", 
                  "PrimProd59M03par (Pollock, Pcod, Arrowtooth) ", "Full59M03par (Pollock, Pcod, Arrowtooth)")  

scene_59M03par <- rsim.runplot.multi(scene_base, my_runs, plot.species, run_names=run_names)


# ---------------------------------------------------------------------------- #
# Run 4: 59 Parameters Fitting (Vulnerabilities + M04) ####
# ---------------------------------------------------------------------------- #
my_runs <- list(results_59M04par$Base$final_run, results_59M04par$Bioen$final_run, 
                results_59M04par$PrimProd$final_run, results_59M04par$Full$final_run)
run_names <- list("Base59M04par (Pollock, Pcod, Arrowtooth, DWflatfish)  ", 
                  "Bioen59M04par (Pollock, Pcod, Arrowtooth, DWflatfish)  ", 
                  "PrimProd59M04par (Pollock, Pcod, Arrowtooth, DWflatfish)  ", 
                  "Full59M04par (Pollock, Pcod, Arrowtooth, DWflatfish)  ")

scene_59M04par <- rsim.runplot.multi(scene_base, my_runs, plot.species, run_names=run_names)


# ---------------------------------------------------------------------------- #
# Run 5: 57 Parameters Fitting (Vulnerabilities + M05) ####
# ---------------------------------------------------------------------------- #
my_runs <- list(results_57M05par$Base$final_run, results_57M05par$PrimProd$final_run, results_57M05par$Bioen$final_run,
                results_57M05par$Full$final_run)
run_names <- list(
  "Base57M05 (Pollock, Pcod, Arrowtooth, DWflatfish, POP)  ",
  "Bioen57M05par (Pollock, Pcod, Arrowtooth, DWflatfish, POP)  ",
  "PrimProd57M05par (Pollock, Pcod, Arrowtooth, DWflatfish, POP)  ",
  "Full57M05par (Pollock, Pcod, Arrowtooth, DWflatfish, POP)  ")


scene_57M05par <- rsim.runplot.multi(scene_base, my_runs, plot.species, run_names=run_names, inset=c(-0.4,0))


# ---------------------------------------------------------------------------- #
# Best models ####
# ---------------------------------------------------------------------------- #

my_runs <- list(results_59M04par$Base$final_run, 
                #results_57M05par$Base$final_run, 
                #results_57M05par$Bioen$final_run,
                results_59M04par$Bioen$final_run,
                results_61M02par$Bioen$final_run)
run_names <- list(
  "Base59M04 (Pollock, Pcod, Arrowtooth, DWflatfish)",
  #"Base57M05par (Pollock, Pcod, Arrowtooth, DWflatfish, POP)",
  #"Bioen57M05par (Pollock, Pcod, Arrowtooth, DWflatfish, POP)",
  "Bioen59M04par (Pollock, Pcod, Arrowtooth, DWflatfish)",
  "Bioen61M02par (Pollock, Pcod)")


best_models<- rsim.runplot.multi(scene_base, my_runs, plot.species, run_names=run_names)
