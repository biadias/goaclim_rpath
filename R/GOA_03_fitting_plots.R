
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