# Update library
  library(devtools)
  install_github('NOAA-EDAB/Rpath', ref="lab_563_spring2020") 

  library(Rpath)

#####################################################################################

setwd("c:/src/sean_rpath/Rpath_public/inst/labs/lab4")

  base.file       <- "GOA_fit_base.csv"
  diet.file       <- "GOA_fit_diet.csv"
  pedigree.file   <- "GOA_fit_ped.csv"

  unbalanced <- read.rpath.params(base.file, diet.file, pedigree.file)
  check.rpath.params(unbalanced)
  balanced <- rpath(unbalanced, eco.name="Gulf of Alaska", eco.area=291840)
  scene0   <- rsim.scenario(balanced, unbalanced, years=1950:2049)
  run0     <- rsim.run(scene0,method="AB",1950:2015)
  rsim.plot(run0)
  
  # Show predprey parameter table
  get.rsim.predprey(scene0,"P.Cod")

  # Show parameter effects
  runRpathShiny("lab3")
  
  
  scene1 <- scene0

  # Read in biomass and catch time series for fitting
  scene1 <- read.fitting.biomass(scene1, "GOA_biomassdata.csv")
  scene1 <- read.fitting.catch(scene1, "GOA_catchdata.csv")
  run1   <- rsim.run(scene1,method="AB",1950:2014)

  par(mfrow=c(2,2))
  rsim.plot.catch(scene1,run1,"P.Cod")
  rsim.plot.biomass(scene1,run1,"P.Cod")
  rsim.plot.catch(scene1,run1,"W.Pollock")
  rsim.plot.biomass(scene1,run1,"W.Pollock")

  # Apply fitting Catch to fishing Catch
  scene1 <- apply.fit.to.catch(scene1)
  run1   <- rsim.run(scene1,method="AB",1950:2014)

  par(mfrow=c(2,2))
  rsim.plot.catch(scene1,run1,"P.Cod")
  rsim.plot.biomass(scene1,run1,"P.Cod")
  rsim.plot.catch(scene1,run1,"W.Pollock")
  rsim.plot.biomass(scene1,run1,"W.Pollock")

  # Now add Pacific Decadal Oscillation forcing
  PDO <- read.csv("PDO_norm_1950_2014.csv")
  PDO.scale <- 0.3
  scene1$forcing$ForcedBio[1:length(PDO$Value),"Outside"] <- 1.0 + PDO$Value * PDO.scale
  run1   <- rsim.run(scene1,method="AB",1950:2014)
  par(mfrow=c(2,2))
  rsim.plot.catch(scene1,run1,"P.Cod")
  rsim.plot.biomass(scene1,run1,"P.Cod")
  rsim.plot.catch(scene1,run1,"W.Pollock")
  rsim.plot.biomass(scene1,run1,"W.Pollock")

  # A full run including past and future, to show what happens after
  # catch forcing stops
  run1   <- rsim.run(scene1,method="AB",1950:2049)
  par(mfrow=c(2,2))
  rsim.plot.catch(scene1,run1,"P.Cod")
  rsim.plot.biomass(scene1,run1,"P.Cod")
  rsim.plot.catch(scene1,run1,"W.Pollock")
  rsim.plot.biomass(scene1,run1,"W.Pollock")


  # Show Cod Ration figure - bottom up control is key?

  # Show forcing method (state space)
  
  
  
# A mini Management Strategy Evaluation run
  # First get an unfished ecosystem state for reference
  scene.unfished <- adjust.fishing(scene0, "ForcedEffort", "Fishery", sim.year = 1950:2049, value=0)   
  run.unfished   <- rsim.run(scene.unfished, method="AB", years=1950:2049)


  # Setup the full scenario as above
  scene1 <- scene0
  scene1 <- read.fitting.biomass(scene1, "GOA_biomassdata.csv")
  scene1 <- read.fitting.catch(scene1, "GOA_catchdata.csv")  
  scene1 <- apply.fit.to.catch(scene1)
  PDO <- read.csv("PDO_norm_1950_2014.csv")
  PDO.scale <- 0.3
  scene1$forcing$ForcedBio[1:length(PDO$Value),"Outside"] <- 1.0 + PDO$Value * PDO.scale


  # Now run the "hindcast portion of the future run"
  run.future <- rsim.run(scene1,method="AB",1950:2014)
  rsim.plot(run.future)
  par(mfrow=c(2,2))
  rsim.plot.catch(scene1,run.future,"P.Cod")
  rsim.plot.biomass(scene1,run.future,"P.Cod")
  rsim.plot.catch(scene1,run.future,"W.Pollock")
  rsim.plot.biomass(scene1,run.future,"W.Pollock")  

  fished.species <- c("W.Pollock", "P.Cod", "Arrowtooth", "P.Halibut", "Sm.Flatfish")
  # some 2014 statistics to examine
  f_current     <- run.future$annual_Catch["2014",fished.species]/run.future$annual_Biomass["2014",fished.species]
  bfrac_current <- run.future$annual_Biomass["2014",fished.species]/run.unfished$annual_Biomass["2049",fished.species]

  
  # Make a list of the species of interest
  fished.species <- c("W.Pollock", "P.Cod", "P.Halibut", "Arrowtooth")
  #fished.species <- c("W.Pollock", "P.Cod", "Arrowtooth", "P.Halibut", "Sm.Flatfish")

  # Now we run the future one year at a time
  for (y in 2015:2049){
    # Measure the ecosystem state in the previous year  
    b_measured <- 1.1 * run.future$annual_Biomass[as.character(y-1), fished.species]
    c_measured <- run.future$annual_Catch[as.character(y-1), fished.species]
    f_measured <- c_measured/b_measured
    
    # Determine the desired F rate and therefore desired catch
    f_desired <- F_target_final # Added 40% Unfished biomass target
    c_desired  <- 1.1 * f_desired * b_measured  
 
    # Apply Catch to next year
    for (sp in fished.species){
      scene1 <- adjust.fishing(scene1, "ForcedCatch", sp, sim.year = y, value=c_desired[sp])   
    }

    # Run for one year
    run.future <- rsim.step(scene1, run.future, method="AB", year.end=y)     
  }

  # Examining the first results
  rsim.plot(run.future)
  par(mfrow=c(2,2))
  rsim.plot.catch(scene1,run.future,"P.Cod")
  rsim.plot.biomass(scene1,run.future,"P.Cod")
  rsim.plot.catch(scene1,run.future,"W.Pollock")
  rsim.plot.biomass(scene1,run.future,"W.Pollock")  
  scene1$fishing$ForcedCatch[,"P.Cod"]/run.future$annual_Biomass[,"P.Cod"]  




