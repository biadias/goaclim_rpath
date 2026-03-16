
# Update library
  library(devtools)
  install_github('NOAA-EDAB/Rpath', ref="lab_563_spring2020") 

# Move to directory with files
  #setwd()
  
# From previous labs:  load library, read in data
#                      check parameters, balance Ecopath model
  library(Rpath)
  ebs.base.file  <- "EBS_condensed_563_base.csv"
  ebs.diet.file  <- "EBS_condensed_563_diet.csv"
  ebs.unbalanced <- read.rpath.params(ebs.base.file, ebs.diet.file)
  check.rpath.params(ebs.unbalanced)
  ebs.balanced <- rpath(ebs.unbalanced, eco.name="E Bering")

# Create an Ecosim Scenario, setting up parameters and
# creating forcing time series (with default, equlibrium values)
  ebs.scene0 <-  rsim.scenario(ebs.balanced, ebs.unbalanced, years=1990:2089)

# Show lists
  # params
  # start_state
  # fishing
  # forcing

# some useful definitions and lists
  eco.area <- 495000
  fished.species <- c("Pollock","P.cod","Arrowtooth","P.halibut","Sm.flatfish",
                      "Rockfish","Comm.crabs")
  apex.species   <- c("Toothed.mammals","Baleen.whales","Birds","Apex.fish")
  other.species  <- c("Forage.fish","Epifauna","Benthos","Pred.zoop","Copepods",
                      "Pel.microbes","Ben.microbes","Lg.phytoplanton",  
                      "Sm.phytoplankton","Discards","Detritus")     
  all.species    <- c(fished.species,apex.species,other.species)
  all.gear       <- c("Trawl","Cod.pots","Longline","Crab.pots",
                      "State.fisheries") 

# Run the ecosystem scenario for indicated years
  ebs.run0 <- rsim.run(ebs.scene0, method="RK4", years=1990:2089)

# Plot and show ebs.run0 outputs
  rsim.plot(ebs.run0, all.species)

  summary(ebs.run0)
  # Biomass (monthly, annual)
  # Catch (monthly, annual)
  # who eats what
  # Catch comes from Effort (show forced Effort in scenario)
  # CATCH = q(species, gear) * B(species, time) * ForcedEffort(gear, time) +
  #                            B(species, time) * ForcedFRate(species, time) +
  #                                               ForcedCatch(species, time)
  # q(species, gear) is normalized such that Effort = 1 at Ecopath start point
  # q(species,gear) = Ecopath catch(species,gear)/Ecopath Biomass(species)
  # 
 extract.node(ebs.run0,"Pollock")

# Scenario 1:  increase fishing on pollock for 20 years
  ebs.scene1 <- adjust.fishing(ebs.scene0, "ForcedFRate", "Pollock", 
                 sim.year = 1992:2011, value=0.1)
  ebs.run1 <- rsim.run(ebs.scene1, method="RK4", years=1990:2089)
  rsim.plot(ebs.run1, all.species)
  rsim.plot(ebs.run1, other.species)
  rsim.plot(ebs.run1, fished.species)
  rsim.plot(ebs.run1, apex.species)

# Scenario 2: increase fishing on P.cod for 100 years (to new equilibrium)
  ebs.scene2 <- adjust.fishing(ebs.scene0, "ForcedFRate", "P.cod", 
                 sim.year = 1992:2089, value=0.1)
  ebs.run2 <- rsim.run(ebs.scene2, method="RK4", years=1990:2089)
  rsim.plot(ebs.run2, all.species)
  rsim.plot(ebs.run2, other.species)
  rsim.plot(ebs.run2, fished.species)
  rsim.plot(ebs.run2, apex.species)

# Scenario 3: Add noise to primary production through "Outside" prey forcing
  ebs.scene3 <- adjust.forcing(ebs.scene0, "ForcedPrey", "Outside", 
                sim.year=1992:1999, value=rnorm(length(1992:1999),1,0.1))
  ebs.run3 <- rsim.run(ebs.scene3, method="RK4", years=1990:2089)
  rsim.plot(ebs.run3, all.species)
  rsim.plot(ebs.run3, other.species)
  rsim.plot(ebs.run3, fished.species)
  rsim.plot(ebs.run3, apex.species)

# Scenario 4: Combine Scene 2 fishing with Scene 3 Primary Production 
  # Copy full sub-objects (fishing and forcing) to mix & match 
  # fishing scenarios with environmental scenarios
  ebs.scene4 <- ebs.scene0
  ebs.scene4$fishing <- ebs.scene2$fishing
  ebs.scene4$forcing <- ebs.scene3$forcing
  ebs.run4 <- rsim.run(ebs.scene4, method="RK4", years=1990:2089)
  rsim.plot(ebs.run4, all.species)
  rsim.plot(ebs.run4, other.species)
  rsim.plot(ebs.run4, fished.species)
  rsim.plot(ebs.run4, apex.species)

# Using Adams-Basforth instead of RK4 for integration method
# (better speed, less accurate transients)
  ebs.run4 <- rsim.run(ebs.scene4, method="RK4", years=1990:2089)
  ebs.run4a <- rsim.run(ebs.scene4, method="AB", years=1990:2089)

  plot(ebs.run4$out_Biomass[,"Pollock"])
  lines(ebs.run4a$out_Biomass[,"Pollock"],col="red")
  plot(ebs.run4$out_Biomass[,"Sm.phytoplankton"])
  lines(ebs.run4a$out_Biomass[,"Sm.phytoplankton"],col="red")


# An Unfished Ecosystem
  ebs.unfished.scene <- adjust.fishing(ebs.scene0, "ForcedEffort", all.gear, 
                                       sim.year = 1992:2089, value=0)  
  ebs.unfished.run <- rsim.run(ebs.unfished.scene, method="RK4", 
                               years=1990:2089)
  rsim.plot(ebs.unfished.run, all.species)
  rsim.plot(ebs.unfished.run, other.species)
  rsim.plot(ebs.unfished.run, fished.species)
  rsim.plot(ebs.unfished.run, apex.species)

# BIOMASS/YIELD PLOT
# Start by making an unfished system
  ebs.unfished.scene <- adjust.fishing(ebs.scene0, "ForcedEffort", all.gear, 
                                       sim.year = 1992:2089, value=0) 
# Now loop through a set of F rates for Pollock, fishing to end of 
# scenario each time.  Save end biomass and end catch.  
# Uses AB integration method for speed.
  f_levels <- seq(0, 2.0, 0.05)
  pollock_bio <- pollock_catch <- rep(NA,length(f_levels))
  i <- 0
  for (frate in f_levels){
    i <- i+1
    ebs.fscene <- adjust.fishing(ebs.unfished.scene, "ForcedFRate", "Pollock", 
                   sim.year = 1992:2089, value = frate)
    ebs.frun    <- rsim.run(ebs.fscene, method="AB", years=1990:2089)
    pollock_bio[i]   <- ebs.frun$annual_Biomass["2089","Pollock"]
    pollock_catch[i] <- ebs.frun$annual_Catch["2089","Pollock"]
    cat(frate,"\n"); flush.console()
  }
 
  plot(pollock_bio,pollock_catch)

# EXERCISE 1. 
# The Ecopath fishing rate for P.cod is 0.1593405.  Repeat the above loop
# with P.cod at this fishing rate, instead of at 0.


# EXERCISE 2.  Approximately (i.e. F for cod and pollock to the nearest 0.1, or
# by examination of a graph), what are:
#   (a) fishing rates of (cod,pollock) and the 
#   (b) biomass levels of each, as a fraction of unfished biomass, that:
#   - maximizes pollock catch
#   - maximizes cod catch
#   - maximizes total catch
#   - If pollock bring $0.10/kg and cod bring $0.50/kg, maximize dollars
#     (include the actual total of dollars)
#   - Additionally, if the cost of fishing pollock is $20,000,000 per
#     unit F, and cod is $80,000,000 per unit F.  Maximize net profit.


# HOMEWORK
# maximize total catch and net profit if a third species, arrowtooth, is
# included.  Arrowtooth Price is $0.05/kg, Unit F cost is $25,000,000 
# Hint:  You can try optim if you are familiar with it (or even if not!).


