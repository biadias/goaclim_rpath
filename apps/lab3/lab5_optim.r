
# Update library
  library(devtools)
  install_github('NOAA-EDAB/Rpath', ref="lab_563_spring2020") 

library(Rpath)

# Function to maximize for Optim.  First argument f_rates is a vector
# of parameters to maximize function over.  Prices should be price/kg,
# and costs should be costs/f_rate.  The order of prices and costs must
# be the same order as f_rates.  The conversion scale (495000000)
# is hard-coded for the Bering Sea. 
# Return ourput is $ in the final year (assumed to be "in eqiuilibrium")

  fishery_output <- function(f_rates, prices, costs, years, rsim.scene, verbose=F){

    species_list <- names(f_rates); 
    final_year   <- as.character(years[length(years)])

    # Loop through species and apply f rates of each to the scenario 
    for (p in 1:length(f_rates)){
      rsim.scene <- adjust.fishing(rsim.scene, "ForcedFRate", species_list[p], 
                                   sim.year = years, value=f_rates[p])
    }

    # Run the model
    res <- rsim.run(rsim.scene, method="AB", years=years)

    # Extract landings from results
    landings <- res$annual_Catch[final_year,species_list]

    # This is the output to maximize
    output <- sum (495000000*(landings*prices) - (f_rates*costs))

    # verbose=TRUE returns the full run time series (default False)
    # otherwise just return the function value
    if (verbose){return(res)} else{return(output)}
  }


# Ecopath Setup
  setwd("inst/labs/lab3")
  ebs.base.file  <- "EBS_condensed_563_base.csv"
  ebs.diet.file  <- "EBS_condensed_563_diet.csv"
  ebs.unbalanced <- read.rpath.params(ebs.base.file, ebs.diet.file)
  check.rpath.params(ebs.unbalanced)
  ebs.balanced <- rpath(ebs.unbalanced, eco.name="E Bering")
  all.species <- ebs.balanced$Group[1:ebs.balanced$NUM_LIVING]

# Ecosim scenario setup (setting up unfished system)
  ebs.scene0   <- rsim.scenario(ebs.balanced, ebs.unbalanced, years=1990:2089)
# Set Effort of all gear to 0
  all.gear <- c("Trawl", "Cod.pots", "Longline", "Crab.pots", "State.fisheries") 
  ebs.unfished.scene <- adjust.fishing(ebs.scene0, "ForcedEffort", all.gear, 
                                       sim.year = 1990:2089, value=0)   

  run.unfished <- rsim.run(ebs.unfished.scene, method="AB", years=1990:2089)
  rsim.plot(run.unfished,all.species)


# TEST fishery_output function
# Vector of F rates to use for maximization, with species names and starting values
# Vectors of prices and costs ($/kg, and $/F) in same order as f_rates
  ebs.f_rates <- c(Pollock=0.2, P.cod=0.2)
  ebs.prices  <- c(Pollock=0.1, P.cod=0.5)
  ebs.costs   <- c(Pollock=20000000, P.cod=80000000) 

# Test function output
  fishery_output(ebs.f_rates, ebs.prices, ebs.costs, 1990:2089, ebs.unfished.scene)

# Get full run details from output
  run0 <- fishery_output(ebs.f_rates, ebs.prices, ebs.costs, 1990:2089, ebs.unfished.scene, verbose=T)
  rsim.plot(run0,all.species)
  barplot(run0$annual_Catch["2089",],las=2)


# Run optim #######
# fnscale = -1 makes this a maximization not a minimization 
# first two arguments are vector to maximize, and function to use
# Other inputs to fishery_output (prices, costs) must be given names here. 
  ebs.f_rates <- c(Pollock=0.2, P.cod=0.2)
  ebs.prices  <- c(Pollock=0.1, P.cod=0.5)
  ebs.costs   <- c(Pollock=20000000, P.cod=80000000) 
  result1 <- optim(ebs.f_rates, fishery_output, prices=ebs.prices, costs=ebs.costs, 
                            years=1990:2089, rsim.scene=ebs.unfished.scene, 
                            control=list(fnscale = -1))
  result1

  # test applying result's found parameters (result$par) to function
  fishery_output(result1$par, ebs.prices, ebs.costs, 1990:2089, ebs.unfished.scene)

  # plot resulting ecosystem time series
  run1 <- fishery_output(result1$par, ebs.prices, ebs.costs, 1990:2089, ebs.unfished.scene, verbose=T)
  rsim.plot(run1,all.species) 
  catch1 <- 495000 * run1$annual_Catch["2089",]  # catch in tons in final year
  barplot(catch1,las=2)   

                     
# Arrowtooth to maximization 
  ebs.f_rates <- c(Pollock=0.2, P.cod=0.2, Arrowtooth=0.0)
  ebs.prices  <- c(Pollock=0.1, P.cod=0.5, Arrowtooth=0.05)
  ebs.costs   <- c(Pollock=20000000, P.cod=80000000, Arrowtooth=25000000) 
  result2 <- optim(ebs.f_rates, fishery_output, prices=ebs.prices, costs=ebs.costs, 
                            years=1990:2089, rsim.scene=ebs.unfished.scene, 
                            control=list(fnscale = -1))
  result2
  ## OOPS NEED TO FIX SOMETHING!


# Attempt 2, adding lower bound for f rates
  ebs.f_rates <- c(Pollock=0.2, P.cod=0.2, Arrowtooth=0.0)
  ebs.prices  <- c(Pollock=0.1, P.cod=0.5, Arrowtooth=0.05)
  ebs.costs   <- c(Pollock=20000000, P.cod=80000000, Arrowtooth=25000000)
  result2 <- optim(ebs.f_rates, fishery_output, prices=ebs.prices, costs=ebs.costs, 
                            years=1990:2089, rsim.scene=ebs.unfished.scene, 
                            control=list(fnscale = -1), lower=0)
  result2
  run2 <- fishery_output(result2$par, ebs.prices, ebs.costs, 1990:2089, ebs.unfished.scene, verbose=T)
  rsim.plot(run2,all.species)
  catch2 <- 495000 * run2$annual_Catch["2089",] 
  barplot(catch2,las=2)   

  data.frame(catch1,catch2)
  # NOTE:  Compare to actual (2M MT is good reference for Bering)
  # 2020 Overfishing Limit OFL:  Pollock 3,914,000         TAC 1,397,000
  #                                  Cod   216,000         TAC   166,475
  #                                  ATF    82,939 (BSAI)  TAC     8,000   
   

# CLASS EXERCISE
# Maximize total yield (in tons) for all of the following species simultaneously:
# Pollock, P.cod, Arrowtooth, P.halibut, Sm.flatfish, Forage.fish, Pred.zoop
# (Pred.zoop represents a krill fishery)

