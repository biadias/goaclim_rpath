
# Ecopath Setup
  setwd("c:/src/sean_rpath/rpath_public/inst/labs/lab3")
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

# get 20% limit
  ebs.limits <- 0.2 * run.unfished$annual_Biomass["2089",all.species]


# function with limit added
  fishery_output_limit <- function(f_rates, prices, costs, limits, years, rsim.scene, verbose=F){

    species_list <- names(f_rates); 
    final_year   <- as.character(years[length(years)])

    # Loop through species and apply f rates of each to the scenario
    if (length(f_rates)>1){ 
     for (p in 1:length(f_rates)){
      rsim.scene <- adjust.fishing(rsim.scene, "ForcedFRate", species_list[p], 
                                   sim.year = years, value=f_rates[p])
     }
    }
    else{
      rsim.scene <- adjust.fishing(rsim.scene, "ForcedFRate", species_list, 
                                   sim.year = years, value=f_rates)
    }
    # Run the model
    res <- rsim.run(rsim.scene, method="AB", years=years)

    # Extract landings from results
    landings <- res$annual_Catch[final_year,species_list]

    limit_list  <- names(limits)
    end_biomass <- res$annual_Biomass[final_year,limit_list]
    penalty <- 1e15 * sum(1 - 1/(1+exp(-100*(end_biomass-limits))))

    # This is the output to maximize
    output <- sum (495000000*(landings*prices) - (f_rates*costs)) - penalty

    # verbose=TRUE returns the full run time series (default False)
    # otherwise just return the function value
    if (verbose){return(res)} else{return(output)}
  }

  ebs.f_rates <- c(Pollock=0, P.cod=0, Arrowtooth=0, P.halibut=0, Sm.flatfish=0, Forage.fish=0, Pred.zoop=0)
  ebs.prices  <- c(Pollock=1, P.cod=1, Arrowtooth=1, P.halibut=1, Sm.flatfish=1, Forage.fish=1, Pred.zoop=1)
  ebs.costs   <- c(Pollock=0, P.cod=0, Arrowtooth=0, P.halibut=0, Sm.flatfish=0, Forage.fish=0, Pred.zoop=0)
  result3 <- optim(ebs.f_rates, fishery_output_limit, prices=ebs.prices, costs=ebs.costs, limits=ebs.limits,
                            years=1990:2089, rsim.scene=ebs.unfished.scene, 
                            control=list(fnscale = -1), lower=0)

  result3
  run3 <- fishery_output_limit(result3$par, ebs.prices, ebs.costs, ebs.limits, 1990:2089, ebs.unfished.scene, verbose=T)
  outval <- fishery_output_limit(result3$par, ebs.prices, ebs.costs, ebs.limits, 1990:2089, ebs.unfished.scene)
  rsim.plot(run3,all.species)

  unfished.biomass <- run.unfished$annual_Biomass["2089",all.species]
  fished.biomass <- run3$annual_Biomass["2089",all.species]
  data.frame(unfished.biomass,fished.biomass,fished.biomass/unfished.biomass)

  result_combo <- list()
  for (sp in c("Pollock","P.cod","Arrowtooth","P.halibut","Sm.flatfish","Forage.fish","Pred.zoop")){
    ebs.f_rates <- 0;  ebs.prices <- 1;  ebs.costs <- 0;
    names(ebs.f_rates) <- sp
    result_combo[[sp]] <- optim(ebs.f_rates, fishery_output_limit, prices=ebs.prices, costs=ebs.costs, 
                            limits=ebs.limits, years=1990:2089, rsim.scene=ebs.unfished.scene, 
                            control=list(fnscale = -1), lower=0)
  }


  ebs.f_rates <- 0;  ebs.prices <- 1;  ebs.costs <- 0;
  names(ebs.f_rates) <- "Pollock"
  outval <- fishery_output_limit(ebs.f_rates, ebs.prices, ebs.costs, ebs.limits, 1990:2089, ebs.unfished.scene)

