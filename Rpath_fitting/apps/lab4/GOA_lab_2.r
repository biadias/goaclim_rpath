
# Find Frates to get close to target biomass
  target_F <- function(f_rates, b_target, years, rsim.scene, verbose=F){

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
    b_actual <- res$annual_Biomass[final_year,species_list]

    # This is the output to minimize
    output <- sum((b_actual-b_target) * (b_actual-b_target))

    # verbose=TRUE returns the full run time series (default False)
    # otherwise just return the function value
    if (verbose){return(res)} else{return(output)}
  }

# First Attempt
fished.species <- c("W.Pollock", "P.Cod", "P.Halibut", "Arrowtooth", "Sm.Flatfish")
goa.f_rates    <- rep(0,length(fished.species)); names(goa.f_rates) <- fished.species
goa.targets    <- 0.40 * run.unfished$annual_Biomass["2049",fished.species]
target_F(goa.f_rates, b_target=goa.targets, years=1950:2049, rsim.scene=scene.unfished)

result  <- optim(goa.f_rates, target_F, 
                 b_target=goa.targets, years=1950:2049, rsim.scene=scene.unfished,
                 control=list(fnscale = 1), lower=0)

run.result <- target_F(result$par, b_target=goa.targets, years=1950:2049, rsim.scene=scene.unfished, verbose=T)
bfrac <- run.result$annual_Biomass["2049",fished.species]/run.unfished$annual_Biomass["2049",fished.species]


# Second Attempt
fished.species <- c("W.Pollock", "P.Cod", "P.Halibut", "Arrowtooth")
goa.f_rates    <- rep(0,length(fished.species)); names(goa.f_rates) <- fished.species
goa.targets    <- 0.40 * run.unfished$annual_Biomass["2049",fished.species]
target_F(goa.f_rates, b_target=goa.targets, years=1950:2049, rsim.scene=scene.unfished)

result  <- optim(goa.f_rates, target_F, 
                 b_target=goa.targets, years=1950:2049, rsim.scene=scene.unfished,
                 control=list(fnscale = 1), lower=0)

run.result <- target_F(result$par, b_target=goa.targets, years=1950:2049, rsim.scene=scene.unfished, verbose=T)
bfrac <- run.result$annual_Biomass["2049",fished.species]/run.unfished$annual_Biomass["2049",fished.species]


fished.species <- c("W.Pollock", "P.Cod", "P.Halibut", "Arrowtooth")

# MEASUREMENT ERROR (surveys)

# IMPLEMENTATION ERROR (catch)

# PROCESS ERROR (CLIMATE)


F_target_final <- result$par

# Change Line 96 from
# fished.species <- c("W.Pollock", "P.Cod", "Arrowtooth", "P.Halibut", "Sm.Flatfish")
#
#to
# fished.species <- c("W.Pollock", "P.Cod", "P.Halibut", "Arrowtooth")

# Change Line 106 in GOA_lab.r from
#    f_desired  <- f_measured
# to
#   f_desired <- F_target_final 




