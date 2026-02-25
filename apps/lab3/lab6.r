
# A bigger maximization 
# (requires code in lab5_optim.r)
  ebs.f_rates <- c(Pollock=0, P.cod=0, Arrowtooth=0, P.halibut=0, Sm.flatfish=0, Forage.fish=0, Pred.zoop=0)
  ebs.prices  <- c(Pollock=1, P.cod=1, Arrowtooth=1, P.halibut=1, Sm.flatfish=1, Forage.fish=1, Pred.zoop=1)
  ebs.costs   <- c(Pollock=0, P.cod=0, Arrowtooth=0, P.halibut=0, Sm.flatfish=0, Forage.fish=0, Pred.zoop=0)
  result3 <- optim(ebs.f_rates, fishery_output, prices=ebs.prices, costs=ebs.costs, 
                            years=1990:2089, rsim.scene=ebs.unfished.scene, 
                            control=list(fnscale = -1), lower=0)

  result3
  run3 <- fishery_output(result3$par, ebs.prices, ebs.costs, 1990:2089, ebs.unfished.scene, verbose=T)
  rsim.plot(run3,all.species)
  catch3 <- 495000 * run3$annual_Catch["2089",] 
  barplot(catch3,las=2)   
  data.frame(catch1,catch2, catch3)

  bio.unfished <- run.unfished$annual_Biomass["2089",]
  biofrac1 <- run1$annual_Biomass["2089",]/bio.unfished
  biofrac2 <- run2$annual_Biomass["2089",]/bio.unfished
  biofrac3 <- run3$annual_Biomass["2089",]/bio.unfished

  data.frame(bio.unfished,biofrac1,biofrac2,biofrac3)
# 
# HOMEWORK:  Perform above "bigger maximization", subject to the constraint
# that no species drop below 20% of its unfished value.
# (hint, you will need to add input terms to the fishery_output function)
# Show the resulting run - what species is the choke (limiting) group?
#

runRpathShiny("lab3")









