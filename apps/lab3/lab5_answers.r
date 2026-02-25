

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


# Preliminaries from previous labs
  library(Rpath)
  ebs.base.file  <- "EBS_condensed_563_base.csv"
  ebs.diet.file  <- "EBS_condensed_563_diet.csv"
  ebs.unbalanced <- read.rpath.params(ebs.base.file, ebs.diet.file)
  check.rpath.params(ebs.unbalanced)
  ebs.balanced <- rpath(ebs.unbalanced, eco.name="E Bering")

# Create an Ecosim Scenario, setting up parameters and
# creating forcing time series (with default, equlibrium values)
  ebs.scene0 <-  rsim.scenario(ebs.balanced, ebs.unbalanced, years=1990:2089)  

# Start by making an unfished system by turning off all gear effort
#  (we're assuming no fishing other than cod/pollock)
  all.gear       <- c("Trawl", "Cod.pots", "Longline", "Crab.pots", "State.fisheries") 
  ebs.unfished.scene <- adjust.fishing(ebs.scene0, "ForcedEffort", all.gear, 
                                       sim.year = 1992:2089, value=0)   
    
# create vectors of fishing levels
# cod are less productive than pollock so using lower F rates for cod
  eco.area <- 495000 # square km in the EBS model
  f_levels  <- seq(0, 2.0, 0.1)
  cf_levels <- seq(0, 1.0, 0.05)
# create some empty matricies for results
  cod_landings <- pol_landings <- matrix(NA,length(f_levels),length(cf_levels))
  price <- cost <- cod_bio <- pol_bio <- matrix(NA,length(f_levels),length(cf_levels))
  
  j<-0
  for (crate in cf_levels){
    j <- j+1
    # Set Cod F
      ebs.fscene <- adjust.fishing(ebs.unfished.scene, "ForcedFRate", "P.cod", 
                    sim.year = 1992:2089, value = crate)
    i<-0
    for (frate in f_levels){
      i <- i+1
      # Set Pollock F
        ebs.fscene <- adjust.fishing(ebs.fscene, "ForcedFRate", "Pollock", 
                      sim.year = 1992:2089, value = frate)
      # run the model
        ebs.frun    <- rsim.run(ebs.fscene, method="AB", years=1990:2089)
      
      # save results (end year is 2089)
        cod_landings[i,j] <- ebs.frun$annual_Catch["2089","P.cod"]
        pol_landings[i,j] <- ebs.frun$annual_Catch["2089","Pollock"]
        cod_bio[i,j]      <- ebs.frun$annual_Biomass["2089","P.cod"]
        pol_bio[i,j]      <- ebs.frun$annual_Biomass["2089","Pollock"]
        # Units of catch are t/km2.  Convert to kg and multiple by total km2
        # $0.50/kg for cod, $0.10/kg for pollock
          price[i,j] <- 1000*eco.area * (0.5*cod_landings[i,j] + 0.1*pol_landings[i,j])
        #  cost of fishing pollock is $20,000,000 per unit F, cod $80,000,000 per unit F.
           cost[i,j]  <- 20000000*frate + 80000000*crate 
      cat(crate,frate,"\n"); flush.console()
    }
  }
  
  grosstons  <- eco.area * (pol_landings+cod_landings)

  contour(grosstons,xlab="Pollock F", ylab="Cod F")
  contour(price,xlab="Pollock F", ylab="Cod F")
  contour(price-cost,xlab="Pollock F", ylab="Cod F")


