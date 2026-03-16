
library(Rpath)

# Update library
  library(devtools)
  install_github('NOAA-EDAB/Rpath', ref="lab_563_spring2020") 

# Move to directory with files
  setwd()
  
# csv files containing base data and diet matrix
  ebs.base.file  <- "EBS_condensed_563_base.csv"
  ebs.diet.file  <- "EBS_condensed_563_diet.csv"

# Examine csv files (note: type, detritus)
  
# read csv files into an "unbalanced ecopath" object
  ebs.unbalanced <- read.rpath.params(ebs.base.file, ebs.diet.file)

# Check unbalanced ecopath for format (note: test format breakage)
  check.rpath.params(ebs.unbalanced)

# Examine help files (note input file types: pedigree, stanzas)
  ?read.rpath.params  
# Alternate method: 
  # create.rpath.params() creates a skeleton set of parameters (all NAs) 
  # given a list of groups and group types
  
# Examine unbalanced Ecopath object (note: model, diet, pedigree)
  ebs.unbalanced
  
# Balance the model (solve for missing parameters)  
  ebs.balanced   <- rpath(ebs.unbalanced, eco.name="Eastern Bering Sea", eco.area="495000")

# Examine the balanced object
  ebs.balanced
  
# testing balance
   # Look up pollock (row names to be added in later version)
     pollock_index <- which(ebs.unbalanced$model$Group == "Pollock")
     poll_biomass  <- ebs.unbalanced$model$Biomass[pollock_index]

   # Set pollock biomass to be 10% of original
     ebs.unbalanced$model$Biomass[pollock_index] <- poll_biomass/10
     
   # try balancing he model again
     ebs.balanced   <- rpath(ebs.unbalanced, eco.name="Eastern Bering Sea", eco.area="495000")
     
   # put things back the way they were and rebalance
     ebs.unbalanced$model$Biomass[pollock_index] <- poll_biomass
     ebs.balanced   <- rpath(ebs.unbalanced, eco.name="Eastern Bering Sea", eco.area="495000")


# SOME Interesting Outputs          
     # Trophic level
       data.frame(ebs.balanced$Group,ebs.balanced$TL)

     # Fisheries - copied for easy renaming of rows/columns until update
       landings <- ebs.balanced$Landings
         rownames(landings)<-ebs.balanced$Group; 
         colnames(landings)<-ebs.balanced$Group[(ebs.balanced$NUM_GROUPS-ebs.balanced$NUM_GEARS+1):ebs.balanced$NUM_GROUPS]
         
        discards <- ebs.balanced$Discards
         rownames(discards)<-ebs.balanced$Group; 
         colnames(discards)<-ebs.balanced$Group[(ebs.balanced$NUM_GROUPS-ebs.balanced$NUM_GEARS+1):ebs.balanced$NUM_GROUPS]
         
       landingsTL <- colSums(landings * ebs.balanced$TL)/colSums(landings)
       discardsTL <- colSums(discards * ebs.balanced$TL)/colSums(discards)
       totTL      <- colSums((landings+discards) * ebs.balanced$TL)/colSums((landings+discards))

       
# EXERCISE
  # Reset the Rpath model
    ebs.unbalanced <- read.rpath.params(ebs.base.file, ebs.diet.file)
    ebs.balanced   <- rpath(ebs.unbalanced, eco.name="Eastern Bering Sea", eco.area="495000")

# Lookup for species of interest  
    species.index <- which(ebs.balanced$Group == "Pollock")
    
  # Lookups for primary productivity, forage fish, and zoop   
    prim.prod.index <- which(ebs.balanced$type==1)
    forage.index <- which(ebs.balanced$Group == "Forage.fish")
    zoop.index   <- which(ebs.balanced$Group %in% c("Pred.zoop","Copepods")) 
    
  # Set up a vector of biomass values to test  
    base.biomass <- ebs.balanced$Biomass[species.index]
    test.biomass <- seq(0.1, 5.0 , 0.1) * base.biomass
    
  # vectors to store outputs 
    prim.prod <- forage.bio <- zoop.bio <- rep(NA,length(test.biomass))
    
  # Loop through biomasses and calculate balance for each  
    ind <- 0
    for (b in test.biomass){
      ind <- ind + 1
      # Insert new biomass into unbalanced model
        ebs.unbalanced$model$Biomass[species.index] <- b
      # Make a new balanced model
        new.balance <- rpath(ebs.unbalanced)
      # Sum B * PB for Primary Producers
        prim.prod[ind]  <- sum(new.balance$Biomass[prim.prod.index] * new.balance$PB[prim.prod.index])
      # Forage fish biomass
        forage.bio[ind] <- new.balance$Biomass[forage.index]
      # Sum zooplankton biomass
        zoop.bio[ind] <- sum(new.balance$Biomass[zoop.index])
    }
    # Reset unbalanced model
      ebs.unbalanced$model$Biomass[species.index] <- base.biomass
      
# 1.  What is the relationship between pollock biomass and forage fish to support? zooplankton?
  plot(test.biomass,forage.bio)   
  plot(test.biomass,zoop.bio)
  
# 2.  What about primary productivity?
  plot(test.biomass,prim.prod) 
          
# 3.  Are there any EEs that are too high in the final run (with highest pollock)?
  data.frame(new.balance$Group, new.balance$EE)
  
# 4.  What is the highest pollock biomass (2 sig. figures) that has all EEs < 1  
      ebs.unbalanced$model$Biomass[species.index] <- 100
      new.balance <- rpath(ebs.unbalanced)
      ebs.unbalanced$model$Biomass[species.index] <- base.biomass
      data.frame(new.balance$Group, new.balance$EE)
      
# Pick another group near the top (P.cod, P.halibut, Arrowtooth, Birds, Toothed.mammals)
# and repeat 1-4 (if no EEs >= 1 for step 4, increase biomass until one goes above 1)
     
      
# DOES THE ABOVE TELL YOU WHAT HAPPENS WHEN BIOMASS ACTUALLY CHANGES?

     