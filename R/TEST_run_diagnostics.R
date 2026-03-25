
# The inputs needed for the following script are:
#   scene - an rsim scenario 
#   years - the years used to set up the scene

# Kerim-specific (path-specific) loading of scene and years
  #simple unfit model for testing
   # source("R/EBS_fitting_setup.R") 
   #  scene <- scene_fit
   #  years <- hind_years

  # Fit EBS model
    curdir <- getwd()
    setwd("../aclim3_rpath")
    source("R/load_EBS_par34_fits.R")
    setwd(curdir)
    
    scene <- scene_par34M03_full
    years <- all_years  
    
    
# If fitting vectors need to be applied, can add that here
#  scene$params <- rsim.fit.apply(values, species, vartype, scene)

# DIET OUTPUTS #################################################################
  drun <- rsim.run(scene, "AB", years)
  from_species <- scene$params$spname[scene$params$PreyFrom + 1] #+1 here is due to R to C++ array conversion
  to_species   <- scene$params$spname[scene$params$PreyTo   + 1] #+1 here is due to R to C++ array conversion
  
# Get flows to a single predator.  Note: annual_Qlink is a single snapshot at a point each year 
# for each pred/prey link, taken at the measure month (same as annual biomass etc).
  predator <- "pollock_adu"
  
  flow_mat           <- drun$annual_Qlink[ , which(to_species==predator)]
  colnames(flow_mat) <- from_species[which(to_species==predator)]

# Normalize to diet proportions  
  diet_mat <- flow_mat/rowSums(flow_mat)
  
# Flow and Diet proportion matrices over time for that predator, I'll leave plotting to experts...
  #View(flow_mat)
  #View(diet_mat)
  

# OUTPUTTING DERIVATIVE PARTS ##################################################
# This needs to use rsim_step to get annual increments from the state vector
# at the end of each year, since values aren't saved during the run.
  
# Doing this for 1 species only, could be modified to do more/all species
  this.species <- "bairdi"
  
# run the first year using rsim.run
  step_run <- rsim.run(scene, method="AB", years=years[1])
# TODO being lazy and not storing first derivative step, will modify later
  
# step through all years except the first year, saving results from dyt (main derivative) in output data frame  
  output <- data.frame(matrix(NA, length(years), 13)); row.names(output) <- years
  for (yy in years[-1]){ # years[-1] does all years except the first
    step_run <- rsim.step(scene, step_run, method="AB", year.end=yy)
    
  # deriv.table contains the derivative values for all species for a particular
  # timestep - here we're only saving 1 species in output, but others could be 
  # extracted from this table.
    deriv.table <- data.frame(step_run$end_state["Biomass"],
                              step_run$dyt[c("TotGain", "TotLoss", "DerivT", "FoodLoss", "FoodGain", "UnAssimLoss",
                                             "ActiveRespLoss", "DetritalGain", "FishingGain", "MzeroLoss", "FishingLoss", "DetritalLoss")], 
                              row.names=step_run$params$spname)
    
  # save the one species we want in the output data frame
    output[as.character(yy),] <- deriv.table[this.species,]
  }
  colnames(output) <- colnames(deriv.table)
  
# Convert flows to rates by dividing by biomass - otherwise all the trends just look like biomass
  output[,2:13] <- output[,2:13]/output$Biomass
  
# Some Plots
  plot(rownames(output),log(output$Biomass),xlab="year", ylab="log Biomass")
  plot(rownames(output), output$FishingLoss, xlab="year", ylab="F rate")
  plot(rownames(output), output$FoodGain - output$UnAssimLoss - output$ActiveRespLoss, xlab="year", ylab="Net feeding rate")
  
  
# FLOW MATRIX #################################################################
  drun <- rsim.run(scene, "AB", years)
  
# lookup vectors
  from_species <- scene$params$spname[scene$params$PreyFrom + 1] #+1 here is due to R to C++ array conversion
  to_species   <- scene$params$spname[scene$params$PreyTo   + 1] #+1 here is due to R to C++ array conversion
  link_lookup <- matrix(c(from_species,to_species),ncol=2) # lookup maps link num to prey/pred matrix by name
  
# Blank matrix with species names in rows/cols
  na_flow_mat <- matrix(NA, length(scene$params$spname), length(scene$params$spname))
  rownames(na_flow_mat) <- scene$params$spname
  colnames(na_flow_mat) <- scene$params$spname  

# vulnerability list and matrix
  vul <- scene$params$VV
  vlist <- data.frame(from_species,to_species,vul) # human readable list of vuls  
  vulmat <- na_flow_mat
  vulmat[link_lookup] <- vul # pred/prey matrix of vuls
  
# Flow matrix over time
    base_flows <- na_flow_mat
    base_flows[link_lookup] <- drun$annual_Qlink[as.character(years[1]),]
    
  for(yy in years){
    year_flows <- na_flow_mat
    year_flows[link_lookup] <- drun$annual_Qlink[as.character(yy),] 
    flow_anom <- year_flows/base_flows
    image(flow_anom,zlim=c(0,2), col=hcl.colors(30, palette = "Blue-Red"),xaxt="n",yaxt="n",xlab=yy)
  }
    
# as a GIF   
library(magick)
    gif_filebase <- "BS_full"
    for(yy in years){
      year_flows <- na_flow_mat
      year_flows[link_lookup] <- drun$annual_Qlink[as.character(yy),] 
      flow_anom <- year_flows/base_flows
      fp <- file.path("tmp", paste0(gif_filebase, yy, ".png"))
      png(fp, width=600,height=600)
      image(flow_anom,zlim=c(0,2), col=hcl.colors(30, palette = "Blue-Red"),xaxt="n",yaxt="n",xlab=yy)
      dev.off()
    }    

    imgs <- list.files("tmp", pattern=paste0(gif_filebase,"*"),full.names = TRUE)
    img_list <- lapply(imgs, image_read)
    img_joined <- image_join(img_list)
    img_animated <- image_animate(img_joined, fps = 2)
    #img_animated
    image_write(image = img_animated, path = paste0("images/", gif_filebase, ".gif"))

        
# AGE STRUCTURE #################################################################
    
    # run the first year using rsim.run
    step_run <- rsim.run(scene, method="AB", years=years[1])
    # TODO being lazy and not storing first derivative step, will modify later
    
    # step through all years except the first year, saving results from dyt (main derivative) in output data frame  
    output <- data.frame(matrix(NA, length(years), 13)); row.names(output) <- years
    for (yy in years[-1]){ # years[-1] does all years except the first
      step_run <- rsim.step(scene, step_run, method="AB", year.end=yy)
      
      # deriv.table contains the derivative values for all species for a particular
      # timestep - here we're only saving 1 species in output, but others could be 
      # extracted from this table.
      deriv.table <- data.frame(step_run$end_state["Biomass"],
                                step_run$dyt[c("TotGain", "TotLoss", "DerivT", "FoodLoss", "FoodGain", "UnAssimLoss",
                                               "ActiveRespLoss", "DetritalGain", "FishingGain", "MzeroLoss", "FishingLoss", "DetritalLoss")], 
                                row.names=step_run$params$spname)
      
      # save the one species we want in the output data frame
      #output[as.character(yy),] <- deriv.table[this.species,]
    }
    #colnames(output) <- colnames(deriv.table)
    
    
    
    