
# The inputs needed for the following script are:
#   scene - an rsim scenario 
#   years - the years used to set up the scene
  #source("Rpath_fitting/R/EBS_fitting_setup.r")
  source("code/00_setup_forecast.R")
  best_model <- readRDS("Rpath_fitting/GOA/GOA_fit_results_59M04par.rds") # Loads the best model scenario and sets up scene_bioen
  
  best_model <- scene_gfdl_persist
  
  #scene <- scene_bioen
  scene <- scene_gfdl_persist
  scene <- scene_gfdl_126
  scene <- scene_gfdl_245
  scene <- scene_gfdl_585
  scene <- scene_gfdl_persist_cons
  scene <- scene_gfdl_persist_res
  #years <- hind_years
  years <- all_years

# If fitting vectors need to be applied, can add that here
#  scene$params <- rsim.fit.apply(values, species, vartype, scene)

# DIET OUTPUTS #################################################################
  drun <- rsim.run(scene, "AB", years)
  from_species <- scene$params$spname[scene$params$PreyFrom + 1] #+1 here is due to R to C++ array conversion
  to_species   <- scene$params$spname[scene$params$PreyTo   + 1] #+1 here is due to R to C++ array conversion
  
# Get flows to a single predator.  Note: annual_Qlink is a single snapshot at a point each year 
# for each pred/prey link, taken at the measure month (same as annual biomass etc).
  #predator <- "walleye_pollock_adult"
  predator <- "arrowtooth_flounder_adult"
  
  flow_mat           <- drun$annual_Qlink[ , which(to_species==predator)]
  colnames(flow_mat) <- from_species[which(to_species==predator)]

# Normalize to diet proportions  
  diet_mat <- flow_mat/rowSums(flow_mat)
  
# Flow and Diet proportion matrices over time for that predator, I'll leave plotting to experts...
  View(flow_mat)
  View(diet_mat)
  

# OUTPUTTING DERIVATIVE PARTS ##################################################
# This needs to use rsim_step to get annual increments from the state vector
# at the end of each year, since values aren't saved during the run.
  
# Doing this for 1 species only, could be modified to do more/all species
  this.species <- "arrowtooth_flounder_adult"
  
# run the first year using rsim.run
  #step_run <- rsim.run(scene, method="AB", years=years[1])
# TODO being lazy and not storing first derivative step, will modify later

# step through all years except the first year, saving results from dyt (main derivative) in output data frame  
  output <- data.frame(matrix(NA, length(years), 13)); row.names(output) <- years
  
  #for (yy in years[-1]){ # years[-1] does all years except the first
  for (i in 1:length(years)){
    current_years <- years[1:i]
   #step_run <- rsim.step(scene, step_run, method="AB", year.end=yy)
    step_run <- rsim.run(scene, method="AB", years=current_years)
  # deriv.table contains the derivative values for all species for a particular
  # timestep - here we're only saving 1 species in output, but others could be 
  # extracted from this table.
    deriv.table <- data.frame(step_run$end_state["Biomass"],
                              step_run$dyt[c("TotGain", "TotLoss", "DerivT", "FoodLoss", "FoodGain", "UnAssimLoss",
                                             "ActiveRespLoss", "DetritalGain", "FishingGain", "MzeroLoss", "FishingLoss", "DetritalLoss")], 
                              row.names=step_run$params$spname)
    
  # save the one species we want in the output data frame
    #output[as.character(yy),] <- deriv.table[this.species,]
    output[as.character(years[i]),] <- deriv.table[this.species,]
  }
  colnames(output) <- colnames(deriv.table)
  
# Convert flows to rates by dividing by biomass - otherwise all the trends just look like biomass
  output[,2:13] <- output[,2:13]/output$Biomass
  
# Some Plots
  plot(rownames(output),log(output$Biomass),xlab="year", ylab="log Biomass", main=this.species)
  plot(rownames(output), output$FishingLoss, xlab="year", ylab="F rate", main=this.species)
  plot(rownames(output), output$FoodGain - output$UnAssimLoss - output$ActiveRespLoss, xlab="year", ylab="Net feeding rate", main=this.species )
  
# ---------------------------------------------------------------------------- #
# several species
# ---------------------------------------------------------------------------- #
  
  
  target_species <- c("walleye_pollock_adult", "pacific_cod_adult", "arrowtooth_flounder_adult")
  
  output_list <- list()
  for(sp in target_species){
    output_list[[sp]] <- data.frame(matrix(NA, length(years), 13)); row.names(output_list[[sp]]) <- years}
    for (i in 1:length(years)){
      
      current_years <- years[1:i]
      step_run <- rsim.run(scene, method="AB", years=current_years)
      
      deriv.table <- data.frame(step_run$end_state["Biomass"],
                                step_run$dyt[c("TotGain", "TotLoss", "DerivT", 
                                               "FoodLoss", "FoodGain", "UnAssimLoss",
                                               "ActiveRespLoss", "DetritalGain", 
                                               "FishingGain", "MzeroLoss",
                                               "FishingLoss", "DetritalLoss")], 
                                row.names=step_run$params$spname)
      
      for (sp in target_species){
      output_list[[sp]][as.character(years[i]),] <- deriv.table[sp,]
      }
    }

  # ---------------------------------------------------------------------------- #
  # several species plots
  # ---------------------------------------------------------------------------- #
 
  par(mfrow= c(length(target_species), 3))
  par(mar= c(4,4,3,1))
  
  for(sp in target_species){
    out <- output_list[[sp]]
    
    colnames(out) <- colnames(deriv.table)
    out[,2:13] <- out[,2:13]/out$Biomass
    
    plot(rownames(out),log(out$Biomass),xlab="year", ylab="log Biomass", main=sp)
    
    plot(rownames(out), out$FishingLoss, xlab="year", ylab="F rate", main=sp)
  
    plot(rownames(out), out$FoodGain - out$UnAssimLoss - out$ActiveRespLoss, xlab="year", 
         ylab="Net feeding rate", main=sp )
  }
  
  par(mfrow=c(1,1))
  
  

  # ---------------------------------------------------------------------------- #
  # flow matrix ####
  # ---------------------------------------------------------------------------- #
  
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
  vlist <- data.frame(from_species,to_species,vul) # human readable list of v
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
    