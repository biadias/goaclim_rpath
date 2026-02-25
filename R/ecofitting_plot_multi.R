#------------------------------------------------------------------------------#
#AUTHORS: Bia Dias (Original code by Kerim Aydin)
#AFFILIATIONS: Alaska Fisheries Science Center / CICOES University of Washington
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#
# This function creates the scenario comparison plots for biomass fitting for all the species.
#------------------------------------------------------------------------------#

################################################################################
rsim.plot.fitbio.small.bd <- function(scene, run_list, species, datasource, cols) {
  bio.obj <- rsim.fit.obj(scene, run_list[[1]])$Biomass
  qdat <- bio.obj[bio.obj$Group == species &
                    bio.obj$Source %in% datasource, ]
  par(mar = c(3, 2, 2, 1))
  mn   <- qdat$obs          / qdat$survey_q
  sdlog <- sqrt(log(1.0 + (qdat$sd / qdat$survey_q) * (qdat$sd / qdat$survey_q) / (mn * mn)))
  up <- mn * exp(1.96 * sdlog)
  dn <- mn / exp(1.96 * sdlog)
  #up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
  #dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q
  #dealing with NAs from species that don't have time series
  all_ests <- unlist(lapply(run_list, function(r) r$annual_Biomass[, species]))
  val_pool <- c(all_ests, up)
  val_pool <- val_pool[!is.na(val_pool) & is.finite(val_pool)]

  
  if (length(val_pool) == 0) {
    max_y <- 1 
 } else {
   max_y <- max(val_pool)
 }
  
  
#  max_y <- max(up, na.rm=TRUE)
#  for(r in 1:length(run_list)){
#    max_y <- max(c(max_y,run_list[[r]]$annual_Biomass[, species]), na.rm=TRUE)
#  }
  
  for(r in 1:length(run_list)) {
    est <- run_list[[r]]$annual_Biomass[, species] # * qdat$survey_q
    years <- as.numeric(rownames(run_list[[r]]$annual_Biomass))
    if(r==1){
      plot(
        years,
        est,
        type = "l",
        col= cols[r],
        ylim = c(0, max_y),
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        bty = "n",
        lwd=2)
      }else{
        lines(years,
      est,
      col = cols[r], lwd=2)
      }
  }
  
  
  #est <- run$annual_Biomass[, species] # * qdat$survey_q
  #tot  <- sum(qdat$fit * qdat$wt)
  #plot(as.numeric(rownames(run$annual_Biomass)),est,type="l",
  #     ylim=c(0,max(up[!is.na(up)],est[!is.na(est)])),xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
    axis(1,
         mgp = c(3, 0.0, 0),
         tck = -0.04,
         cex.axis = 0.6)
    axis(2,
         mgp = c(3, 0.2, 0),
         tck = -0.04,
         cex.axis = 0.8)
  
  #mtext(paste(datasource,species),     side=1, cex=0.6, line=0.8, adj=0)
  multiline_text <- paste(species, ifelse(is.na(datasource), "", datasource), sep="\n")
  mtext(
    paste(multiline_text),
    side = 1,
    cex = 0.6,
    line = 1.8, 
    adj = 0
  )
  if (nrow(qdat) > 0) {
    tot  <- sum(qdat$fit * qdat$wt)
    #mtext(
    #  sprintf(
    #    "nl:%.3g  q:%.3g (%s)",
    #    tot,
    #    qdat$survey_q
    #  ),
    #  side = 1,
    #  line = 1.8,
    #  cex = 0.6,
    #  adj = 0
    #)
    points(as.numeric(qdat$Year), mn, pch=16, cex=0.9, col="grey60")
    segments(as.numeric(qdat$Year), y0 = up, y1 = dn, col="grey60")
    #mtext(sprintf("  q: %.3g", ), side=1, line=3, cex=0.8)
  }
  abline(h = scene$params$B_BaseRef[species],
         col = "#800000",
         lty = 4)
  
}

################################################################################
##### modified by BD 1/6/2026 ####
rsim.runplot.multi <- function(scene, run_list, species, cols = NULL, run_names=NULL) {
  
  if (is.null(cols)) {
    cols <- c("#000000", "#e69f00", "#56b4e9", "#009e73", "#0072b2", "#d55e00","#cc79a7", "#f0e442")
  }
  if (!is.list(run_list[[1]])) {
    run_list <- list(run_list)
  } 
  if (is.null(run_names)) {
    run_names <- paste("Run", 1:length(run_list))
  }
  legend_labels <- sapply(1:length(run_list), function(i){
    total_nll <- round(rsim.fit.obj(scene, run_list[[i]], FALSE),2)
  paste0(run_names[i], " (nll: ", total_nll, ")")
  })
  #run.out   <- rsim.run(scene, method='AB', years=years)
  
  #run.out <- run #rsim.fit.run(values,groups,vartypes,scene,"AB",years,T)

  
  #all_sources <- rsim.fit.list.bio.series(scene_fit)$sources
  all_series <- strsplit(rsim.fit.list.bio.series(scene)$all, ":")
  obs_species <- sapply(all_series, "[[", 1)
  obs_sources <- sapply(all_series, "[[", 2)
  
 #relevant_sources <- obs_sources[obs_species %in% species]
 #total_plots <- length(species) + length(relevant_sources[relevant_sources != ""])
 if (length(species) <= 4) {
   my_grid <- c(2, 2)
 } else if (length(species) <= 12) {
   my_grid <- c(3, 4)
   my_cex <- 0.8     # Medium size
   my_mar <- c(3, 3, 2, 1)
 } else {
   my_grid <- c(9, 11) # Fallback to your original dense grid
 }
  
  
  #fup()
  X11()
  #this_survey <- all_sources #"race_wgoa"

  par(mfrow = my_grid)
  #par(mfrow = c(9, 11))#, oma = c(0, 0, 3, 0)) #oma= outside of margin area for the text
  for (i in 1:length(species)) {
    sp_name <- species[i]
    sources <- obs_sources[obs_species == sp_name]
    if (length(sources) == 0) {
      rsim.plot.fitbio.small.bd(scene, run_list, species[i], NA, cols)
    } else{
      for (j in 1:length(sources)) {
        rsim.plot.fitbio.small.bd(scene, run_list, species[i], sources[j], cols)
      }
    }
  }
  
  #mtext("Multi-Run Biomass Comparison", outer=TRUE, side=3, line=1, cex=1.2)
  #return(run.out)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #plot.new()
  legend("bottomright", 
         legend = legend_labels, 
         col = cols[1:length(run_list)], 
         lty = 1, 
         lwd = 4, 
         inset=c(-0.2,0),
         #horiz = TRUE, 
         cex = 1.2, #1.2
         bty = "n")
}