################################################################################

render.show.fit <- function(bal, scene_fit, run_fit, output="test"){
  wdir <- paste(getwd(), "/html", sep="")
  #print(wdir)
  #print(output)
  
  rmarkdown::render("fit_display.Rmd", output_file=output, output_dir=wdir,  
                    params=list(bal=bal, scene_fit=scene_fit, run_fit=run_fit)
  )
  viewer <- getOption("viewer")
  viewer(paste(wdir,"/",output,".html",sep=''))
}

################################################################################
##### modified by BD 1/6/2026 ####
rsim.runplot <- function(scene, run, species, scene_name=NULL){
  if (is.null(scene_name)) {
    display_name <- " " 
  } else {
    display_name <- scene_name
  }
  #run.out   <- rsim.run(scene, method='AB', years=years)
  run.out <- run #rsim.fit.run(values,groups,vartypes,scene,"AB",years,T)
  
  #all_sources <- rsim.fit.list.bio.series(scene_fit)$sources
  all_series <- strsplit(rsim.fit.list.bio.series(scene)$all,":")
  obs_species <- sapply(all_series, "[[", 1)
  obs_sources <- sapply(all_series, "[[", 2)
  #fup()
  X11()
  #this_survey <- all_sources #"race_wgoa"  
  par(mfrow=c(9,11), oma=c(0,0,3,0)) #oma= outside of margin area
  for (i in 1:length(species)){
    sources <- obs_sources[obs_species==species[i]]
    if (length(sources)==0){    
      rsim.plot.fitbio.small(scene, run.out, species[i], NA)
    } else{
      for(j in 1:length(sources)){
        rsim.plot.fitbio.small(scene, run.out, species[i], sources[j])
      }
    }
  }
  mtext(paste("Fitting results for:", display_name, paste(", nll:", round(rsim.fit.obj(scene, run, FALSE),2))), outer=TRUE, side=3, line=1, cex=1.5, font=1)

  X11()
  par(mfrow=c(9,11), oma=c(0,0,3,0))
  for (sp in species){
    rsim.plot.fitcatch.small(scene, run.out, sp)
  }
  mtext(paste("Catch fitting results for:", display_name,paste(", nll:", round(rsim.fit.obj(scene, run, FALSE),2))), outer=TRUE, side=3, line=1, cex=1.5, font=1)
 #return(run.out)
  
}

################################################################################
rsim.plot.fitbio <- function(scene, run, species, datasource){
  bio.obj <- rsim.fit.obj(scene,run)$Biomass
  qdat <- bio.obj[bio.obj$Group==species & bio.obj$Source==datasource,]
  par(mar=c(3,2,2,1))
  mn   <- qdat$obs          / qdat$survey_q 
  sdlog <- sqrt(log(1.0 + (qdat$sd/qdat$survey_q) * (qdat$sd/qdat$survey_q) /(mn*mn)))
  up <- mn*exp(1.96*sdlog)
  dn <- mn/exp(1.96*sdlog)
  #up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
  #dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q 
  
  est <- run$annual_Biomass[,species] # * qdat$survey_q 
  tot  <- sum(qdat$fit * qdat$wt)
  
  #par(mai=c(0.1,0.1,0.1,0.1))
  # FIX NA issue in plot
  plot(as.numeric(rownames(run$annual_Biomass)),est,type="l",
       ylim=c(0,max(up[!is.na(up)],est[!is.na(est)])),xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
  axis(1,mgp=c(3, 0.0, 0),tck=-0.04, cex.axis=0.6)
  axis(2,mgp=c(3, 0.2, 0),tck=-0.04, cex.axis=0.8)
  #title(paste(species,"biomass",sprintf("   nll: %.3g",tot)), cex.main=0.9, line=0.5, adj=0)
  mtext(paste(datasource,species,"biomass",sprintf("   nll: %.3g",tot)), side=1, cex.main=0.9, line=0.8, adj=0)
  if (nrow(qdat)>0){
    mtext(sprintf("%s  q: %.3g", unique(qdat$Type), qdat$survey_q), side=1, line=1.8, cex=0.9, adj=0)
    #mtext(sprintf("  q: %.3g", ), side=1, line=3, cex=0.8)
  }
  abline(h=scene$params$B_BaseRef[species],col="red",lty=3)
  points(as.numeric(qdat$Year),mn)
  segments(as.numeric(qdat$Year),y0=up,y1=dn)  
}
################################################################################
################################################################################
rsim.plot.fitbio.small <- function(scene, run, species, datasource){
  bio.obj <- rsim.fit.obj(scene,run)$Biomass
  qdat <- bio.obj[bio.obj$Group==species & bio.obj$Source %in% datasource,]
  par(mar=c(3,2,2,1))
  mn   <- qdat$obs          / qdat$survey_q 
  sdlog <- sqrt(log(1.0 + (qdat$sd/qdat$survey_q) * (qdat$sd/qdat$survey_q) /(mn*mn)))
  up <- mn*exp(1.96*sdlog)
  dn <- mn/exp(1.96*sdlog)
  #up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
  #dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q 
  
  est <- run$annual_Biomass[,species] # * qdat$survey_q 
  tot  <- sum(qdat$fit * qdat$wt)
  
  #par(mai=c(0.1,0.1,0.1,0.1))
  # FIX NA issue in plot
  plot(as.numeric(rownames(run$annual_Biomass)),est,type="l",
       ylim=c(0,max(up[!is.na(up)],est[!is.na(est)])),xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
  axis(1,mgp=c(3, 0.0, 0),tck=-0.04, cex.axis=0.6)
  axis(2,mgp=c(3, 0.2, 0),tck=-0.04, cex.axis=0.8)
  #mtext(paste(datasource,species),     side=1, cex=0.6, line=0.8, adj=0)
  mtext(paste(species, datasource),     side=1, cex=0.6, line=0.8, adj=0)
  if (nrow(qdat)>0){
    mtext(sprintf("nl:%.3g  q:%.3g (%s)", tot, qdat$survey_q, unique(qdat$Type)), side=1, line=1.8, cex=0.6, adj=0)
    #mtext(sprintf("  q: %.3g", ), side=1, line=3, cex=0.8)
  }
  abline(h=scene$params$B_BaseRef[species],col="red",lty=3)
  points(as.numeric(qdat$Year),mn)
  segments(as.numeric(qdat$Year),y0=up,y1=dn)  
}
#######################################
rsim.plot.fitcatch.small <- function(scene, run, species){
catch.obj <- rsim.fit.obj(scene,run)$Catch
qdat <- catch.obj[catch.obj$Group==species,]
par(mar=c(3,2,2,1))
mn   <- qdat$obs
sdlog <- sqrt(log(1.0 + (qdat$sd*qdat$sd)/(mn*mn)))
up <- mn*exp(1.96*sdlog)
dn <- mn/exp(1.96*sdlog)
#up   <- mn + 1.96*qdat$sd
#dn   <- mn - 1.96*qdat$sd
est  <- run$annual_Catch[,species]
tot <- sum(qdat$fit * qdat$wt)
plot(as.numeric(rownames(run$annual_Catch)),est,type="l",
     ylim=c(0,max(up[!is.na(up)],est[!is.na(est)])),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
axis(1,mgp=c(3, 0.0, 0),tck=-0.04, cex.axis=0.6)
axis(2,mgp=c(3, 0.2, 0),tck=-0.04, cex.axis=0.8)
#title(paste(species,"catch",sprintf("   nll: %.3g",tot)), cex.main=0.9, line=0.5, adj=0)
mtext(paste(species,sprintf("   nll: %.3g",tot)), side=1, cex=0.6, line=0.8, adj=0)
points(as.numeric(qdat$Year),mn)
segments(as.numeric(qdat$Year),y0=up,y1=dn)
slist <- scene$params$spname[scene$params$FishFrom+1]
tcatch <- sum((scene$params$FishQ * scene$params$B_BaseRef[slist])[slist==species])
abline(h=tcatch,col="red",lty=3)
}
#################################################################################
rsim.plot.full <- function(scene, run, species){

  #layout(matrix(c(1,1,2,3), 2, 2, byrow = T)) 

# Biomass plotting  
# Remove TMP here - don't forget.
  bio.obj <- rsim.fit.obj(scene,run)$Biomass
  sdat <- bio.obj[bio.obj$Group==species,]
  #sid <- paste(sdat$Source, sdat$Group, sep=":")
  sidset <- unique(sdat$Source)
  
  par(mfrow=c(1,length(sidset)+1))
  par(oma=c(0.0,0.0,0.0,0.0))
  par(mar=c(3,2,2,1))  
  
for (S in sidset){
  qdat <- sdat[sdat$Source==S,] 
  #survey_q <- 1
  mn   <- qdat$obs          / qdat$survey_q 
  sdlog <- sqrt(log(1.0 + (qdat$sd/qdat$survey_q) * (qdat$sd/qdat$survey_q) /(mn*mn)))
  up <- mn*exp(1.96*sdlog)
  dn <- mn/exp(1.96*sdlog)
  #up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
  #dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q 
  
  est <- run$annual_Biomass[,species] # * qdat$survey_q 
  tot  <- sum(qdat$fit * qdat$wt)
  
  #par(mai=c(0.1,0.1,0.1,0.1))
# FIX NA issue in plot
  plot(as.numeric(rownames(run$annual_Biomass)),est,type="l",
       ylim=c(0,max(up[!is.na(up)],est[!is.na(est)])),xaxt="n",yaxt="n",xlab="",ylab="", bty="n")
  axis(1,mgp=c(3, 0.0, 0),tck=-0.04, cex.axis=0.6)
  axis(2,mgp=c(3, 0.2, 0),tck=-0.04, cex.axis=0.8)
  #title(paste(species,"biomass",sprintf("   nll: %.3g",tot)), cex.main=0.9, line=0.5, adj=0)
  mtext(paste(S,species,"biomass",sprintf("   nll: %.3g",tot)), side=1, cex.main=0.9, line=0.8, adj=0)
  if (nrow(qdat)>0){
    mtext(sprintf("%s  q: %.3g", unique(qdat$Type), qdat$survey_q), side=1, line=1.8, cex=0.9, adj=0)
    #mtext(sprintf("  q: %.3g", ), side=1, line=3, cex=0.8)
  }
  abline(h=scene$params$B_BaseRef[species],col="red",lty=3)
  points(as.numeric(qdat$Year),mn)
  segments(as.numeric(qdat$Year),y0=up,y1=dn)  
}
  
# Catch plotting  
# REMOVE TMP here
  catch.obj <- rsim.fit.obj(scene,run)$Catch
  qdat <- catch.obj[catch.obj$Group==species,]
  mn   <- qdat$obs
  sdlog <- sqrt(log(1.0 + (qdat$sd*qdat$sd)/(mn*mn)))
  up <- mn*exp(1.96*sdlog)
  dn <- mn/exp(1.96*sdlog)
  #up   <- mn + 1.96*qdat$sd
  #dn   <- mn - 1.96*qdat$sd
  est  <- run$annual_Catch[,species]
  tot <- sum(qdat$fit * qdat$wt)
  plot(as.numeric(rownames(run$annual_Catch)),est,type="l",
       ylim=c(0,max(up[!is.na(up)],est[!is.na(est)])),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  axis(1,mgp=c(3, 0.0, 0),tck=-0.04, cex.axis=0.6)
  axis(2,mgp=c(3, 0.2, 0),tck=-0.04, cex.axis=0.8)
  #title(paste(species,"catch",sprintf("   nll: %.3g",tot)), cex.main=0.9, line=0.5, adj=0)
  mtext(paste(species,"catch",sprintf("   nll: %.3g",tot)), side=1, cex.main=0.9, line=0.8, adj=0)
  points(as.numeric(qdat$Year),mn)
  segments(as.numeric(qdat$Year),y0=up,y1=dn)
  slist <- scene$params$spname[scene$params$FishFrom+1]
  tcatch <- sum((scene$params$FishQ * scene$params$B_BaseRef[slist])[slist==species])
  abline(h=tcatch,col="red",lty=3)
  

}
#################################################################################
rsim.plot.ylim <- function(Rsim.output, spname, indplot = F, ...){
  opar <- par(no.readonly = T)
  if(indplot == F){
    # KYA April 2020 this seems incorrect?
    #biomass <- Rsim.output$out_Biomass[, 2:ncol(Rsim.output$out_Biomass)]
    biomass <- Rsim.output$out_Biomass[, spname]
    n <- ncol(biomass)
    start.bio <- biomass[1, ]
    start.bio[which(start.bio == 0)] <- 1
    rel.bio <- matrix(NA, dim(biomass)[1], dim(biomass)[2])
    for(isp in 1:n) rel.bio[, isp] <- biomass[, isp] / start.bio[isp]
  }
  if(indplot == T){
    spnum <- which(Rsim.output$params$spname == spname)
    biomass <- Rsim.output$out_Biomass[, spnum]
    n <- 1
    rel.bio <- biomass / biomass[1]
  }
  
  ymax <- max(rel.bio) + 0.1 * max(rel.bio)
  ymin <- min(rel.bio) - 0.1 * min(rel.bio)
  ifelse(indplot, xmax <- length(biomass), xmax <- nrow(biomass))
  
  #Plot relative biomass
  opar <- par(mar = c(4, 6, 2, 0))
  
  #Create space for legend
  plot.new()
  l <- legend(0, 0, bty='n', spname, 
              plot=FALSE, fill = line.col, cex = 0.6)
  # calculate right margin width in ndc
  w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
  
  par(omd=c(0, 1-w, 0, 1))
  plot(0, 0, xlim = c(0, xmax), 
       axes = F, xlab = '', ylab = '', type = 'n', ...)
  axis(1)
  axis(2, las = T)
  box(lwd = 2)
  mtext(1, text = 'Months', line = 2.5, cex = 1.8)
  mtext(2, text = 'Relative Biomass', line = 3, cex = 1.8)
  
  line.col <- rainbow(n)
  for(i in 1:n){
    if(indplot == T) lines(rel.bio,      col = line.col[i], lwd = 3)
    if(indplot == F) lines(rel.bio[, i], col = line.col[i], lwd = 3)
  }
  
  legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
         spname, fill = line.col, cex = 0.6)
  
  par(opar)
}

################################################################################
#'@export
rsim.plot.catch <- function(scene, run, species){
  qdat <- scene$fitting$Catch[scene$fitting$Catch$Group==species,]
  mn   <- qdat$obs
  up   <- mn + 1.96*qdat$sd
  dn   <- mn - 1.96*qdat$sd 
  tot <- 0 #sum(qdat$fit)
  plot(as.numeric(rownames(run$annual_Catch)),run$annual_Catch[,species],type="l",
       ylim=c(0,max(up,run$annual_Catch[,species])),xlab=tot,ylab="")
  mtext(side=2, line=2.2, paste(species,"catch"), font=2, cex=1.0)
  points(as.numeric(qdat$Year),mn)
  segments(as.numeric(qdat$Year),y0=up,y1=dn)
}


################################################################################
#'@export
rsim.plot.biomass <- function(scene, run, species){
  bio.obj <- rsim.fit.obj(scene,run)$Biomass 
  qdat <- bio.obj[bio.obj$Group==species,]
  #survey_q <- 1
  mn   <- qdat$obs/qdat$survey_q #qdat$obs_scaled #* qdat$survey_q   #/survey_q
  up   <- mn + 1.96*qdat$sd / qdat$survey_q #/survey_q
  dn   <- mn - 1.96*qdat$sd / qdat$survey_q #/survey_q 
  tot  <- sum(qdat$fit)
  #par(mai=c(0.1,0.1,0.1,0.1))
  plot(as.numeric(rownames(run$annual_Biomass)),run$annual_Biomass[,species],type="l",
       ylim=c(0,max(up,run$annual_Biomass[,species])),xlab="",ylab="")
  mtext(side=2, line=2.2, paste(species,"biomass"), font=2, cex=0.8)
  #cat(species, tot, "\n");
  if (nrow(qdat)>0){
    mtext(sprintf("NLL: %.3g", tot),           side=1, line=2, cex=0.8)
    mtext(sprintf("  q: %.3g", qdat$survey_q), side=1, line=3, cex=0.8)
  }
  abline(h=scene$params$B_Initial[species],col="red",lty=3)
  points(as.numeric(qdat$Year),mn)
  segments(as.numeric(qdat$Year),y0=up,y1=dn)
}

