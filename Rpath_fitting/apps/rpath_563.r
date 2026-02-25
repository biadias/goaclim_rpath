#####################################################################################
#'@useDynLib Rpath
#'@export
runRpathShiny <- function(lab) {
  # locate all the shiny app examples that exist
  validLabs <- list.files(system.file("labs", package = "Rpath"))
  validLabsMsg <-paste0("Valid Labs: '",paste(validLabs, collapse = "', '"),"'")

  # if an invalid example is given, throw an error
  if (missing(lab) || !nzchar(lab) ||
      !lab %in% validLabs) {stop(validLabsMsg, call. = F)}

  # find and launch the app
  appDir <- system.file("labs", lab, package = "Rpath")
  shiny::runApp(appDir, display.mode = "normal")
}

#####################################################################################
#'@export
rsim.deriv.q <- function(Rpath.scenario, year=0, month=0, tstep=0){
  scene <- copy(Rpath.scenario)
  rout <- deriv_vector(scene$params,  scene$start_state, 
                       scene$forcing, scene$fishing,
                       scene$stanzas, year, month, tstep)
  
  rtab <- data.frame(scene$params$spname[scene$params$PreyTo+1], 
                     scene$params$spname[scene$params$PreyFrom+1],
                     rout$Qlink)
  colnames(rtab)<-c("Predator","Prey","Q")
  return(rtab)
}

#####################################################################################
#'@export
get.rsim.predprey <- function(scene,predator="all",prey="all"){  
  pd <- scene$params$spname[scene$params$PreyTo+1]
  py <- scene$params$spname[scene$params$PreyFrom+1]
  if (predator=="all") predator <- scene$params$spname
  if (prey=="all")     prey     <- scene$params$spname
  link <- which((pd %in% predator) & (py %in% prey))
  predname <- pd[link]
  preyname <- py[link]
  QQ <- scene$params$QQ[link]
  DD <- scene$params$DD[link]
  VV <- scene$params$VV[link]
  PredPredWeight <- scene$params$PredPredWeight[link]
  PreyPreyWeight <- scene$params$PreyPreyWeight[link]
  HandleSwitch <- scene$params$HandleSwitch[link]
  return(data.frame(link,predname,preyname,QQ,VV,DD,
                    HandleSwitch,PredPredWeight,PreyPreyWeight))
}  

#####################################################################################
#'@export
get.rsim.predprey.link <- function(scene,predator="all",prey="all"){  
  pd <- scene$params$spname[scene$params$PreyTo+1]
  py <- scene$params$spname[scene$params$PreyFrom+1]
  if (predator=="all") predator <- scene$params$spname
  if (prey=="all")     prey     <- scene$params$spname
  link <- which((pd %in% predator) & (py %in% prey))
  return(link)
}  

#####################################################################################
#'@export
adjust.rsim.predprey <- function(scene,parameter,
                                 predator="all",prey="all",values){  
  pd <- scene$params$spname[scene$params$PreyTo+1]
  py <- scene$params$spname[scene$params$PreyFrom+1]
  if (predator=="all") predator <- scene$params$spname
  if (prey=="all")     prey     <- scene$params$spname
  link <- which((pd %in% predator) & (py %in% prey))
  
  if (parameter %in% c('QQ', 'VV', 'DD', 'HandleSwitch', 'PredPredWeight', 
                       'PreyPreyWeight')){
    scene$params[[parameter]][link] <- values
  }
  else {(warning(parameter," is not a valid predator/prey parameter."))}
  
  return(scene)
}   

#################################################################
#'@export
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

