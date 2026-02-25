################################################################################
#'@useDynLib Rpath
#'@export
read.fitting.biomass <- function(SCENE, cdat){
  
  # TODO KYA 7/24/24: add intelligent warning if missing column (e.g. Type)
  
  # Base variables
  SIM   <- SCENE
  years <- as.numeric(row.names(SCENE$fishing$ForcedFRate))
  species <- SIM$params$spname
  
  
  #cdat  <-read.csv(filename)
  missing_sp <- unique(cdat$Group[which(!(cdat$Group %in% species))])
  if(length(missing_sp)>0){
    warning("Following species in BIOMASS fit data not in model, dropped: ",missing_sp)
  }
  
  # drop Groups not found in model
  cmdat <- cdat[cdat$Group %in% species,]
  # drop lines with NAs in Value or Stdev or Scale
  c0dat <- cmdat[!is.na(cmdat$Value) & !is.na(cmdat$Stdev) & !is.na(cmdat$Scale),]
  # drop lines not in model run years
  c1dat <- c0dat[c0dat$Year %in% years,]
  # drop lines where Value, StDev or Scale are <=0
  ccdat <- c1dat[c1dat$Value>0 & c1dat$Stdev>0 & c1dat$Scale>0,]
  # replace NAs in type with index (maybe with future warning)
  ccdat$Type[is.na(ccdat$Type)] <- "index"
  
  #type <- as.character(rep("absolute",length(ccdat$YEAR)))
  ccdat$Year  <- as.character(ccdat$Year)
  ccdat$Group <- as.character(ccdat$Group)
  ccdat$Source <- as.character(ccdat$Source)
  ccdat$ID    <- paste(as.character(ccdat$Group),as.character(ccdat$Source),sep=":")
  
  obs  <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(ccdat$Value),
                 as.numeric(ccdat$Value) * as.numeric(ccdat$Scale))   
  sd   <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(ccdat$Stdev),
                 as.numeric(ccdat$Stdev) * as.numeric(ccdat$Scale))   
  wt   <- rep(1,length(obs))
  initial_q   <- rep(1,length(obs))
  SIM$fitting$Biomass <- cbind(ccdat,obs,sd,initial_q,wt)
  
  return(SIM)
  
}

################################################################################

read.fitting.biomass.cv <- function(SCENE, cdat){
  
  # TODO KYA 7/24/24: add intelligent warning if missing column (e.g. Type)
  
  # Base variables
  SIM   <- SCENE
  years <- as.numeric(row.names(SCENE$fishing$ForcedFRate))
  species <- SIM$params$spname
  
  # check if CV column exists to prevent script errors
  if(!"CV" %in% colnames(cdat)) cdat$CV <- NA
  
  #cdat  <-read.csv(filename)
  missing_sp <- unique(cdat$Group[which(!(cdat$Group %in% species))])
  if(length(missing_sp)>0){
    warning("Following species in BIOMASS fit data not in model, dropped: ",missing_sp)
  }
  
  # drop Groups not found in model
  cmdat <- cdat[cdat$Group %in% species,]
  # drop lines with NAs in Value or Stdev or Scale or CV
  c0dat <- cmdat[!is.na(cmdat$Value) & !is.na(cmdat$Scale) & 
                   (!is.na(cmdat$Stdev) | !is.na(cmdat$CV)), ]
  # drop lines not in model run years
  c1dat <- c0dat[c0dat$Year %in% years,]
  # drop lines where Value, StDev or Scale are <=0
  ccdat <- c1dat[c1dat$Value>0 & c1dat$Scale!=0,]
  # replace NAs in type with index (maybe with future warning)
  ccdat$Type[is.na(ccdat$Type)] <- "index"
  
  #type <- as.character(rep("absolute",length(ccdat$YEAR)))
  ccdat$Year  <- as.character(ccdat$Year)
  ccdat$Group <- as.character(ccdat$Group)
  ccdat$Source <- as.character(ccdat$Source)
  ccdat$ID    <- paste(as.character(ccdat$Group),as.character(ccdat$Source),sep=":")
  
  # Select/Calculate SD
  # give priotity to Stdev; if NA or 0, use Value * CV
  use_cv <- (is.na(ccdat$Stdev) | ccdat$Stdev <= 0) & (!is.na(ccdat$CV) & ccdat$CV > 0)
  if(any(use_cv)){
    cv_sources <- unique(ccdat$ID[use_cv])
    message("Note: Calculating Stdev from CV for: ", paste(cv_sources, collapse=", "))
  }
  final_stdev <- ifelse(use_cv, 
                        ccdat$Value * ccdat$CV, 
                        ccdat$Stdev)
  
  obs  <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(ccdat$Value),
                 as.numeric(ccdat$Value) * as.numeric(ccdat$Scale))   
  sd   <- ifelse(as.numeric(ccdat$Scale)<0, as.numeric(final_stdev),
                 as.numeric(final_stdev) * as.numeric(ccdat$Scale)) 
  # adding a safety net to find if any sd are NA or <=0 after calculation
  valid_indices <- !is.na(sd) & sd > 0
  ccdat <- ccdat[valid_indices, ]
  obs   <- obs[valid_indices]
  sd    <- sd[valid_indices]
  
  wt   <- rep(1,length(obs))
  initial_q   <- rep(1,length(obs))
  SIM$fitting$Biomass <- cbind(ccdat,obs,sd,initial_q,wt)
  
  return(SIM)
  
}

################################################################################


rsim.fit.remove.bio.timeseries <- function(scene, timeseries){
  scene$fitting$Biomass <- scene$fitting$Biomass[!(scene$fitting$Biomass$ID %in% timeseries),]
  return(scene)
}
  
################################################################################
#'@export
read.fitting.catch <- function(SCENE, cdat){
  
  # TODO KYA 7/24/24: add intelligent warning if missing column 
  
  SIM   <- SCENE
  years <- as.numeric(row.names(SCENE$fishing$ForcedFRate))
  # Columns needed
  #  Group	Year	Value	SD	Scale   
  #cdat  <- read.csv(filename)
  missing_sp <- unique(cdat$Group[which(!(cdat$Group %in% SIM$params$spname))])
  if(length(missing_sp)>0){
    warning("Following species in CATCH fit data not in model, dropped: ",missing_sp)
  }
  # drop Groups not found in model
  cmdat <- cdat[cdat$Group %in% SIM$params$spname,]
  # drop lines with NAs in Value or with years outside scenario years range
  ccdat <- cmdat[!is.na(cmdat$Value) & cmdat$Year %in% years,] 
  ccdat$Year  <- as.character(ccdat$Year)
  ccdat$Group <- as.character(ccdat$Group) 
  obs  <- as.numeric(ccdat$Value) * as.numeric(ccdat$Scale)   
  sd   <- as.numeric(ccdat$Stdev) * as.numeric(ccdat$Scale)  
  wt   <- rep(1,length(obs))
  SIM$fitting$Catch <- cbind(ccdat,obs,sd,wt)
  #sdat  <- aggregate(as.numeric(ccdat$Value)*as.numeric(ccdat$Scale),list(ccdat$Year,ccdat$Group),"sum")
  #sd    <- 0.1*sdat$x
  #colnames(SIM$fitting$CATCH) <- c("year","species","obs","sd","wt")

  # Apply fit fishing to matrix
  #SIM$fishing$ForcedEffort[] <- 0
  #SIM$fishing$ForcedCatch[matrix(c(SIM$fitting$Catch$Year, SIM$fitting$Catch$Group),
  #                        length(SIM$fitting$Catch$Year),2)] <- SIM$fitting$Catch$obs
  return(SIM)
}

################################################################################
# Performs 3 steps:  Sets effort to 0, sets default F, sets catch
#'@export
fitcatch.to.forcecatch <- function(Rsim.scenario, Rpath){

# TODO: figure out detritus
  scene <- Rsim.scenario
  bal   <- Rpath
  
  # zero out all gear effort
  scene$fishing$ForcedEffort[] <- 0 
  
  # Set forced Frate to Ecopath F
  splist <- c(rpath.living(bal),rpath.detrital(bal))
  flist  <- (rowSums(bal$Landings) + rowSums(bal$Discards))/bal$Biomass
  for (sp in splist){
    scene$fishing$ForcedFRate[,sp] <- flist[sp]
  }
  
  #Put supplied catch in ForcedCatch, and zero Frate for those years
  catchmat <- matrix(c(scene$fitting$Catch$Year, scene$fitting$Catch$Group),
                     length(scene$fitting$Catch$Year),2)
  
  scene$fishing$ForcedCatch[catchmat] <- scene$fitting$Catch$obs
  scene$fishing$ForcedFRate[catchmat] <- 0
  
  return(scene)

}

#################################################################################
#
#'@export
rsim.fit.list.catch.series <- function(Rsim.scenario){
  
  scene <- Rsim.scenario
  out   <- list()
  out$groups  <- unique(scene$fitting$Catch$Group)
  return(out)  
  
}

#################################################################################
#
#'@export
rsim.fit.list.bio.bysurvey <- function(Rsim.scenario){
  
  scene <- Rsim.scenario
  out   <- list()
  sources <- unique(scene$fitting$Biomass$Source)
  for (i in 1:length(sources)){
   out[[i]] <- unique(scene$fitting$Biomass$Group[scene$fitting$Biomass$Source==sources[i]])
   names(out)[i] <- sources[i]
  }
  
  return(out)  
  
}

#################################################################################
#
#'@export
rsim.fit.list.bio.series <- function(Rsim.scenario){
  
  scene <- Rsim.scenario
  out   <- list()
  out$groups  <- unique(scene$fitting$Biomass$Group)
  out$sources <- unique(scene$fitting$Biomass$Source)
  out$all     <- unique(paste(scene$fitting$Biomass$Group,scene$fitting$Biomass$Source,sep=':'))
  return(out)  
  
}

#################################################################################
#
#'@export
rsim.fit.set.catch.wt <- function(Rsim.scenario, group, wt){
  
  scene  <- Rsim.scenario
  series <- scene$fitting$Catch$Group %in% group  
  
  if (sum(series)==0){
    warning("No catch fitting data found for ",group) 
    return(scene)
  }
  
  scene$fitting$Catch[series,"wt"] <- wt
  
  return(scene)
  
}

#################################################################################
#
#'@export
rsim.fit.set.bio.wt <- function(Rsim.scenario, group, source, wt){

  scene <- Rsim.scenario
  series <- scene$fitting$Biomass$Source %in% source & 
            scene$fitting$Biomass$Group %in% group  

  if (sum(series)==0){
    warning("No biomass fitting data found for ",source," ",group) 
    return(scene)
  }
  
  scene$fitting$Biomass[series,"wt"] <- wt
  
  return(scene)
  
}


#################################################################################
#
#'@export
rsim.fit.set.q <- function(Rsim.scenario, group, source, q=NULL, years=NULL, type=NULL){
  
  scene <- Rsim.scenario
  series <- scene$fitting$Biomass$Source %in% source & 
            scene$fitting$Biomass$Group %in% group
  if(is.null(type)){type="fixed"}
  
  if (sum(series)==0){
    warning("No biomass fitting data found for ",source," ",group) 
    return(scene)
  }
  
  # If both q and years are null, set initial_q to 1.0 and Type to "index"
  if (is.null(q) & is.null(years)){
    scene$fitting$Biomass[series,"initial_q"] <- 1.0
    scene$fitting$Biomass[series,"Type"]      <- "index"
    return(scene)
  }
  
  # If a non-NULL q is supplied, use that
  if (!is.null(q)){
    qq <- as.numeric(q)
    if (!is.na(qq) & qq>0){
      scene$fitting$Biomass[series,"initial_q"] <- qq
      scene$fitting$Biomass[series,"Type"]      <- type
      return(scene)
    } else {
      warning("Supplied q for ",source," ",group," is NA or non-positive - q unchanged.")
      return(scene)
    }
  }
      
  # At this point, use years
  lookup <- scene$fitting$Biomass$Source %in% source & 
            scene$fitting$Biomass$Group %in% group &  
            scene$fitting$Biomass$Year %in% years
  
  dat <- scene$fitting$Biomass[lookup,]
  est <- scene$params$B_BaseRef[group]
  obs <- mean(dat$Value * dat$Scale)
  qq  <- obs/est
  
  if (is.na(qq) | is.nan(qq) | qq<=0){
    warning(source," ",group," survey data in ", years, " is NA, NaN or non-positive - q unchanged.")
  } else {
    scene$fitting$Biomass[series,"initial_q"] <- qq
    scene$fitting$Biomass[series,"Type"]      <- type
  }
  
  return(scene)
  
}

#################################################################################
#'@export
rsim.fit.obj <- function(SIM, RES, verbose=TRUE){
  FLOGTWOPI <- 0.5*log(2*pi) #0.918938533204672
  epsilon <- 1e-36
  
  OBJ <- list()
  OBJ$tot <- 0
  
  # BIOMASS to NON-RESCALED "Actual" biomass estimate
  est <- RES$annual_Biomass[matrix(c(as.character(SIM$fitting$Biomass$Year),as.character(SIM$fitting$Biomass$Group)),
                                   ncol=2)] + epsilon
  obs <- SIM$fitting$Biomass$obs + epsilon
  sd  <- SIM$fitting$Biomass$sd  + epsilon
  wt  <- SIM$fitting$Biomass$wt
  initial_q <- SIM$fitting$Biomass$initial_q
  # Series id (sid) is Source and Group columns combined
  sid <- paste(SIM$fitting$Biomass$Source, SIM$fitting$Biomass$Group, sep=":")
  
  # We need to get variance-weighted survey means by species, for
  # calculating mean values needed for setting best-fit q
  
  # Formula for weighted average q: 
  # q = exp(sum(w * log(obs/est))/sum(w)) where w is wt/sd  
  logdiff       <- log(obs/est)
  sdlog         <- sqrt(log(1.0+sd*sd/(obs*obs))) # sigma^2 of lognormal dist 
  wt_sd_inverse <- wt/sdlog# sd
  wt_logdiffsum <- tapply(logdiff*wt_sd_inverse, sid ,sum)
  wt_sum        <- tapply(wt_sd_inverse,         sid ,sum)
  q_est         <- exp(wt_logdiffsum/wt_sum) # need ifelse here for 0 weights?
  survey_q      <- ifelse(SIM$fitting$Biomass$Type=="index", 
                          q_est[sid], initial_q)
  survey_q      <-ifelse(is.na(survey_q) | is.nan(survey_q),initial_q,survey_q)
  
  ## Jan 2023 incorrect code
  #inv_var <- 1.0/(sd*sd)
  #obs_sum <- tapply(obs*inv_var*wt, as.character(SIM$fitting$Biomass$Group),sum)
  #inv_sum <- tapply(inv_var*wt,     as.character(SIM$fitting$Biomass$Group),sum)
  #obs_mean <- obs_sum/inv_sum
  #est_mean <- tapply(est,as.character(SIM$fitting$Biomass$Group),mean)
  #survey_q <- ifelse(SIM$fitting$Biomass$Type=="absolute", 1.0,
  #            #(obs_mean/est_mean)[as.character(SIM$fitting$Biomass$Group)])
  #            (est_mean/obs_mean)[as.character(SIM$fitting$Biomass$Group)])
  #obs_scaled <-obs*survey_q 
  #sdlog  <- sqrt(log(1.0+sd*sd*survey_q*survey_q/(obs_scaled*obs_scaled)))
  sdiff  <- log((obs/survey_q)/est)/sdlog
  fit    <- wt * (FLOGTWOPI + log(sdlog) + 0.5*sdiff*sdiff)
  
  if (verbose){
    obs_scaled  <- obs/survey_q
    OBJ$Biomass <- cbind(SIM$fitting$Biomass,est,survey_q,obs_scaled,sdiff,fit)
  } else {
    OBJ$tot <- OBJ$tot + sum(fit)
  }
  
  # Catch compared (assumes all catch is clean, absolute values)
  est <- RES$annual_Catch[matrix(c(as.character(SIM$fitting$Catch$Year),as.character(SIM$fitting$Catch$Group)),
                                 ncol=2)] + epsilon
  obs <- SIM$fitting$Catch$obs + epsilon
  sd  <- SIM$fitting$Catch$sd  + epsilon
  sdlog  <- sqrt(log(1.0+sd*sd/(obs*obs)))
  sdiff  <- log(obs/est)/sdlog
  fit    <- SIM$fitting$Catch$wt * (log(sdlog) + FLOGTWOPI + 0.5*sdiff*sdiff)
  if (verbose){
    OBJ$Catch <- cbind(SIM$fitting$Catch,est,sdiff,fit)
  } else {
    OBJ$tot <- OBJ$tot + sum(fit)
  }
  
  # # RATION
  # obs <- SIM$fitting$ration$obs + epsilon
  # sd  <- SIM$fitting$ration$sd  + epsilon
  # inv_var <- (1.0/sd)*(1.0/sd)
  # obs_sum <- tapply(obs*inv_var,as.character(SIM$fitting$ration$Group),sum)
  # inv_sum <- tapply(inv_var,as.character(SIM$fitting$ration$Group),sum)
  # obs_mean <- obs_sum/inv_sum
  # est <- RES$annual_QB[matrix(c(as.character(SIM$fitting$ration$Year),as.character(SIM$fitting$ration$Group)),
  #                             ncol=2)] + epsilon
  # est_mean <- tapply(est,as.character(SIM$fitting$ration$Group),mean)
  # survey_q <- (obs_mean/est_mean)[as.character(SIM$fitting$ration$Group)]
  # est_scaled <-est*survey_q 
  # sdlog  <- sqrt(log(1.0+sd*sd/(obs*obs)))
  # sdiff  <- (log(obs)-log(est_scaled))/sdlog
  # fit    <- SIM$fitting$ration$wt * (log(sdlog) + FLOGTWOPI + 0.5*sdiff*sdiff)
  # OBJ$ration <- cbind(GOA_SIM$fitting$ration,est,survey_q,est_scaled,sdiff,fit)
  # 
  # # Diet proportions estimation
  # linklook   <- matrix(c(as.character(SIM$fitting$diets$Year),as.character(SIM$fitting$diets$simlink)),ncol=2)
  # totlook    <- matrix(c(as.character(SIM$fitting$diets$Year),as.character(SIM$fitting$diets$pred)),ncol=2) 
  # dietTot    <- tapply(RES$annual_Qlink[linklook],list(SIM$fitting$diets$Year,SIM$fitting$diets$pred),sum)
  # dietProp   <- RES$annual_Qlink[linklook]/dietTot[totlook]
  # logest     <- log(dietProp)
  # #NEGATIVE log likelihood now
  # fit        <- -SIM$fitting$diets$wt * (SIM$fitting$diets$log_diff + SIM$fitting$diets$alphaM1*logest)  
  # OBJ$diet   <- cbind(SIM$fitting$diets,dietProp,logest,fit)
  
  # Final summation and return
  if(verbose){
    OBJ$tot <- sum(OBJ$Biomass$fit, OBJ$Catch$fit)# , OBJ$ration$fit, OBJ$diet$fit)
    return(OBJ)
  }
  else{
    return(OBJ$tot)
  }
}

#################################################################################
#'@export
rsim.fit.table <- function(SIM,RES){
  fitobj  <- rsim.fit.obj(SIM,RES,verbose=T)
  Btmp <- tapply(fitobj$Biomass$fit,fitobj$Biomass$Group,sum)
  Ctmp <- tapply(fitobj$Catch$fit,fitobj$Catch$Group,sum)
  out <- rep(NA,length(SIM$params$spname)); names(out)<- SIM$params$spname
  Biomass <- out; Biomass[names(Btmp)] <- Btmp
  Catch <- out;   Catch[names(Ctmp)] <- Ctmp
  return(data.frame(Biomass,Catch))
}
#################################################################################
#'@export
rsim.fit.obj.species <- function(SIM,RES,species=NULL){
  OBJ <- list()
  fitobj <- rsim.fit.obj(SIM,RES,verbose=T)
  OBJ$Biomass <- fitobj$Biomass[fitobj$Biomass$Group%in%species,] 
  OBJ$Catch   <- fitobj$Catch[fitobj$Catch$Group%in%species,]
  return(OBJ)
}

#################################################################################
#Internal Only
rsim.fit.apply <- function(values, species, vartype, scene.params){

# Mzero
  mzerodiff <- values[vartype=="mzero"]
  mzero.sp  <- species[vartype=="mzero"]
  scene.params$MzeroMort[mzero.sp] <- mzerodiff

# Predator contribution to predprey vulnerability   
  predvuls <- values[vartype=="predvul"]
  names(predvuls) <- species[vartype=="predvul"]   
  preddiff <- as.numeric(predvuls[scene.params$spname[scene.params$PreyTo+1]])
  preddiff[is.na(preddiff)] <- 0
  
# Prey contribution to predprey vulnerability
  preyvuls <- values[vartype=="preyvul"]
  names(preyvuls) <- species[vartype=="preyvul"]   
  preydiff <- as.numeric(preyvuls[scene.params$spname[scene.params$PreyFrom+1]])
  preydiff[is.na(preydiff)] <- 0
  
# Apply pred and prey vuls above to actual pred/prey VV in scenario
  scene.params$VV <- (1 + exp(log(scene.params$VV-1) + preddiff + preydiff))
  
  return(scene.params)
}
#################################################################################
#'@export
rsim.fit.run <- function(values,
                         species,
                         vartype,
                         scene,
                         run_method,
                         run_years,
                         verbose = F,
                         penalty_weight = 0,
                         nll_details=FALSE,
                         ...) {
  if (!all(is.na(values))){
    scene$params <- rsim.fit.apply(values, species, vartype, scene$params)
  }
  run.out <- rsim.run(scene, method = run_method, years = run_years)
  if (!verbose) { #return(rsim.fit.obj(scene, run.out, FALSE))
    base_error <- rsim.fit.obj(scene, run.out, FALSE)
    penalty_term <- 0
    total_nll <- base_error
    if(!all(is.na(values))){
      penalty_term <- penalty_weight * sum(values^2)
      total_nll <- base_error + penalty_term
      
      if(nll_details){
        return(list(
          total_nll = total_nll,
          base_error = base_error,
          penalty_term = penalty_term
        ))
      }
      
    }
    return(total_nll)
  }
  else{
    return(run.out)
  }
}
#################################################################################
#'@export
rsim.fit.update <- function(values, species, vartype, scene){
  scene$params <- rsim.fit.apply(values, species, vartype, scene$params) 
  return(scene)
}

#################################################################################
rsim.predprey.table <- function(scene, pred=NULL, prey=NULL){
  pp_all <- data.frame(
    from = scene$params$spname[scene$params$PreyFrom+1],
    to   = scene$params$spname[scene$params$PreyTo+1],
    QQ   = scene$params$QQ,
    VV   = scene$params$VV,
    DD   = scene$params$DD)
  if (!is.null(pred)){pp_all <- pp_all[pp_all$to   %in% pred,]}
  if (!is.null(prey)){pp_all <- pp_all[pp_all$from %in% prey,]}
  return(pp_all)
}

#################################################################################
test<-function(){

#Group	Year	Value	SD	Scale

# DATA from CATCH time series (Angie provided)   


# APPLY FISHING TO FITTING
  SIM$fishing$EFFORT[]<-0
  colnames(SIM$fishing$CATCH)<-SIM$params$spname[1:(SIM$params$NUM_BIO+1)]
  #rownames(SIM$fishing$CATCH)<-c(years,end_year+1)
  SIM$fishing$CATCH[matrix(c(as.character(SIM$fitting$CATCH$year),as.character(SIM$fitting$CATCH$species)),
                          length(SIM$fitting$CATCH$year),2)] <- SIM$fitting$CATCH$obs 

# diet composition  
  dfiles <- c("data/HMC_GOA_pollockdiet.csv","data/HMC_GOA_coddiet.csv","data/HMC_GOA_atfdiet.csv","data/HMC_GOA_halibutdiet.csv")
  dcdat <- read_diet_alphas(SIM,dfiles)
  SIM$fitting$diets <- dcdat[dcdat$year %in% years,]
# total ration index
  qdat <- NULL
  for (f in dfiles){
    ddat <- read.csv(f)
    qdat <- rbind(qdat,data.frame(ddat$year,ddat$pred,ddat$cperwMean,ddat$cperwSD,rep(1,length(ddat[,1]))))
  }
  colnames(qdat)<-c("year","species","obs","sd","wt")
  SIM$fitting$ration<-qdat[qdat$year %in% years,]


}

