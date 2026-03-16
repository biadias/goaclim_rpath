
#library(dplyr)
#library(Rpath)
#source("../R/xml_convert.r")

# The '-Basic estimates.csv' file is created for a model in EwE as follows:
#   1. Load the .eiixml (or access database) file into EwE.
#   2. On the 'Model parameters' tab, set 'General options -> Relevant decimal digits' to 8. 
#   3. On the 'Output->Basic estimates' tab, click export icon in upper right corner. 

# simple color output function  
  catcol <- function(txt, col=31){cat(paste0("\033[0;", col, "m",txt,"\033[0m","\n"))}  

# Function assumes model name is a "directory/base name" combo that points
# where the directory/base name includes a -Basic estimates.csv and an .eiixml
rpath.ewe.comptable <- function(model_name, dir_name){   
  #csv <- "Tampa Bay-Basic estimates.csv"
  eiifile     <- paste(dir_name,"/",model_name, ".eiixml", sep="")
  csvfile     <- paste(dir_name,"/",model_name, "-Basic estimates.csv", sep="")
  unbal <- Rpath::create.rpath.from.eiixml(eiifile)
  #if(unbal$stanzas$NStanzaGroups>0){
    unbal <- rpath.stanzas(unbal)
  #}
  Rpath::check.rpath.params(unbal)
  #bal <- Rpath::rpath(unbal) 
  bal <- rpath(unbal)
  csv.out <- read.csv(csvfile)
  # If there's no stanzas, the total mortality column is missing from CSV - add a dummy
    if(length(csv.out)<13){csv.out <- data.frame(append(csv.out, list(tm=NA), after=6))}
  names(csv.out)<-  c("X", "Group.name", "Trophic.level", "Hab.area",
                    "Biomass.in.habitat.area", "Biomass", "Total.mortality",
                    "Production.biomass", "Consumption.biomass",   
                    "Ecotrophic.Efficiency", "Production.consumption",
                    "Biomass.accumulation", "BA.rate")  
  
  # drop placeholder lines that head stanza groups
  csv.out <- csv.out[!is.na(csv.out$X),]
  
  # Clean names to rpath standard
  csv.out$group <- janitor::make_clean_names(csv.out$Group.name)
  
  # data frame of rpath output, removing fisheries
  rpath.dat <- data.frame(group=bal$Group,type=bal$type, tl=bal$TL, biomass=bal$Biomass,
                          pb=bal$PB,qb=bal$QB, EE=bal$EE, pc=bal$GE, 
                          ba=bal$BA)[c(rpath.living(bal), rpath.detrital(bal)),]
  
  # Joining rpath and ewe
  rpath_ewe_table <- rpath.dat %>% 
    dplyr::left_join(csv.out, by="group") %>%
    dplyr::mutate(
      tl_test      = abs(tl-Trophic.level)/Trophic.level, 
      biomass_test = abs(biomass-Biomass)/Biomass, 
      pb_test      = abs(pb-ifelse(is.na(Production.biomass),Total.mortality,Production.biomass))/
                            ifelse(is.na(Production.biomass),Total.mortality,Production.biomass),       
      qb_test      = abs(qb-Consumption.biomass)/Consumption.biomass,
      ee_test      = abs(EE-Ecotrophic.Efficiency),#/Ecotrophic.Efficiency,
      pc_test      = abs(pc-Production.consumption)/Production.consumption,
      ba_test      = abs(ba-Biomass.accumulation)/Biomass.accumulation
  )
  row.names(rpath_ewe_table) <- rpath_ewe_table$group
  return(rpath_ewe_table)
}  

rpath.ewe.check.table <- function(rp.ewe){
  parameter <- c("TL","B","PB","QB","EE","PC") 
  max_group <- max_diff <- flag <- rep(NA,length(parameter))
  vals <- rp.ewe %>%
    select(tl_test, biomass_test, pb_test, qb_test, ee_test, pc_test)
  for (v in 1:length(vals)){
    dat <- vals[,v]
    ind <- which.max(dat)
    if(length(ind)>0){
      max_group[v] <- rp.ewe$group[which.max(dat)]
      max_diff[v] <- max(dat, na.rm=T)
      flag[v]=ifelse((max_diff[v]>=0.01),"***",ifelse(max_diff[v]>=0.0001,"*","-"))
    }
    else{
      max_group[v] <- "UNKNOWN"
      max_diff[v] <- -1
      flag[v] <- "***"
    }
    
  }
  return(data.frame(parameter,max_group,max_diff,flag))
}

# # Single use test
#   model_name <- "Tampa Bay"
#   cat("\n",model_name,'\n')
#   rpath_ewe <- rpath.ewe.comptable(model_name,"xml_examples")
#   rpath.ewe.checks(rpath_ewe)
#   
# # Get a list of all files in the xml_examples directory that have
# # a '-Basic estimates.csv' output present.  It's assumed there's
# # a .eiixml file of the same name in the same directory.
# ewe_csvs <- list.files("xml_examples", pattern="-Basic estimates.csv")
# model_names <- stringr::str_remove_all(ewe_csvs,"-Basic estimates.csv")  
# 
# #broken_models <- c("Antarctic", "Barents Sea", "Denmark, Faroe Islands")
# #models_filtered <- model_names[!(model_names %in% broken_models)]
# 
# models_filtered <- model_names
# 
# tests<-list()
# for (m in models_filtered){
#   tryCatch({
#     cat("\n",m,'\n')
#     tests[[m]] <- rpath.ewe.comptable(m, "xml_examples")
#     rpath.ewe.checks(tests[[m]])
#   },error=function(e){catcol(paste("ERROR :",conditionMessage(e)))})    
# }  

rpathTMP <- function(Rpath.params, eco.name = NA, eco.area = 1) {
  #Need to define variables to eliminate check() note about no visible binding
  Type <- Group <- DetInput <- ProdCons <- PB <- QB <- noB <- noEE <- alive <- NULL
  BEE <- Biomass <- Q <- BioAcc <- BioQB <- diag.a <- EEa <- B <- M0 <- NULL
  QBloss <- Unassim <- Ex <- NULL
  
  # Model Parameters - Basic parameters, detritus fate, catch, discards in that order
  model <- copy(Rpath.params$model)
  
  #Diet Parameters - diet matrix, predators as columns, prey as rows - include
  #producers as predators even though they do not consume any groups
  diet <- copy(Rpath.params$diet)
  
  #Check that all columns of model are numeric and not logical
  if(length(which(sapply(model, class) == 'logical')) > 0){
    logic.col <- which(sapply(model, class) == 'logical')
    for(i in 1:length(logic.col)){
      set(model, j = logic.col[i], value = as.numeric(model[[logic.col[i]]]))
    }
  }
  
  #Remove first column if names (factor or character)
  if(sapply(diet, class)[1] == 'factor')    diet[, 1 := NULL]
  if(sapply(diet, class)[1] == 'character') diet[, 1 := NULL]
  
  #Adjust diet comp of mixotrophs
  mixotrophs <- which(model[, Type] > 0 & model[, Type] < 1)
  mix.Q <- 1 - model[mixotrophs, Type]
  for(i in seq_along(mixotrophs)){
    new.dc <- diet[, mixotrophs[i], with = F] * mix.Q[i]
    diet[, mixotrophs[i] := new.dc]
  }
  
  #Convert NAs to zero in diet matrix
  diet[is.na(diet)] <- 0
  
  # Get number of groups, living, dead, and gear
  ngroups <- nrow(model)
  nliving <- nrow(model[Type <  2, ])
  ndead   <- nrow(model[Type == 2, ])
  ngear   <- nrow(model[Type == 3, ])
  
  nodetrdiet <- diet[1:nliving, ]
  model[is.na(DetInput), DetInput := 0]
  
  # fill in GE(PQ), QB, or PB from other inputs
  GE   <- ifelse(is.na(model[, ProdCons]), model[, PB / QB],       model[, ProdCons])
  QB.1 <- ifelse(is.na(model[, QB]),       model[, PB / GE],       model[, QB])
  PB.1 <- ifelse(is.na(model[, PB]),       model[, ProdCons * QB], model[, PB])
  model[, QB := QB.1]
  model[, PB := PB.1]
  
  # define landings, discards, necessary sums
  landmat     <- model[, (10 + ndead + 1):(10 + ndead + ngear), with = F]
  discardmat  <- model[, (10 + ndead + 1 + ngear):(10 + ndead + (2 * ngear)), with = F]
  totcatchmat <- landmat + discardmat
  
  if (is.data.frame(totcatchmat)){
    totcatch <- rowSums(totcatchmat)
    landings <- rowSums(landmat)    
    discards <- rowSums(discardmat)  
    gearland <- colSums(landmat,   na.rm = T)
    geardisc <- colSums(discardmat, na.rm = T)
  }else{
    totcatch <- totcatchmat
    landings <- landmat    
    discards <- discardmat 
    gearland <- sum(landmat,    na.rm = T)
    geardisc <- sum(discardmat, na.rm = T)                     
  }   
  
  geartot <- gearland + geardisc
  model[, landings := landings]
  model[, discards := discards]
  model[, totcatch := totcatch]
  
  # flag missing pars and subset for estimation
  model[, noB   := 0]
  model[, noEE  := 0]
  model[, alive := 0]
  model[, BEE   := 0]
  model[is.na(Biomass), noB   := 1]
  model[is.na(EE),      noEE  := 1]
  model[Type < 2,       alive := 1]
  model[noB == 0 & noEE == 0, BEE := 1]
  
  # define detritus fate matrix
  detfate <- model[, (10 + 1):(10 + ndead), with = F]
  detdetfate <- model[Type==2, (10 + 1):(10 + ndead), with = F]
  
  # set up and solve the system of equations for living group B or EE
  living  <- model[alive == 1, ]
  
  #Set up right hand side b
  living[, Ex := totcatch + BioAcc]
  living[, BioQB := Biomass * QB]
  cons  <- as.matrix(nodetrdiet) * living$BioQB[col(as.matrix(nodetrdiet))]
  living[, b := Ex + rowSums(cons, na.rm = T)] 
  
  #Set up A matrix
  living[noEE == 1, diag.a := Biomass * PB]
  living[noEE == 0, diag.a := PB * EE]
  
  #Special case where B and EE are known then need to solve for BA
  #living[BEE == 1, b := b - (Biomass * PB * EE)]
  #living[BEE  == 1, diag.a := 0] #Need to work on this solution
  
  A       <- matrix(0, nliving, nliving)
  diag(A) <- living[, diag.a]
  QBDC    <- as.matrix(nodetrdiet) * living$QB[col(as.matrix(nodetrdiet))]
  dimnames(QBDC) <- list(NULL, NULL)
  QBDC[is.na(QBDC)] <- 0
  #Flip noB flag for known B and EE
  #living[BEE == 1, noB := 1]
  QBDCa <- as.matrix(QBDC) * living$noB[col(as.matrix(QBDC))]
  A     <- A - QBDCa 
  #Switch flag back
  #living[BEE == 1, noB := 0]
  
  # Generalized inverse does the actual solving
  #Invert A and multiple by b to get x (unknowns)
  x <- MASS::ginv(A, tol = .Machine$double.eps) %*% living[, b]
  
  #Assign unknown values
  living[, EEa := x * noEE]
  living[is.na(EE), EE := EEa]
  
  living[, B := x * noB]
  living[is.na(Biomass), Biomass := B]
  
  # detritus EE calcs
  living[, M0 := PB * (1 - EE)]
  living[, QBloss := QB]
  living[is.na(QBloss), QBloss := 0]
  #KYA fix Aug 2025
  #loss <- c((living[, M0] * living[, Biomass]) + 
  #            (living[, Biomass] * living[, QBloss] * living[, Unassim]),
  #          model[Type ==2, DetInput], 
  #          geardisc)
  #detinputs1  <- colSums(loss * detfate)
  #detinputs1 is "first pass" at det inputs, final detinputs is after initial EE
  loss <- c((living[, M0] * living[, Biomass]) + 
              (living[, Biomass] * living[, QBloss] * living[, Unassim]),
            rep(0,ndead), 
            geardisc) 
  detinputs1  <- colSums(loss * detfate + model[, DetInput])
  ## end fix
  detdiet    <- diet[(nliving + 1):(nliving + ndead), ]
  BQB        <- living[, Biomass * QB]
  detcons    <- as.matrix(detdiet) * BQB[col(as.matrix(detdiet))]
  detoutputs <- rowSums(detcons, na.rm = T)
  det_unused <- ifelse(detinputs1>detoutputs, detinputs1-detoutputs, 0.0)
  detinputs  <- detinputs1 + colSums(det_unused*detdetfate)
  EE         <- c(living[, EE], as.vector(detoutputs / detinputs))
  
  # added by kya
  # if a detritus biomass is put into the spreadsheet, use that and 
  # calculate PB.  If no biomass, but a PB, use that pb with inflow to 
  # calculate biomass.  If neither, use default PB=0.5, Bio = inflow/PB  
  # This is done because Ecosim requires a detrital biomass.
  
  Default_Detrital_PB <- 0.5 
  inDetPB <- model[(nliving + 1):(nliving + ndead), PB] 
  inDetB  <- model[(nliving + 1):(nliving + ndead), Biomass]
  DetPB   <- ifelse(is.na(inDetPB), Default_Detrital_PB, inDetPB)
  DetB    <- ifelse(is.na(inDetB), detinputs / DetPB, inDetB)
  DetPB   <- detinputs / DetB
  
  # Trophic Level calcs
  b             <- rep(1, ngroups)
  TLcoeff       <- matrix(0, ngroups, ngroups)
  diag(TLcoeff) <- rep(1, ngroups)
  gearcons      <- as.matrix(totcatchmat) / geartot[col(as.matrix(totcatchmat))]
  dimnames(gearcons) <- list(NULL, NULL)
  gearcons[is.na(gearcons)] <- 0
  dietplus <- as.matrix(diet)
  dimnames(dietplus) <- list(NULL, NULL)
  
  #Adjust for mixotrophs (partial primary producers) - #Moved this code up so that
  #it also impacted the EE calculation
  # mixotrophs <- which(model[, Type] > 0 & model[, Type] < 1)
  # mix.Q <- 1 - model[mixotrophs, Type]
  # for(i in seq_along(mixotrophs)){
  #   dietplus[, mixotrophs[i]] <- dietplus[, mixotrophs[i]] * mix.Q[i]
  # }
  #Adjust for diet import (Consumption outside model)
  import <- which(dietplus[nrow(diet), ] > 0)
  for(i in seq_along(import)){
    import.denom <- 1 - dietplus[nrow(diet), import[i]]
    dietplus[, import[i]] <- dietplus[, import[i]] / import.denom
  }
  dietplus <- dietplus[1:(nliving + ndead), ]
  dietplus <- rbind(dietplus, matrix(0, ngear, nliving))
  dietplus <- cbind(dietplus, matrix(0, ngroups, ndead), gearcons)
  TLcoeffA <- TLcoeff - dietplus
  TL       <- solve(t(TLcoeffA), b)     
  
  #kya changed these following four lines for detritus, and removing NAs
  #to match header file format (replacing NAs with 0.0s)
  Bplus  <- c(living[, Biomass], DetB, rep(0.0, ngear))
  
  PBplus <- model[, PB] 
  PBplus[(nliving + 1):(nliving + ndead)] <- DetPB
  PBplus[is.na(PBplus)] <- 0.0
  
  EEplus <- c(EE, rep(0.0, ngear))
  
  QBplus <- model[, QB]
  QBplus[is.na(QBplus)] <- 0.0
  
  GE[is.na(GE)] <- 0.0
  
  RemPlus <- model[, totcatch]
  RemPlus[is.na(RemPlus)] <- 0.0
  
  balanced <- list(Group    = model[, Group], 
                   TL       = TL, 
                   Biomass  = Bplus, 
                   PB       = PBplus, 
                   QB       = QBplus, 
                   EE       = EEplus, 
                   GE       = GE, 
                   Removals = RemPlus)
  
  M0plus  <- c(living[, M0], as.vector(detoutputs / detinputs))
  gearF   <- as.matrix(totcatchmat) / living[, Biomass][row(as.matrix(totcatchmat))]
  newcons <- as.matrix(nodetrdiet)  * BQB[col(as.matrix(nodetrdiet))]
  predM   <- as.matrix(newcons) / living[, Biomass][row(as.matrix(newcons))]
  predM   <- rbind(predM, detcons)
  morts   <- list(Group = model[Type < 3, Group], 
                  PB    = model[Type < 3, PB], 
                  M0    = M0plus, 
                  F     = gearF[1:(nliving + ndead), ], 
                  M2    = predM)
  
  # convert from levels to characters
  gnames <- as.character(balanced$Group)
  
  # cleanup before sending to sim -- C code wants 0 as missing value, not NA
  balanced$Biomass[is.na(balanced$Biomass)] <- 0
  balanced$PB[is.na(balanced$PB)]     <- 0
  balanced$QB[is.na(balanced$QB)]     <- 0
  balanced$EE[is.na(balanced$EE)]     <- 0
  balanced$GE[is.na(balanced$GE)]     <- 0
  model$BioAcc[is.na(model$BioAcc)]   <- 0
  model$Unassim[is.na(model$Unassim)] <- 0
  dietm                               <- as.matrix(diet)
  dimnames(dietm)                     <- list(c(gnames[1:(nliving+ndead)],"Import"), gnames[1:nliving])
  dietm[is.na(dietm)]                 <- 0
  landmatm                            <- as.matrix(landmat)
  dimnames(landmatm)                  <- list(gnames, gnames[(ngroups-ngear+1):ngroups])
  landmatm[is.na(landmatm)]           <- 0
  discardmatm                         <- as.matrix(discardmat)
  dimnames(discardmatm)               <- list(gnames, gnames[(ngroups-ngear+1):ngroups])
  discardmatm[is.na(discardmatm)]     <- 0
  detfatem                            <- as.matrix(detfate)
  dimnames(detfatem)                  <- list(gnames, gnames[(nliving+1):(nliving+ndead)])
  detfatem[is.na(detfatem)]           <- 0
  
  # Add names for output list
  out.Group   <- gnames;           names(out.Group) <- gnames
  out.type    <- model[, Type];    names(out.type) <- gnames
  out.TL      <- TL;               names(out.TL) <- gnames
  out.Biomass <- balanced$Biomass; names(out.Biomass) <- gnames
  out.PB      <- balanced$PB;      names(out.PB) <- gnames
  out.QB      <- balanced$QB;      names(out.QB) <- gnames
  out.EE      <- balanced$EE;      names(out.EE) <- gnames
  out.BA      <- model[, BioAcc];  names(out.BA) <- gnames
  out.Unassim <- model[, Unassim]; names(out.Unassim) <- gnames
  out.GE      <- balanced$GE;      names(out.GE) <- gnames    
  
  # list structure for sim inputs
  path.model <- list(NUM_GROUPS = ngroups,
                     NUM_LIVING = nliving,
                     NUM_DEAD   = ndead,
                     NUM_GEARS  = ngear,
                     Group      = out.Group,
                     type       = out.type,
                     TL         = out.TL,
                     Biomass    = out.Biomass,
                     PB         = out.PB,
                     QB         = out.QB,
                     EE         = out.EE,
                     BA         = out.BA,
                     Unassim    = out.Unassim,
                     GE         = out.GE,
                     DC         = dietm,
                     DetFate    = detfatem,
                     Landings   = landmatm,
                     Discards   = discardmatm)      
  
  #Define class of output
  class(path.model) <- 'Rpath'
  attr(path.model, 'eco.name') <- eco.name
  attr(path.model, 'eco.area') <- eco.area
  
  return(path.model)
}





  
  
