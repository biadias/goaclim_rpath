
#library(dplyr)
#library(Rpath)

source("../R/xml_convert.r")

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
  unbal <- create.rpath.object.from.xml(eiifile)
  if(unbal$stanzas$NStanzaGroups>0){
    unbal <- rpath.stanzas(unbal)
  }
  Rpath::check.rpath.params(unbal)
  bal <- Rpath::rpath(unbal) 
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
  rpath.dat <- data.frame(group=bal$Group,tl=bal$TL, biomass=bal$Biomass,
                          pb=bal$PB,qb=bal$QB, EE=bal$EE, pc=bal$GE, 
                          ba=bal$BA)[c(rpath.living(bal), rpath.detrital(bal)),]
  
  # Joining rpath and ewe
  rpath_ewe_table <- rpath.dat %>% 
    dplyr::left_join(csv.out, by="group") %>%
    dplyr::mutate(
      tl_test      = abs(tl-Trophic.level)/Trophic.level, 
      biomass_test = abs(biomass-Biomass.in.habitat.area)/Biomass.in.habitat.area, 
      pb_test      = abs(pb-ifelse(is.na(Production.biomass),Total.mortality,Production.biomass))/
                            ifelse(is.na(Production.biomass),Total.mortality,Production.biomass),       
      qb_test      = abs(qb-Consumption.biomass)/Consumption.biomass,
      ee_test      = abs(EE-Ecotrophic.Efficiency)/Ecotrophic.Efficiency,
      pc_test      = abs(pc-Production.consumption)/Production.consumption,
      ba_test      = abs(ba-Biomass.accumulation)/Biomass.accumulation
  )
  return(rpath_ewe_table)
}  

rpath.ewe.checks <- function(rp.ewe){
  cat("Max differences (proportion)\n")
  v=rp.ewe$tl_test; catcol(paste("TL",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
  v=rp.ewe$biomass_test; catcol(paste("B ",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
  v=rp.ewe$pb_test; catcol(paste("PB",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
  v=rp.ewe$qb_test; catcol(paste("QB",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
  v=rp.ewe$ee_test; catcol(paste("EE",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
  v=rp.ewe$pc_test; catcol(paste("PC",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
  #v=rp.ewe$ba_test; catcol(paste("BA",rp.ewe$group[which.max(v)],":",max(v,na.rm=T)),ifelse(max(v,na.rm=T)>=0.01,31,29))
}

# Single use test
  model_name <- "Tampa Bay"
  cat("\n",model_name,'\n')
  rpath_ewe <- rpath.ewe.comptable(model_name,"xml_examples")
  rpath.ewe.checks(rpath_ewe)
  
# Get a list of all files in the xml_examples directory that have
# a '-Basic estimates.csv' output present.  It's assumed there's
# a .eiixml file of the same name in the same directory.
ewe_csvs <- list.files("xml_examples", pattern="-Basic estimates.csv")
model_names <- stringr::str_remove_all(ewe_csvs,"-Basic estimates.csv")  

#broken_models <- c("Antarctic", "Barents Sea", "Denmark, Faroe Islands")
#models_filtered <- model_names[!(model_names %in% broken_models)]

models_filtered <- model_names

tests<-list()
for (m in models_filtered){
  tryCatch({
    cat("\n",m,'\n')
    tests[[m]] <- rpath.ewe.comptable(m, "xml_examples")
    rpath.ewe.checks(tests[[m]])
  },error=function(e){catcol(paste("ERROR :",conditionMessage(e)))})    
}  






  
  
