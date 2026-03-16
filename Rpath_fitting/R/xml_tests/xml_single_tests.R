# Test scripts for bugfixes related to parsing of XML EwE models
# Done by KYA started 14-Aug-2025

library(Rpath)
library(dplyr)

# Specific path to KA's machine, for fitting test scripts
fup<-function(){
#source("R/xml_convert.r")
source("R/xml_test_functions.r")
}; fup()


#####################
# Testing a single model
#m <- "WGOA_9May2025"
m <- "Gulf of California"
eii_file <- paste("xml_examples/",m,".eiixml",sep="")
ewe_tables <- Rpath::import.eiixml(eii_file, verbose=T)
unbal <- Rpath::create.rpath.from.eiixml(eii_file, verbose=T)
Rpath.params <- rpath.stanzas(unbal)


bal <- rpath(Rpath.params)
unbal_test   <- rpath.ewe.comptable(m, "xml_examples")
test_results <- rpath.ewe.check.table(unbal_test)
for(r in 1:nrow(test_results)){
  if (test_results$flag[r]=="***"){
    warning(test_results$parameter[r], " of ", test_results$max_group[r],
            " diff: ", test_results$max_diff[r],
            " (TYPE ",unbal_test[test_results$max_group[r],"type"],")"
            )
  }
}

#Rpath.params$model$Biomass
bal <- rpath(Rpath.params)
scene <- rsim.scenario(bal,Rpath.params)
run0 <- rsim.run(scene,"AB",1:100)
X11(); rsim.plot(run0)
rsim.plot(run0,c("pelagic_detritus","benthic_detritus"))
dstep <- rsim.deriv(scene,sim.year=1)
det.frame <- data.frame(from=scene$params$spname[scene$params$DetFrom+1],
                        to=  scene$params$spname[scene$params$DetTo+1],
                        detfrac=scene$params$DetFrac)


mz_loss <- scene$params$MzeroMort * scene$params$B_BaseRef
ua_loss <- scene$params$FtimeQBOpt * scene$params$B_BaseRef * scene$params$UnassimRespFrac

det.frame$mz <- det.frame$detfrac * mz_loss[det.frame$from]
det.frame$ua <- det.frame$detfrac * ua_loss[det.frame$from]


###########################################
m <- "Campeche"
Rpath.params <- unbals[[m]]
eco.name = NA; eco.area = 1
#bal <- rpath(Rpath.params)
bal <- rpath(Rpath.params)
scene <- rsim.scenario(bal,Rpath.params)
run0 <- rsim.run(scene,"AB",1:5)
rsim.plot(run0,c("detritus","phytoplankton"))

####### Empty param set to test
utest <-  Rpath::create.rpath.params(
  group = c("Predators","Herbivores","Phytoplankton","Detritus","Fleet"),
  type =  c(0,0,1,2,3),
  stgroup = NA
)

