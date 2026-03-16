#library(devtools)
#install_github('NOAA-EDAB/Rpath', ref="dev") 

library(Rpath)
library(dplyr)
source("R/xml_convert.r")

xml_unbal <- function(eiifile){
  unbal <- import.eiixml(eiifile)
  # Remove Biomass and add PB to detritus to estimate flows
  det_names <- unbal$model[Type==2]$Group
  unbal$model[Group%in%det_names, "Biomass"] <- NA
  unbal$model[Group%in%det_names, "PB"]      <- 1.0 # 1.0 means once a year
  unbal <- rpath.stanzas(unbal)
  check.rpath.params(unbal)
  return(unbal)
}

# File names for models
  WGOA_EwE_file <- "GOA/WGOA_9May2025.eiixml"
  #EGOA_EwE_file <- "rpath_files/EGOA_20250317_simpleDet.eiixml"
  #NBS_EwE_file  <- "rpath_files/aclim_nbs.eiixml"
  # EBS  
  #  Sbase <- "rpath_files/ebs_aclim3_76bio_base.csv"  # Base biomass, production, fishing, etc.
  #  Sdiet <- "rpath_files/ebs_aclim3_76bio_diet.csv"  # Diet matrix
  #  Sped  <- "rpath_files/ebs_aclim3_76bio_pedigree.csv"  # Data pedigree = quality of input data
  #  Sstz  <- "rpath_files/ebs_aclim3_76bio_stanzas.csv"  # Stanzas
  #  Sstg  <- "rpath_files/ebs_aclim3_76bio_stanza_groups.csv" # Stanza groups
    
  w.unbal <- xml_unbal(WGOA_EwE_file)
  
  w.bal   <- rpath(w.unbal) # balanced
  
  w.scene0 <- rsim.scenario(w.bal, w.unbal, years=1990:2089)
  
  w.run0   <- rsim.run(w.scene0, method="AB", years = 1990:2089)
  
  
  
  
  
  
  png("results/w_test.png",width=1067,height=600, type="cairo")  
  rsim.plot(w.run0)
  title("WGOA", cex.main=3, line=-1, adj=0.85)
  dev.off()  
  
  e.unbal <- xml_unbal(EGOA_EwE_file)
  n.unbal <- xml_unbal(NBS_EwE_file)
  s.unbal <- rpath.stanzas(read.rpath.params(Sbase, Sdiet, Sped, Sstg, Sstz)) # unbalanced
  check.rpath.params(s.unbal)
  
  n.bal   <- rpath(n.unbal) # balanced
  e.bal   <- rpath(e.unbal) # balanced
  s.bal   <- rpath(s.unbal) # balanced  
  w.bal   <- rpath(w.unbal) # balanced

  w.scene0 <- rsim.scenario(w.bal, w.unbal, years=1990:2089)
  e.scene0 <- rsim.scenario(e.bal, e.unbal, years=1990:2089)
  s.scene0 <- rsim.scenario(s.bal, s.unbal, years=1990:2089)  
  n.scene0 <- rsim.scenario(n.bal, n.unbal, years=1990:2089)
    
  w.run0   <- rsim.run(w.scene0, method="AB", years = 1990:2089)
  e.run0   <- rsim.run(e.scene0, method="AB", years = 1990:2089)  
  s.run0   <- rsim.run(s.scene0, method="AB", years = 1990:2089)
  n.run0   <- rsim.run(n.scene0, method="AB", years = 1990:2089)  
  
# FLATLINES
  png("results/w_test.png",width=1067,height=600, type="cairo")  
  rsim.plot(w.run0)
  title("WGOA", cex.main=3, line=-1, adj=0.85)
  dev.off()
  png("results/e_test.png",width=1067,height=600, type="cairo")
  rsim.plot(e.run0)
  title("EGOA", cex.main=3, line=-1, adj=0.85)
  dev.off()
  png("results/s_test.png",width=1067,height=600, type="cairo")
  rsim.plot(s.run0)
  title("EBS", cex.main=3, line=-1, adj=0.85)
  dev.off()
  png("results/n_test.png",width=1067,height=600, type="cairo")
  rsim.plot(n.run0)
  title("NBS", cex.main=3, line=-1, adj=0.85)
  dev.off()

# POLLOCK BUMP    
  w.scene1 <- adjust.forcing(w.scene0, "ForcedRecs", "walleye_pollock_adult", sim.year=1992, value=10)
  w.run1   <- rsim.run(w.scene1, method="AB", years = 1990:2039)
  png("results/w_pol.png",width=1067, height=600, type="cairo")
  rsim.plot(w.run1)
  title("WGOA", cex.main=3, line=-1, adj=0.85)
  dev.off()
  
  e.scene1 <- adjust.forcing(e.scene0, "ForcedRecs", "walleye_pollock_adult", sim.year=1992, value=10)
  e.run1   <- rsim.run(e.scene1, method="AB", years = 1990:2039)
  png("results/e_pol.png",width=1067, height=600, type="cairo")
  rsim.plot(e.run1)
  title("EGOA", cex.main=3, line=-1, adj=0.85)
  dev.off()

  s.scene1 <- adjust.forcing(s.scene0, "ForcedRecs", "pollock_adu", sim.year=1992, value=10)
  s.run1   <- rsim.run(s.scene1, method="AB", years = 1990:2039)
  png("results/s_pol.png",width=1067, height=600, type="cairo")
  rsim.plot(s.run1)
  title("EBS", cex.main=3, line=-1, adj=0.85)
  dev.off()
  
  n.scene1 <- adjust.forcing(n.scene0, "ForcedRecs", "pollock_adu", sim.year=1992, value=10)
  n.run1   <- rsim.run(n.scene1, method="AB", years = 1990:2039)
  png("results/n_pol.png",width=1067, height=600, type="cairo")
  rsim.plot(n.run1)
  title("NBS", cex.main=3, line=-1, adj=0.85)
  dev.off()
  
  source("code/rsim.plot.interactive.R")
  rsim.plot.interactive(w.run1)
  rsim.plot.interactive(e.run1)  
  rsim.plot.interactive(s.run1) 
  rsim.plot.interactive(n.run1) 

  
fup<-function(){source("code/FourSystems_functions.R")}
fup()  

steps <- 3650

w.markov <- flowmat(w.bal)
w.markov$state[,"small_phytoplankton"] <- 1
w.out <- matrix(0, nrow=steps, ncol=length(w.markov$state)); colnames(w.out) <- colnames(w.markov$state) 
for (i in 1:steps){
   w.out[i,] <- w.markov$state
   w.markov$state <- w.markov$state %*% w.markov$transition
}
  
e.markov <- flowmat(e.bal)
e.markov$state[,"small_phytoplankton"] <- 1
e.out <- matrix(0, nrow=steps, ncol=length(e.markov$state)); colnames(e.out) <- colnames(e.markov$state) 
for (i in 1:steps){
  e.out[i,] <- e.markov$state
  e.markov$state <- e.markov$state %*% e.markov$transition
}

s.markov <- flowmat(s.bal)
s.markov$state[,"sm_phytoplankton"] <- 1
s.out <- matrix(0, nrow=steps, ncol=length(s.markov$state)); colnames(s.out) <- colnames(s.markov$state) 
for (i in 1:steps){
  s.out[i,] <- s.markov$state
  s.markov$state <- s.markov$state %*% s.markov$transition
}

n.markov <- flowmat(n.bal)
n.markov$state[,"sm_phytoplankton"] <- 1
n.out <- matrix(0, nrow=steps, ncol=length(n.markov$state)); colnames(n.out) <- colnames(n.markov$state) 
for (i in 1:steps){
  n.out[i,] <- n.markov$state
  n.markov$state <- n.markov$state %*% n.markov$transition
}

png("results/dye_eup.png",width=1067, height=600, type="cairo")
plot(w.out[1:365,"euphausiids"],type="l",lwd=3,col="purple",
     xlab="months",ylab="Proportion dye in Euphausiids",cex.axis=1.5,cex.lab=1.5)
lines(e.out[1:365,"euphausiids"],lwd=3,col="red")
lines(s.out[1:365,"euphausiids"],lwd=3,col="green")
lines(n.out[1:365,"euphausiids"],lwd=3,col="blue")
legend("topright",legend=c("WGOA","EGOA","EBS","NBS"),lwd=3,col=c("purple","red","green","blue"),lty=1,cex=3)
dev.off()

png("results/dye_pol.png",width=1067, height=600, type="cairo")
plot(w.out[1:3650,"walleye_pollock_adult"],type="l",lwd=3,col="purple",ylim=c(0,0.06),
     xlab="months",ylab="Proportion dye in Adult Pollock",cex.axis=1.5,cex.lab=1.5)
lines(e.out[1:3650,"walleye_pollock_adult"],lwd=3,col="red")
lines(s.out[1:3650,"pollock_adu"],lwd=3,col="green")
lines(n.out[1:3650,"pollock_adu"],lwd=3,col="blue")
legend("topright",legend=c("WGOA","EGOA","EBS","NBS"),lwd=3,col=c("purple","red","green","blue"),lty=1,cex=3)
dev.off()

plot(w.out[1:3650,"pacific_cod_adult"],type="l",col="purple",ylim=c(0,.04))
lines(e.out[1:3650,"pacific_cod_adult"],col="red")
lines(s.out[1:3650,"pcod_adu"],col="green")
lines(n.out[1:3650,"pcod_adu"],col="blue")
legend("topright",legend=c("WGOA","EGOA","EBS","NBS"),col=c("purple","red","green","blue"),lty=1)


plot(w.out[1:365,"large_copepods"]+w.out[1:365,"small_copepods"],type="l",col="purple")
lines(e.out[1:365,"large_copepods"]+e.out[1:365,"small_copepods"],col="red")
lines(s.out[1:365,"copepods"],col="green")
lines(n.out[1:365,"copepods"],col="blue")
legend("topright",legend=c("WGOA","EGOA","EBS","NBS"),col=c("purple","red","green","blue"),lty=1)


w.out[,"euphausiids"]


write.csv(w.out,"wout.csv",row.names=F) 
write.csv(e.out,"eout.csv",row.names=F) 
write.csv(s.out,"sout.csv",row.names=F) 
write.csv(n.out,"nout.csv",row.names=F) 

# Checking group names #########################################################
# Groups in the west not the east
  w.unbal$model$Group[which(!(w.unbal$model$Group %in% e.unbal$model$Group))]

# Groups in the east not the west 
  e.unbal$model$Group[which(!(e.unbal$model$Group %in% w.unbal$model$Group))]
  
  
# Making a comparative table for basic balance stats ###########################
  e.groups <- data.frame(e.bal$Group,e.bal$TL,e.bal$Biomass,e.bal$PB,e.bal$QB,e.bal$EE, e.bal$GE)
  w.groups <- data.frame(w.bal$Group,w.bal$TL,w.bal$Biomass,w.bal$PB,w.bal$QB,w.bal$EE, w.bal$GE) 
  
  ew.groups <- e.groups %>%
               full_join(w.groups, by=join_by(e.bal.Group==w.bal.Group))  
  
  #plot(ew.groups$e.bal.PB/ew.groups$w.bal.PB)
  #text(ew.groups$e.bal.PB,ew.groups$w.bal.PB,ew.groups$e.bal.Group)
  
  write.csv(ew.groups,"ewgroups_balance_comp_14Mar25.csv",row.names=F)
  
# Some plots ###################################################################
  
  dat <- ew.groups %>% filter(!( is.na(e.bal.Biomass) | is.na(w.bal.Biomass) ))
  plot(dat$e.bal.TL,dat$w.bal.TL)
  abline(0,1)
  points(dat$e.bal.TL[dat$e.bal.Group=="euphausiids"], dat$w.bal.TL[dat$e.bal.Group=="euphausiids"],col="red", pch=15)
  
  plot(dat$e.bal.TL,dat$w.bal.TL, type='n')
  abline(0,1)
  text(dat$e.bal.TL,dat$w.bal.TL, dat$e.bal.Group, cex=0.5)
  
  
# Testing rsim #################################################################
w.scene0 <- rsim.scenario(w.bal, w.unbal, years=1990:2089)
e.scene0 <- rsim.scenario(e.bal, e.unbal, years=1990:2089)
  
w.run0   <- rsim.run(w.scene0, method="RK4", years = 1990:2089)
rsim.plot(w.run0)

e.run0   <- rsim.run(e.scene0, method="RK4", years = 1990:2089)
rsim.plot(e.run0)

w.scene1 <- adjust.fishing(w.scene0, "ForcedFRate", "sablefish_adult", sim.year=1995:2000, value=0.2)
w.run1   <- rsim.run(w.scene1, method="RK4", years = 1990:2089)
rsim.plot(w.run1)

e.scene1 <- adjust.fishing(e.scene0, "ForcedFRate", "sablefish_adult", sim.year=1995:2000, value=0.2)
e.run1   <- rsim.run(e.scene1, method="RK4", years = 1990:2089)
rsim.plot(e.run1)

# Interactive plots. Use the 

source("code/rsim.plot.interactive.R")
rsim.plot.interactive(n.run0, spname = "all", indplot = FALSE, palette="b_palette")


test<-rsim.plot.interactive(w.run1, spname = "all", indplot = FALSE, palette=colorspace::rainbow_hcl)
b_palette <- colorspace::rainbow_hcl 
b2_palette <- colorRamps::primary.colors
rsim.plot.interactive(w.run1, spname = "all", indplot = FALSE, palette = "b2_palette" )

test


  