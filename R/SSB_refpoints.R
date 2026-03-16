# NEED TO INSTALL OR PULL-COMPILE BRANCH
# install_github('NOAA-EDAB/Rpath', ref="ssb_output")

source("R/EBS_fitting_setup.R")
run1 <-rsim.run(scene_fit,"AB",hind_years)

# Using the base model only, no forcing, fitting etc
  run.years <- 1990:2089
  scene <- rsim.scenario(bal, unbal, years = run.years)
# zero out all gear effort
  scene$fishing$ForcedEffort[] <- 0 
# Set forced Frate to Ecopath F
  splist <- c(rpath.living(bal),rpath.detrital(bal))
  flist  <- (rowSums(bal$Landings) + rowSums(bal$Discards))/bal$Biomass
  for (sp in splist){
    scene$fishing$ForcedFRate[,sp] <- flist[sp]
  }
# Currently forced catch and forced F do not go to detritus, so need to keep
# that group constant.
  scene$forcing$ForcedBio[,"discards_offal"] <- bal$Biomass["discards_offal"]
  
# This only works for AB right now, does NOT work for RK4
  run1 <- rsim.run(scene, "AB", run.years)

  
this.species <- "pollock_adu" #"pcod_adu"

PB_ref <- bal$PB[this.species]
Frates <- seq(0, PB_ref*3, 0.01)
SSB <- Catch <- AdultBiomass <- AdultRecW <- AdultRecN <- rep(NA,length(Frates))
t <- 0
for (F in seq(0, PB_ref*3, 0.01)){
  cat(t,F,"\n"); flush.console()
  t <- t+1
  sceneA <- adjust.fishing(scene, "ForcedFRate", this.species, sim.year=run.years, value=F)  
  run.out <- rsim.run(sceneA, "AB", run.years)
  # tail function takes last n rows of matrix
  
    SSB[t]          <- tail(run.out$out_SSB, n=1)[,this.species]
    Catch[t]        <- tail(run.out$out_Catch,n=1)[,this.species] * 12
    AdultRecN[t]    <- tail(run.out$out_Nrec, n=1)[,this.species] 
    AdultRecW[t]    <- tail(run.out$out_Wrec, n=1)[,this.species] 
    AdultBiomass[t] <- tail(run.out$out_Biomass, n=1)[,this.species]
}

res <- data.frame(AdultBiomass,SSB,Catch,AdultRecN,AdultRecW)

par(mfrow=c(1,4))
plot(res$AdultBiomass,res$Catch,type="l")
lines(res$SSB,res$Catch,col="red")
plot(res$SSB,res$AdultRecN,type='l')
plot(res$SSB,res$AdultRecW,type='l')
plot(res$SSB,res$AdultRecN*res$AdultRecW,type='l')

# NEW MONTHLY OUTPUTS

# These output tables are 1 column per whole stanza species, with column
# labels set to the oldest stanza group for each species.
  run1$out_SSB  # spawning stock biomass
  run1$out_eggs # Eggs produced in SSB units - by default, equal to SSB
  run1$out_Ninf # N of the oldest age class (single oldest month)
  run1$out_Winf # W of the oldest age class (single oldest month)
  
# These output tables are 1 column per stanza (e.g. juv and adu separate
# columns)
  run1$out_Nrec # Number of recruits to youngest month for that stanza
  run1$out_Wrec # Weight of recruits to youngest month for that stanza
  
# Now a plot
  
# Make an age table for indexing SSB to resulting adult recruits N months
# later.  First make a better=labelled table of first adu recruitment month
# (columm #3 specific to EBS) and get the last month of the run.
  adu.rec.age <- scene_fit$stanzas$Age1[,3]
  names(adu.rec.age) <- scene_fit$stanzas$Oldest
  last.month <- dim(run1$out_SSB)[1]
    
this.species <- "pollock_adu"
# Line up SSB with the adult recruits the correct number of months later
# (I bet I get this off by 1 month no matter what I do here...)

Stock.Rec <- data.frame(
               SSB  = run1$out_SSB[1:(last.month-adu.rec.age[this.species]), this.species],
               Nrec = run1$out_Nrec[(1+adu.rec.age[this.species]):last.month, this.species],
               Wrec = run1$out_Wrec[(1+adu.rec.age[this.species]):last.month, this.species]
               )
Stock.Rec$Brec <- Stock.Rec$Nrec * Stock.Rec$Wrec 

par(mfrow=c(1,3))
plot(Stock.Rec$SSB, Stock.Rec$Nrec)
plot(Stock.Rec$SSB, Stock.Rec$Wrec)
plot(Stock.Rec$SSB, Stock.Rec$Brec)

  
  