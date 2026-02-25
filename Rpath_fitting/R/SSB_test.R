# NEED TO INSTALL OR PULL-COMPILE BRANCH
# install_github('NOAA-EDAB/Rpath', ref="ssb_output")

source("R/EBS_fitting_setup.R")

# This only works for AB right now, does NOT work for RK4
run1 <- rsim.run(scene_fit, "AB", hind_years)

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

# These plots are pretty messy because they're plotting monthy recruitment
# during hindcast forcing
  plot(Stock.Rec$SSB, Stock.Rec$Nrec)
  plot(Stock.Rec$SSB, Stock.Rec$Wrec)
  plot(Stock.Rec$SSB, Stock.Rec$Brec)


# Now let's look at reference points!

# Using a climate-free base run here  
  scene_ssb <- rsim.scenario(bal, unbal, years = 1:100)

# Standard F setup - zero out all gear effort and set forced Frates to Ecopath F, fix discards
  scene_ssb$fishing$ForcedEffort[] <- 0 
  splist <- c(rpath.living(bal),rpath.detrital(bal))
  flist  <- (rowSums(bal$Landings) + rowSums(bal$Discards))/bal$Biomass
  for (sp in splist){scene_ssb$fishing$ForcedFRate[,sp] <- flist[sp]}
  scene_ssb$forcing$ForcedBio[,"discards_offal"] <- bal$Biomass["discards_offal"]
  
# Loop through various Frates
  this.species <- "pcod_adu"
  Frates <- seq(0,1,0.01)
  frate <- SSB <- adu_bio <- SSB_prop <- adu_bio_prop <- catch <- rep(NA,length(Frates))
  
  for (ff in Frates){
    cat(ff,"\n"); flush.console()
    ind <- which(Frates==ff)
    scene_ssb$fishing$ForcedFRate[,this.species] <- ff
    frun <- rsim.run(scene_ssb, "AB", 1:100) 
    frate[ind]        <- ff
    SSB[ind]          <- frun$out_SSB[1200, this.species]
    adu_bio[ind]      <- frun$out_Biomass[1200, this.species]
    catch[ind]        <- frun$annual_Catch[100, this.species]
    SSB_prop[ind]     <- SSB[ind]/SSB[1] # SSB[1] is unfished state
    adu_bio_prop[ind] <- adu_bio[ind]/adu_bio[1]
  }
  
  result <- data.frame(frate,SSB,adu_bio,catch,SSB_prop,adu_bio_prop)
    

  plot(result$SSB,result$catch)  
# IF ALL OF THE SSB is contained within the adult pool (i.e. younger stanzas
# don't contribute to SSB), then adult F will equal SSB F, because there's
# age-specific selectivity.  So F_MSY will be the same F for both.
  frate[which.max(catch)]
  
# However the same F doesn't result in the same biomass change, so the SSB_40
# and Adu_bio 40 will be different
  
  # Find Frate closest to target adu_bio_40
  frate[which.min(abs(adu_bio_prop - 0.40))]
    
  # Find Frate closest to targed SSB_40
  frate[which.min(abs(SSB_prop - 0.40))]  
    
  