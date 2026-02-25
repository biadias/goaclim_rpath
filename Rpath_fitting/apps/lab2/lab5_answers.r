
## Answer key  1
   ebs.fscene <- adjust.fishing(ebs.fscene, "ForcedFRate", "P.cod", 
                   sim.year = 1992:2089, value = 0.1593405)
   bio1 <- pollock_bio
   cat1 <- pollock_catch 
   plot(bio1,cat1)



## Answer key  2
  f_levels <- seq(0, 2.0, 0.1)
  cf_levels <- seq(0,1,0.05)
  cod_landings <- pol_landings <- matrix(NA,length(f_levels),length(cf_levels))
  cod_bio <- pol_bio <- matrix(NA,length(f_levels),length(cf_levels))
  j<-0
  for (crate in cf_levels){
    j <- j+1
    ebs.fscene <- adjust.fishing(ebs.unfished.scene, "ForcedFRate", "P.cod", 
                    sim.year = 1992:2089, value = crate)
    i<-0
    for (frate in f_levels){
      i <- i+1
      ebs.fscene <- adjust.fishing(ebs.fscene, "ForcedFRate", "Pollock", 
                    sim.year = 1992:2089, value = frate)
      ebs.frun    <- rsim.run(ebs.fscene, method="AB", years=1990:2089)
      cod_landings[i,j] <- ebs.frun$annual_Catch["2089","P.cod"]
      pol_landings[i,j] <- ebs.frun$annual_Catch["2089","Pollock"]
      cod_bio[i,j]      <- ebs.frun$annual_Biomass["2089","P.cod"]
      pol_bio[i,j]      <- ebs.frun$annual_Biomass["2089","Pollock"]
      cat(crate,frate,"\n"); flush.console()
    }
  }
  

  fmat <- cfmat <- matrix(NA,length(f_levels),length(cf_levels))

  fmat[] <- f_levels[row(fmat)]
  cfmat[] <- f_levels[col(cfmat)]

  netvalue <- 1000*eco.area*(pol_landings*0.10 + cod_landings*0.50) -
              20000000*fmat - 80000000*cfmat
  image(netvalue)



  hist(1000*eco.area*(pol_landings*0.10 + cod_landings*0.50))
