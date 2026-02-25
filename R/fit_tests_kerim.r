# These tests assume you have run a fitting setup script (EBS, GOA or other), 
# and that you have unbal, bal, and scene_fit loaded.

library(foreach)
library(doParallel)

testfun <- function(vul_level, scene, years){
  library(Rpath)
  result=rep(NA, scene$params$NumPredPreyLinks+1)
  for (L in 1:(scene$params$NumPredPreyLinks+1)){
    #cat(L,vul_level,"\n"); flush.console()
    v_orig <- scene$params$VV[L]
    v_log  <- log(v_orig-1.0)
    scene$params$VV[L] <- 1.0 + exp(v_log + vul_level)
    result[L] <- rsim.fit.run(NA,NA,NA,
                              scene=scene, run_method="AB", run_years=years, verbose=F)
    scene$params$VV[L] <- v_orig
  }
  return(result)
}


# How many cores does your CPU have
n_cores <- detectCores()
n_cores

# 8 cores, so doing 7 levels
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)

scene_test <- scene_fit

v_levels <- c(-3, -2, -1, 0, 1, 2, 3)
results <- list()
results <- foreach(v=1:length(v_levels)) %dopar% {
  results[v] <- testfun(v_levels[v],scene_test,hind_years)
}

stopCluster(cl = cluster)

write.csv(
data.frame(from=scene_test$params$spname[scene_test$params$PreyFrom+1],
           to=scene_test$params$spname[scene_test$params$PreyTo+1],
           V1=results[[1]], V2=results[[2]], V3=results[[3]],
           V4=results[[4]], V5=results[[5]],V6=results[[6]],
           V7=results[[7]]),
"testout.csv",
row.names=F)

for (L in 1:scene_test$params$NumPredPreyLinks){
  cat(L,"\n"); flush.console()
  v_orig <- scene_test$params$VV[L]
  v_log <- log(v_orig-1.0)
  for (V in 1:length(v_levels)){
    scene_test$params$VV[L] <- 1.0 + exp(v_log + v_levels[V])
    result[L,V] <- rsim.fit.run(NA,NA,NA,
      scene=scene_test, run_method="AB", run_years=hind_years, verbose=F)
    scene_test$params$VV[L] <- v_orig
  }
}

data.frame(to=scene_test$params$spname[scene_test$params$PreyTo+1],
from=scene_test$params$spname[scene_test$params$PreyFrom+1],result)