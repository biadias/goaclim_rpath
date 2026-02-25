# Test scripts for bugfixes related to parsing of XML EwE models
# Done by KYA started 14-Aug-2025

library(Rpath)
library(data.table)

options(warn = 1)

unbal.base <- copy(AB.params)

bal.base   <- rpath(unbal.base)
sp <- "cod"
BB.base <- as.numeric(bal.base$Biomass[sp])
PB.base <- as.numeric(bal.base$PB[sp])
QB.base <- as.numeric(bal.base$QB[sp])
EE.base <- as.numeric(bal.base$EE[sp])
GE.base <- as.numeric(bal.base$GE[sp])
BA.base <- as.numeric(bal.base$BA[sp])

outs <- rbind(NULL,c(2,2,2,2,2,as.numeric(c(bal.base$Biomass[sp], 
                         bal.base$PB[sp],
                         bal.base$QB[sp],
                         bal.base$EE[sp],
                         bal.base$GE[sp])), NA, NA, NA ))
colnames(outs)<-c("isB","isPB","isQB","isEE","isPC","B","PB","QB","EE","PC","Status","Warns","Err")

for (i in 0:31){
  unbal <- copy(unbal.base)
  unbal$model[Group==sp, Biomass:=  ifelse(bitwAnd(i, 1), BB.base, NA) ]
  unbal$model[Group==sp, PB:=       ifelse(bitwAnd(i, 2), PB.base, NA) ]
  unbal$model[Group==sp, QB:=       ifelse(bitwAnd(i, 4), QB.base, NA) ]
  unbal$model[Group==sp, EE:=       ifelse(bitwAnd(i, 8), EE.base, NA) ]
  unbal$model[Group==sp, ProdCons:= ifelse(bitwAnd(i,16), GE.base, NA) ]
  cat(i,"\n"); flush.console()
  outv <- as.numeric(bitwAnd(2^(0:4), i)>0)

  status <- NA
  status <- capture.output(check.rpath.params(unbal), type="output")
  warns <- NA
  warns <- capture.output(check.rpath.params(unbal), type="message")
  
  noerr <- TRUE; emsg <- NA
  tryCatch({
   bal.test <- rpath(unbal)
   res <- as.numeric(c(bal.test$Biomass[sp], 
                    bal.test$PB[sp],
                    bal.test$QB[sp],
                    bal.test$EE[sp],
                    bal.test$GE[sp]))
  }, error=function(e){
      cat(paste("\nERROR :",conditionMessage(e),"\n"))
      emsg  <<- conditionMessage(e)
      noerr <<- FALSE
      res   <<- rep(NA,5)
    }
  )
  
  outs <- rbind(outs,c(outv,res,status,paste(warns,collapse=" "),paste(emsg,collapse=" ")))
}

write.csv(outs,"xml_tests/all_param_combos2.csv",row.names=F)


