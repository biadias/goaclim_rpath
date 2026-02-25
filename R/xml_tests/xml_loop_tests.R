# Test scripts for bugfixes related to parsing of XML EwE models
# Done by KYA started 14-Aug-2025

library(Rpath)
library(dplyr)

# Specific path to KA's machine, for fitting test scripts
fup<-function(){
source("R/xml_convert.r")
source("R/xml_test_functions.r")
}; fup()

ewe_csvs <- list.files("xml_examples", pattern="-Basic estimates.csv")
model_names <- stringr::str_remove_all(ewe_csvs,"-Basic estimates.csv")
models_filtered <- model_names
#models_filtered <- c("Albatross Bay","Tampa Bay")
tests<-unbals<-list()
options(warn=2)
noerr<-vector()
for (m in models_filtered){
  cat("\n##",m,"\n")
  tryCatch({
    tables <- import.eiixml(paste("xml_examples/",m,".eiixml",sep=""))
    unbals[[m]] <- create.rpath.object.from.xml(paste("xml_examples/",m,".eiixml",sep=""))
    tests[[m]]   <- rpath.ewe.comptable(m, "xml_examples")
    test_results <- rpath.ewe.check.table(tests[[m]])
    for(r in 1:nrow(test_results)){
      if (test_results$flag[r]=="***"){
        warning(test_results$parameter[r], " of ", test_results$max_group[r],
                " diff (proportion of EwE): ", test_results$max_diff[r])
      }
    }
    print(test_results) #print(knitr::kable(test_results,"simple")); cat("\n")
    noerr<-append(noerr,m)
  },error=function(e){cat(paste("\nERROR :",conditionMessage(e),"\n"))
  })    
}  

EE_list <- EE_yes <-NULL
for(mind in names(tests)){
  cat (mind,"\n"); flush.console()#cat(models_filtered[mind],"\n"); flush.console()
  tryCatch({  
  mod <- unbals[[mind]]
  if(mod$stanzas$NStanzaGroups == 0){
    Rpath.params <- mod
  } else {
    Rpath.params <- rpath.stanzas(mod)
  }
  #; eco.name = NA; eco.area = 1
  baltest <- rpath(Rpath.params)
  ind <-(tests[[mind]]$Ecotrophic.Efficiency==0)
  if(sum(ind)>0){
   EE_list <- rbind(EE_list, 
    data.frame(model=mind,
               group=tests[[mind]]$group[ind], 
               rpathEE=tests[[mind]]$EE[ind]))
  }
  EE_yes <- rbind(EE_yes,
                  data.frame(model=mind,
                             group=tests[[mind]]$group[!ind], 
                             rpathEE=tests[[mind]]$EE[!ind],
                             eweEE  =tests[[mind]]$Ecotrophic.Efficiency[!ind]))
  },error=function(e){warning(paste("\nERROR :",conditionMessage(e),"\n"))
  }) 
}

