library(tidyverse)

create.rpath.aggtable <- function(bal){return(data.frame(old=bal$Group, new=bal$Group, type=bal$type))}
rpath.aggtable.add    <- function(aggtable, old, new){
  if(length(unique(aggtable$type[aggtable$old %in% old])) !=1){
    warning("Group types in aggregation don't match - not added.")
    return(aggtable)
  }
  else{
     aggtable$new[aggtable$old %in% old] <- new
     return(aggtable)
  }   
}

# I'm coding it so the process starts with a BALANCED model, and returns
# an UNBALANCED one where all the unknowns are in EE and PC.  Note that this
# procedure may have numeric mistakes (unflagged) if the model is out-of-balance
# to start, particularly for detritus fates.

# Create a lookup table for old names and new names (new column is all old names to start)
  aggtab <- create.rpath.aggtable(w.bal) 
  
# Add some aggregations to the lookup table
  aggtab <- rpath.aggtable.add(aggtab, 
            old = c("benthic_detritus","pelagic_detritus"), 
            new = "detritus")
  
  aggtab <- rpath.aggtable.add(aggtab, 
            old = c("piscivorous_surface_birds","planktivorous_surface_birds",
                   "piscivorous_diving_birds","planktivorous_diving_birds"), 
            new = "birds")

# Look at the lookup table - pretty basic you can modify by hand too.
# The rpath.aggtable.add is "safer" to use because it flags/stops you
# from combining groups of different trophic types.
  view(aggtab)
  
####################################
#rpath.aggregate <- function(bal, aggtable)
bal <- w.bal
aggtable <- aggtab

#------
atab <- cbind(aggtable, data.frame(
              B=bal$Biomass, PB=bal$PB, QB=bal$QB, EE=bal$EE, BA=bal$BA, Unassim=bal$Unassim,
              DetFlow=bal$QB * bal$Biomass * bal$Unassim + (1.0-bal$EE) * bal$PB * bal$Biomass))            

# Unique keeps groups in order of first appearance
  n.groups <- unique(atab$new)

n.agg <- atab %>% 
  group_by(new) %>%
  summarize(type=mean(type),
            BB=sum(B),
            PP=sum(B*PB),
            QQ=sum(B*QB),
            PB=PP/BB,
            QB=QQ/BB,
            BA=sum(BA),
            UA=sum(B*QB*Unassim)/sum(B*QB),
            .groups="keep")

# o.agg is ordered agg (using n.groups to order values)
o.agg <- data.frame(n.groups) %>%
  left_join(n.agg,by=c("n.groups"="new"))

n.unbal <- create.rpath.params(o.agg$n.groups,o.agg$type)
 n.unbal$model$Biomass <- ifelse(o.agg$type<2, o.agg$BB, NA)
 n.unbal$model$PB      <- ifelse(o.agg$type<3, o.agg$PB, NA)
 n.unbal$model$QB      <- ifelse(o.agg$type<2, o.agg$QB, NA)
 n.unbal$model$EE      <- NA
 n.unbal$model$PC      <- NA
 n.unbal$model$BioAcc  <- ifelse(o.agg$type<3, o.agg$BA, NA)
 n.unbal$model$Unassim <- ifelse(o.agg$type<3, o.agg$UA, NA)
 n.unbal$model$DetInput<- ifelse(o.agg$type==2, 0, NA)

  
# DIETS
  QQ <- (bal$Biomass * bal$QB)[1:bal$NUM_LIVING]
  PP <- (bal$Biomass * bal$PB)[1:bal$NUM_LIVING]

  diet_flows <- t(QQ*t(bal$DC))

  diet_frame <- data.frame(                         
    from = rownames(diet_flows)[row(diet_flows)],
    to   = colnames(diet_flows)[col(diet_flows)],
    flow = as.vector(diet_flows)) %>%
      left_join(atab,by=c("from"="old")) %>%
      left_join(atab,by=c("to"="old")) %>%
      rename(n.from=new.x, n.to=new.y) %>%
      select(from,to,n.from,n.to,flow) %>%
      mutate(n.from=ifelse(from=="Import","Import",n.from)) %>%
      group_by(n.from,n.to) %>%
      summarize(flow=sum(flow),.groups="keep")

  

