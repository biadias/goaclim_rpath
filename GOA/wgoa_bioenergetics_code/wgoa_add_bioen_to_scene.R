source("GOA/wgoa_bioenergetics_code/bioenergetic_projections.r")

wgoa_add_bioenergetics <- function(ssp){
#hind_years <- 1991:2020 # just for knowledge, it won't be used. 
prefix <- "tdc"
cons <- TRUE
resp <- TRUE
#ssp <- "ssp126" # or ssp245, ssp585, etc.
scene_bd <- rsim.scenario(bal, unbal, years=1990:2089)
#run_bd  <- rsim.run(scene_bd, method="AB", years = 1990:2089)

#–– Bioenergetics fitting ––#
# KYA 23-Jun-25 below var needs to be globally available to other functions - exporting 
bioen_sp_noceph <<- setdiff(bioen_sp, c("octopus", "squids"))

if (cons) {
  ##–– Consumption multiplier ––##
  # how many hindcast rows do we have?
  n_hind_tdc <- nrow(tdc_hind_bt)
  # how many rows does the model matrix have?
  n_fs_tdc  <- nrow(scene_bd$forcing$ForcedSearch)
  # only fill as many as will fit
  n_use_tdc <- min(n_hind_tdc, n_fs_tdc)
  
  # only use species that actually exist in ForcedSearch
  use_sp_tdc <- intersect(bioen_sp_noceph,
                          colnames(scene_bd$forcing$ForcedSearch))
  if (length(use_sp_tdc)==0) stop("No matching cons columns in ForcedSearch")
  
  # hindcast
  scene_bd$forcing$ForcedSearch[1:n_use_tdc, use_sp_tdc] <-
    tdc_hind_bt[1:n_use_tdc, use_sp_tdc]
  
  # projection
  proj1   <- bioen_params(ssp)[[1]]
  n_proj1 <- nrow(proj1)
  # cap end row so we don’t overshoot
  end1    <- min(n_use_tdc + n_proj1, n_fs_tdc)
  
  scene_bd$forcing$ForcedSearch[
    (n_use_tdc + 1):end1,
    use_sp_tdc
  ] <- proj1[1:(end1 - n_use_tdc), use_sp_tdc]
}

if (resp) {
  ##–– Respiration multiplier ––##
  n_hind_tdr <- nrow(tdr_hind_bt)
  n_fs_tdr   <- nrow(scene_bd$forcing$ForcedActresp)
  n_use_tdr  <- min(n_hind_tdr, n_fs_tdr)
  
  use_sp_tdr <- intersect(bioen_sp_noceph,
                          colnames(scene_bd$forcing$ForcedActresp))
  if (length(use_sp_tdr)==0) stop("No matching resp columns in ForcedActresp")
  
  # hindcast
  scene_bd$forcing$ForcedActresp[1:n_use_tdr, use_sp_tdr] <-
    tdr_hind_bt[1:n_use_tdr, use_sp_tdr]
  
  # projection
  proj2   <- bioen_params(ssp)[[2]]
  n_proj2 <- nrow(proj2)
  end2    <- min(n_use_tdr + n_proj2, n_fs_tdr)
  
  scene_bd$forcing$ForcedActresp[
    (n_use_tdr + 1):end2,
    use_sp_tdr
  ] <- proj2[1:(end2 - n_use_tdr), use_sp_tdr]
}

 return(scene_bd)
}
#------------------------------------------------------------------------------#  
