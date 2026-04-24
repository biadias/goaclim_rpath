source("code/bioenergetic_projections_v2.R")

wgoa_add_bioenergetics_v2 <- function(ssp){
  #hind_years <- 1991:2020 # just for knowledge, it won't be used. 
  prefix <- "tdc"
  cons <- TRUE
  resp <- TRUE
  #ssp <- "ssp126" # or ssp245, ssp585, etc.
  scene_bd <- rsim.scenario(bal, unbal, years=1991:2099)
  
  #–– Bioenergetics fitting ––#
  # KYA 23-Jun-25 below var needs to be globally available to other functions - exporting 
  bioen_sp_noceph <<- setdiff(bioen_sp, c("octopus", "squid"))
  
  if (cons) {
    n_hind_tdc <- nrow(tdc_hind)
    n_fs_tdc   <- nrow(scene_bd$forcing$ForcedSearch)
    n_use_tdc  <- min(n_hind_tdc, n_fs_tdc)
    
    use_sp_tdc <- intersect(bioen_sp_noceph, colnames(scene_bd$forcing$ForcedSearch))
    
    # --- 1. HINDCAST BULK ASSIGNMENT ---
    # tdc_hind already contains surface temps for pelagics and bottom temps for demersals
    scene_bd$forcing$ForcedSearch[1:n_use_tdc, use_sp_tdc] <- tdc_hind[1:n_use_tdc, use_sp_tdc]
    
    # --- 2. PROJECTION BULK ASSIGNMENT ---
    proj1   <- bioen_params(ssp)[[1]]
    end1    <- min(n_use_tdc + nrow(proj1), n_fs_tdc)
    scene_bd$forcing$ForcedSearch[(n_use_tdc + 1):end1, use_sp_tdc] <- proj1[1:(end1 - n_use_tdc), use_sp_tdc]
  }
  
  if (resp) {
    n_hind_tdr <- nrow(tdr_hind)
    n_fs_tdr   <- nrow(scene_bd$forcing$ForcedActresp)
    n_use_tdr  <- min(n_hind_tdr, n_fs_tdr)
    
    use_sp_tdr <- intersect(bioen_sp_noceph, colnames(scene_bd$forcing$ForcedActresp))
    
    # --- 1. HINDCAST BULK ASSIGNMENT ---
    scene_bd$forcing$ForcedActresp[1:n_use_tdr, use_sp_tdr] <- tdr_hind[1:n_use_tdr, use_sp_tdr]
    
    # --- 2. PROJECTION BULK ASSIGNMENT ---
    proj2   <- bioen_params(ssp)[[2]]
    end2    <- min(n_use_tdr + nrow(proj2), n_fs_tdr)
    scene_bd$forcing$ForcedActresp[(n_use_tdr + 1):end2, use_sp_tdr] <- proj2[1:(end2 - n_use_tdr), use_sp_tdr]
  }
  
  return(scene_bd)
}
