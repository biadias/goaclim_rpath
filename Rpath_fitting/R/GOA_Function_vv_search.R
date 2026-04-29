#------------------------------------------------------------------------------#
#AUTHORS: Bia Dias, Andy Whitehouse
#AFFILIATIONS: Alaska Fisheries Science Center / CICOES University of Washington
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#
# GOA Vulnerability Search Function
#------------------------------------------------------------------------------#




run_all_vv_searches <- function(scenarios, bal, hind_years, run.base.nll) {
  all_results <- list()
  
  for (scen_name in names(scenarios)) {
    message(paste("\n--- Starting VV Sensitivity Search for:", scen_name, "---"))
    
    current_scene <- scenarios[[scen_name]]
    
    # 1. Run baseline to get base likelihood table
    run_base <- rsim.run(current_scene, method = "AB", years = hind_years)
    base.like <- rsim.fit.table(current_scene, run_base)
    
    # 2. Run the sensitivity loop (using your sourced function)
    ptm <- proc.time()
    vv_loop_out <- vv_fit_loop(bal, current_scene, hind_years, run.base.nll, base.like)
    time_taken <- proc.time() - ptm
    
    # 3. Process results: Exponentiate VVs for interpretation
    vv_table <- vv_loop_out[[1]]
    vv_table$Vprey <- 1 + exp(vv_table$preyvul)
    vv_table$Vpred <- 1 + exp(vv_table$predvul)
    
    # 4. Save individual CSVs
    write.csv(vv_table, paste0("GOA_vv_search_", scen_name, ".csv"))
    
    # Store in master list
    all_results[[scen_name]] <- list(
      summary = vv_table,
      dnll_matrix = vv_loop_out[[2]],
      time = time_taken
    )
    
    message(paste("Completed", scen_name, "in", round(time_taken[3]/60, 2), "minutes."))
  }
  
  return(all_results)
}