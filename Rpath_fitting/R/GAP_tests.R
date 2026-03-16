
# Library is installed with:
#  devtools::install_github("afsc-gap-products/gapindex")
library(gapindex)

# Connect to Oracle - long with username and pw, needs to be done every time
  channel <- gapindex::get_connected()

# Fetch data - see function help for details
  gapindex_data <- gapindex::get_data(
    year_set = c(1982:2024),
    survey_set = "EBS", #One of c("GOA", "AI", "EBS", "NBS", "BSS")
    spp_codes = c(471,10285,21740),   
    haul_type = 3, abundance_haul = "Y", pull_lengths = T, channel = channel)  

# Calculate haul-level CPUE at species level, including filling in zeros
  cpue <- gapindex::calc_cpue(gapdata = gapindex_data)

# Calculate stratum-level biomass, population abundance, mean CPUE and 
# associated variances from the data pull and the cpue calcs
  biomass_stratum <- gapindex::calc_biomass_stratum(
    gapdata = gapindex_data, cpue = cpue)

# Now check it for our own purposes (not part of gap library)
library(tidyverse)

  ebs_biomass_results <- biomass_stratum %>%
    group_by(YEAR,SPECIES_CODE) %>%
    summarize(bio_tons_tot    = sum(BIOMASS_MT),
              bio_tons_varest = sum(BIOMASS_VAR),
              bio_tons_se     = sqrt(bio_tons_varest),
              bio_tons_cv     = bio_tons_se/bio_tons_tot,
              .groups="keep"
    )
  
  
## WGOA testing####################
  # Fetch data - see function help for details
  gapindex_data_goa <- gapindex::get_data(
    year_set = c(1982:2024),
    survey_set = "GOA", #One of c("GOA", "AI", "EBS", "NBS", "BSS")
    spp_codes = c(21740, 10110,21921,420),   
    haul_type = 3, abundance_haul = "Y", pull_lengths = T, channel = channel)    
  
# Calculate haul-level CPUE at species level, including filling in zeros
  cpue_goa <- gapindex::calc_cpue(gapdata = gapindex_data_goa)
  
  biomass_stratum_goa <- gapindex::calc_biomass_stratum(
    gapdata = gapindex_data_goa, cpue = cpue_goa) %>%
    mutate(BIOMASS_VAR=replace_na(BIOMASS_VAR,0))
  
  
 whole_goa_results <- biomass_stratum_goa %>%
   group_by(YEAR,SPECIES_CODE) %>%
   summarize(bio_tons_tot    = sum(BIOMASS_MT),
             bio_tons_varest = sum(BIOMASS_VAR),
             bio_tons_se     = sqrt(bio_tons_varest),
             bio_tons_cv     = bio_tons_se/bio_tons_tot,
             .groups="keep"
   )   
  
 wgoa_domains <-  c("20", "21", "22", "121", "122",    #  "Chirikof_shelf" 
                        "120", "220",                      #  "Chirikof_gully"
                        "221", "320", "420", "520",        #  "Chirikof_slope"
                        "130", "230", "232",               #  "Kodiak_shelf"
                        "30", "31","32", "33", "35","131", "132", "133", "134", #  "Kodiak_gully"
                        "231", "330", "430", "530",        #  "Kodiak_slope"
                        "10", "11", "12", "13","111",      #  "Shumagin_shelf"
                        "110", "112",                      #  "Shumagin_gully"
                        "210", "310", "410", "510"         #  "Shumagin_slope"
 ) 
   
 wgoa_results <- biomass_stratum_goa %>% 
   filter(STRATUM %in% wgoa_domains) %>%
   group_by(YEAR,SPECIES_CODE) %>%
   summarize(bio_tons_tot    = sum(BIOMASS_MT),
             bio_tons_varest = sum(BIOMASS_VAR),
             bio_tons_se     = sqrt(bio_tons_varest),
             bio_tons_cv     = bio_tons_se/bio_tons_tot,
             .groups="keep"
   )     
 
 
 
  
  