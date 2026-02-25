#------------------------------------------------------------------------------#
#AUTHORS: Bia Dias
#AFFILIATIONS: CICOES University of Washington/ Alaska Fisheries Science Center
#E-MAIL OF CORRESPONDENCE AUTHOR: bia.dias@noaa.gov
#
#Script to update the groups names of the bioenergetic_EBS_ACLIM2_bioen.csv file
#to match WGOA names. 
#------------------------------------------------------------------------------#

library(tidyverse)
library(janitor)
library(here)

bioen <- read.csv("data/bioenergetic_EBS_ACLIM2_bioen.csv")
survey_temps <- read.csv("WGOA_source_data/species_weighted_temp_WGOA.csv") %>% 
  filter(year==1990) %>% 
  select(race_group:total_catch_year_kg) %>% 
  rename(Species=race_group, mean_st1990=agg_weighted_surf_temp, mean_bt1990= agg_weighted_bot_temp)

bioen_v2 <- bioen %>%
  mutate(
    Species = case_when(
      Species == "arrowtooth_juv" ~ "arrowtooth_flounder_juvenile" ,
      Species == "arrowtooth_adu" ~ "arrowtooth_flounder_adult" ,
      Species == "atka"           ~ "atka_mackerel" ,
      Species == "fh_sole"        ~ "flathead_sole_adult" ,
      Species == "gr_turbot_juv"  ~ "MISC" ,
      Species == "gr_turbot_adu"  ~ "MISC" ,
      Species == "kamchatka"      ~ "MISC" ,
      Species == "north_rockfish" ~ "MISC" ,
      Species == "octopus"        ~ "octopus" ,
      Species == "oth_flatfish"   ~ "MISC" ,
      Species == "oth_rockfish"   ~ "MISC" ,
      Species == "pcod_juv"       ~ "pacific_cod_juvenile" ,
      Species == "pcod_adu"       ~ "pacific_cod_adult" ,
      Species == "ak_plaice"      ~ "MISC" ,
      Species == "pac_ocean_perch"~ "pacific_ocean_perch_adult" ,
      Species == "pollock_juv"    ~ "walleye_pollock_juvenile" ,
      Species == "pollock_adu"    ~ "walleye_pollock_adult" ,
      Species == "nr_sole"        ~ "MISC" ,
      Species == "rougheye_rock"  ~ "MISC" ,
      Species == "sablefish"      ~ "sablefish_adult" ,
      Species == "sharks"         ~ "pacific_sleeper_shark" ,
      Species == "shortraker_rock"~ "MISC" ,
      Species == "skates"         ~ "MISC" ,
      Species == "yf_sole"        ~ "MISC" ,
      Species == "herring"        ~ "pacific_herring_adult" ,
      Species == "halibut_juv"    ~ "pacific_halibut_juvenile" ,
      Species == "halibut_adu"    ~ "pacific_halibut_adult" ,
      Species == "lg_sculpins"    ~ "MISC" ,
      Species == "squids"         ~ "squid" ,
      Species == "capelin"        ~ "pacific_capelin" ,
      Species == "sandlance"      ~ "pacific_sandlance"
    )
  ) %>% 
  select(Species, Tmax, Topt, Q10) %>% filter(!Species=="MISC")
#  full_join(survey_temps, by= join_by(Species)) # not joining because we don't need the mean_bt1990 here 
  

write.csv(bioen_v2, "WGOA_source_data/WGOA_bioen.csv", row.names=FALSE)
                               