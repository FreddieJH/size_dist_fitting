# About ========================================================================

# This script takes the cleaned data from data/processed and performs further 
# cleaning and summarising

# Packages =====================================================================

library(tidyverse)

# Mean sizes ===================================================================

# for all RLS and CBF species
for(pop in c("species", "ecoregion", "gridcell")){
  
  # using location for the CBF population-level if RLS is either at the ecoregion
  # level or the gridcell level
  cbf_pop <- ifelse(pop %in% c("ecoregion", "gridcell"), "location", pop)
  
  
  meansizes_rls <-
    get(paste0("obsdata_rls_", pop)) %>% 
    select(population, size_class, n) %>% 
    uncount(n) %>% 
    summarise(mean_size = mean(size_class),
              .by = population) %>%
    mutate(data = "rls")
  
  meansizes_cbf <- 
    get(paste0("obsdata_cbf_", cbf_pop)) %>% 
    uncount(n) %>% 
    summarise(mean_size = mean(size_class),
              .by = population) %>%  
    mutate(data = "cbf")
  
  
  
  bind_rows(meansizes_rls,meansizes_cbf) %>% 
    mutate(species = str_extract(population, ".*(?=__)")) %>% 
    assign(x = paste0("meansizes_", pop), 
           value = ., 
           envir = .GlobalEnv)
  
  rm(pop, cbf_pop, meansizes_rls, meansizes_cbf)
}


# Targeted species =============================================================

fished_species <- 
  tibble(species = c(obsdata_cbf_species$species, 
                     obsdata_rls_species$species)) %>% 
  mutate(genus = str_extract(species, ".*(?=\\s)")) %>% 
  filter((species %in% frdc$species)| 
           (genus %in% frdc$genus)|
           (species %in% targeted$species)) %>% 
  pull(species) %>% 
  unique()

rm(frdc, targeted)

# End ==========================================================================