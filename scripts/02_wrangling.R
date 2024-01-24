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


# Plotting data ================================================================


for(pop in c("species", 
             "ecoregion", 
             "gridcell")){
  
  cbf_pop <- ifelse(pop %in% c("gridcell", "ecoregion"), "location", pop)
  
  norm_data_rls <- 
    paste0("output/models/summary/pars/rls_", pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "rls")
  
  lnorm_data_rls <- 
    paste0("output/models/summary/pars/rls_", pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "rls")
  
  norm_data_cbf <- 
    paste0("output/models/summary/pars/cbf_", cbf_pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "cbf")
  
  lnorm_data_cbf <- 
    paste0("output/models/summary/pars/cbf_", cbf_pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "cbf")
  
  norm_data_rls_cv <- 
    paste0("output/models/summary/cv/rls_", pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_norm")) %>% 
    mutate(dat = "rls")
  
  lnorm_data_rls_cv <- 
    paste0("output/models/summary/cv/rls_", pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_lnorm")) %>% 
    mutate(dat = "rls")
  
  norm_data_cbf_cv <- 
    paste0("output/models/summary/cv/cbf_", cbf_pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_norm")) %>% 
    mutate(dat = "cbf")
  
  lnorm_data_cbf_cv <- 
    paste0("output/models/summary/cv/cbf_", cbf_pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_lnorm")) %>% 
    mutate(dat = "cbf")
  
  bimodal_pops_rls <- 
    paste0("input/data/data_cleaning/bimodal_rls_", pop, ".csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    pull(population)
  
  bimodal_pops_cbf <- 
    paste0("input/data/data_cleaning/bimodal_cbf_", cbf_pop, ".csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    pull(population)
  
  
  bind_rows(norm_data_rls, 
            norm_data_cbf) %>% 
    rename(logl_50_norm = logl_50) %>% 
    left_join(bind_rows(lnorm_data_rls, 
                        lnorm_data_cbf) %>% 
                rename(logl_50_lnorm = logl_50), 
              by = join_by(population, dat)) %>% 
    left_join(bind_rows(norm_data_rls_cv, 
                        norm_data_cbf_cv),
              by = join_by(population, dat)) %>% 
    left_join(bind_rows(lnorm_data_rls_cv, 
                        lnorm_data_cbf_cv),
              by = join_by(population, dat)) %>% 
    mutate(species = str_extract(population, "^.*(?=__)")) %>% 
    mutate(fished = species %in% fished_species) %>% 
    mutate(bimodal = population %in% c(bimodal_pops_rls, 
                                       bimodal_pops_cbf)) %>% 
    mutate(normal_better = case_when(
      is.na(mu_50) ~ FALSE,
      is.na(meanlog_50) ~ TRUE,
      logl_50_norm > logl_50_lnorm ~ TRUE, 
      logl_50_lnorm > logl_50_norm ~ FALSE
    )) %>% 
    mutate(better_dist = ifelse(normal_better, "normal", "lognormal")) %>% 
    mutate(cov_pref = ifelse(normal_better, cv_50_norm, cv_50_lnorm)) %>% 
    assign(x = paste0("plotdata_", pop),  
           value = ., 
           envir = .GlobalEnv)
  
}



# End ==========================================================================