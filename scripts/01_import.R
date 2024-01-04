# About ========================================================================

# This script imports and cleans the raw data

# Two data sources: Reef Life Survey (RLS), and Cryptobenthic fishes (CBF)
# RLS has three population levels: species, gridcell (1x1 deg), ecoregion
# CBF has two population levels: species, location

# Packages =====================================================================

library(tidyverse)
library(arrow) # for parquet files

# Parameters ===================================================================

rls_min_bins  <- 4
rls_min_count <- 200
cbf_min_count <- 10
force_run <- FALSE

# Targeted species =============================================================

# Data from https://www.fish.gov.au/reports/species
frdc <- read_csv("input/data/data_cleaning/frdc_fished.csv", 
                 show_col_types = FALSE) 

# From Audzijonyte et al (2020) Nat. Evo. Eco.
targeted <- 
  read_csv("input/data/data_cleaning/fishing_intensity.csv", 
           show_col_types = FALSE) %>% 
  filter(fishing_intensity > 1) # only showing targeted species

# RLS import ===================================================================

# Body size bins used in RLS surveys (=middle of bin)
rls_bin_breaks <- 
  c(2.5, 5.0, 7.5,  10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 
    50.0, 62.5, 75.0, 87.5, 100.0, 112.5, 125.0, 137.5, 150.0, 
    162.5, 175.0, 187.5, 200.0, 250.0, 300.0, 350.0, 400.0)

rls_bin_table <-
  tibble(size_class = c(0, rls_bin_breaks, 500)) %>% 
  mutate(
    size_min = (size_class + lag(size_class))/2,
    size_max = lead(size_min)
  ) %>% 
  filter(size_class %in% c(rls_bin_breaks, 500))


for(pop in c("species", "ecoregion", "gridcell")){
  
  filename <-
    paste0("input/data/processed/rls_obsdata_", 
           pop, "_", 
           rls_min_bins,  "_",
           rls_min_count, ".csv")

  if(!file.exists(filename)|force_run){
    
    raw_data <- 
      read_parquet("input/data/raw/data_obs_cleaned.parquet") %>% 
      left_join(read_parquet("input/data/raw/survey_list_m1_aus.parquet"), 
                by = join_by(survey_id)) %>% 
      rename(species = species_name) %>% 
      mutate(ecoregion = str_replace_all(ecoregion, "/", "-"), 
             ecoregion = str_replace_all(ecoregion, "Great Barrier Reef", "GBR"), 
             ecoregion = str_replace_all(ecoregion, "Central", "C"), 
             ecoregion = str_replace_all(ecoregion, "Northern", "N"), 
             ecoregion = str_replace_all(ecoregion, "Southern", "S"), 
             lat_grid = round(latitude), 
             lon_grid = round(longitude), 
             gridcell = paste(lat_grid, lon_grid, sep = "_"), 
             population = paste(species, !!sym(pop), sep = "__")) 
    
    n_transects <- 
      raw_data %>% 
      select(population, survey_date) %>% 
      distinct() %>% 
      count(population, name = "n_transects")
    
    simple_data <- 
      raw_data %>% 
      count(population, size_class, wt = n) 
    
    filtered_data <- 
      simple_data %>% 
      add_count(population, name = "n_sizebins") %>% 
      add_count(population, wt = n, name = "population_n") %>%
      filter(n_sizebins >= rls_min_bins,
             population_n >= rls_min_count) %>%
      arrange(desc(population_n)) 
    
    clean_data <-
      filtered_data %>% 
      left_join(rls_bin_table, by = join_by(size_class)) %>% 
      left_join(n_transects, by = join_by(population)) %>% 
      mutate(p_obs = n/population_n) %>% 
      mutate(species = str_extract(population, ".*(?=__)")) %>% 
      mutate(genus = str_extract(species, ".*(?=\\s)")) %>% 
      mutate(fished = (species %in% frdc$species)| 
               (genus %in% frdc$genus)|
               (species %in% targeted$species)) %>% 
      select(-genus)
    
    assign(paste0("obsdata_rls_", pop), 
           value = clean_data, 
           envir = .GlobalEnv)
    
    write_csv(clean_data, 
              filename)
    
    rm(pop, raw_data, simple_data, filtered_data, clean_data, n_transects) 
    
  } else {
    
    assign(paste0("obsdata_rls_", pop), 
           value = read_csv(filename, show_col_types = FALSE), 
           envir = .GlobalEnv)
    
    print(paste("RLS obsdata at", pop, "level imported."))
  }
  
  rm(filename)
}

# CBF import ===================================================================

for(pop in c("species", "location")){
  
  filename <-
    paste0("input/data/processed/cbf_obsdata_", 
           pop, "_", 
           cbf_min_count, ".csv")
  
  
  if(!file.exists(filename)|force_run){
    
    raw_data <- 
      read_csv("input/data/raw/cbf_raw.csv", 
             show_col_types = FALSE) %>%
      mutate(sl_cm = SL/10,
             tl_cm = TL/10) %>% # now in cm, not mm
      select(location = Location,
             species = sciname,
             tl_cm
      ) 
    
    simple_data <- 
      raw_data %>%
      mutate(size_class = tl_cm %>% round(1), 
             size_min = size_class - 0.05, 
             size_max = size_class + 0.05) %>% 
      filter(str_detect(species, "^[A-Z]{1}[a-z]+\\s[a-z]+$")) %>% 
      mutate(population = paste(species, !!sym(pop), sep = "__")) %>% 
      count(population, size_class, size_min, size_max)

    filtered_data <- 
      simple_data %>% 
      add_count(population, wt = n,  name = "population_n") %>%
      filter(population_n >= cbf_min_count) %>% 
      arrange(desc(population_n)) 
    
    
    clean_data <-
      filtered_data %>% 
      mutate(p_obs = n/population_n) %>% 
      mutate(species = str_extract(population, ".*(?=__)"),
             genus = str_extract(species, ".*(?=\\s)")) %>% 
      mutate(fished = (species %in% frdc$species)| 
               (genus %in% frdc$genus)|
               (species %in% targeted$species)) %>% 
      select(-genus)
    
    
    
    assign(paste0("obsdata_cbf_", pop), 
           value = clean_data, 
           envir = .GlobalEnv)
    
    write_csv(clean_data, 
              filename)
    
    rm(raw_data, simple_data, filtered_data, clean_data)
    
  } else {
    
    assign(paste0("obsdata_cbf_", pop), 
           value = read_csv(filename, show_col_types = FALSE), 
           envir = .GlobalEnv)
    
    print(paste("CBF obsdata at", pop, "level imported."))
  }
  
  assign(
    paste0("cbf_npops_", pop),
    get(paste0("obsdata_cbf_", pop)) %>% 
      pull(population) %>%
      n_distinct()
  )
  
  rm(pop)
}

# End ==========================================================================