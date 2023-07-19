# ==============================================================================
#  DATA WRANGLING
# ==============================================================================

# RLS data wrangling -----------------------------------------------------------

# Cleaning data into three spatial scales: species level, 
# gridcell level (1x1 degree cells), and ecoregion level.

# Species-level RLS data
RLS_data_obs_species <-
  RLS_data_obs %>% 
  rename(population = species_name) %>% 
  count(population, size_class) %>% 
  add_count(population, name = "n_sizebins") %>% 
  add_count(population, wt = n, name = "population_n") 


# sample of 20 species
RLS_data_obs_20spp <-
  RLS_data_obs_sample %>% 
  count(species_name, size_class)

# Ecoregion-level (Marine Ecoregions Of the World)
RLS_data_obs_ecoregion <-
  RLS_data_obs %>% 
  left_join(RLS_data_survs, by = join_by(survey_id)) %>% 
  mutate(population = paste(species_name, ecoregion, sep = "__")) %>% 
  count(population, ecoregion, species_name, size_class) %>% 
  add_count(ecoregion, species_name, name = "n_sizebins") %>% 
  add_count(ecoregion, species_name, wt = n, name = "population_n")

# Gridcell-level data (1x1 lat-lon degree gridcell)
RLS_data_obs_gridcell <-
  RLS_data_obs %>% 
  left_join(RLS_data_survs, by = join_by(survey_id)) %>% 
  mutate(population = paste(species_name, lat_lon, sep = "__")) %>% 
  count(population, species_name, size_class, wt = n) %>% 
  add_count(population, name = "n_sizebins") %>% 
  add_count(population, wt = n, name = "population_n") 

# No sites go over a lat-lon boundary
RLS_data_obs %>% 
  left_join(RLS_data_survs, by = join_by(survey_id)) %>% 
  select(site_code, 
         lat_lon) %>% 
  distinct() %>% 
  count(site_code) %>% 
  filter(n > 1)

# CBF data wrangling -----------------------------------------------------------

# Cleaning data into two spatial scales: species level, and location level.

# Species-level data
CBF_data_obs_species <- 
  CBF_data_obs %>% 
  mutate(population = paste(species, sep = "_")) %>% 
  select(population, tl) %>% 
  add_count(population, name = "population_n")

# Location-level data
CBF_data_obs_location <- 
  CBF_data_obs %>% 
  mutate(population = paste(location, species, sep = "_")) %>% 
  select(population, tl) %>% 
  add_count(population, name = "population_n")

# Species-level data
CBF_data_obs_species_binned <- 
  CBF_data_obs %>% 
  mutate(population = paste(species, sep = "_")) %>% 
  mutate(size_class = rls_bin(tl*10)) %>% 
  count(population, species, size_class) %>% 
  add_count(population, name = "n_sizebins") %>% 
  add_count(population, wt = n, name = "population_n")


# Location-level data
CBF_data_obs_location_binned <- 
  CBF_data_obs %>% 
  mutate(population = paste(location, species, sep = "_")) %>% 
  mutate(size_class = rls_bin(tl*10)) %>% 
  count(population, species, location, size_class) %>% 
  add_count(population, name = "n_sizebins") %>% 
  add_count(population, wt = n, name = "population_n")

# ==============================================================================
#  END
# ==============================================================================