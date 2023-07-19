# ==============================================================================
#  DATA IMPORT
# ==============================================================================

# RLS data import --------------------------------------------------------------

RLS_data_obs <- 
  "analysis/input/data/data_obs_cleaned.parquet" %>% 
  read_parquet()

# survey-level data
RLS_data_survs <- 
  "analysis/input/data/survey_list_m1_aus.parquet" %>% 
  read_parquet() %>% 
  mutate(lat_grid = round(latitude), 
         lon_grid = round(longitude), 
         lat_lon = paste(lat_grid, lon_grid, sep = "_")) 

# simulated data (from data_simulation.R)
RLS_data_sim <- 
  "analysis/input/data/data_simulated.parquet" %>% 
  read_parquet()

# Subset of real data
set.seed(1)
RLS_data_obs_sample <-
  RLS_data_obs %>% 
  filter(species_name %in% sample(unique(RLS_data_obs$species_name), 20))


# Cryptobenthic data import ----------------------------------------------------

# Raw Cryptobenthic fish (CBF) data 
CBF_raw <- 
  read_csv("analysis/input/data/crypto.size.data.csv", 
           show_col_types = FALSE)

# Observational level CBF data
CBF_data_obs <- 
  CBF_raw %>%
  mutate(sl = SL/10,
         tl = TL/10) %>% # now in cm, not mm
  select(location = Location,
         species = sciname,
         sl,
         tl,
         wt = W) %>%
  filter(str_detect(species, "^[A-Z]{1}[a-z]+\\s[a-z]+$")) 

# Species information from CBF data
CBF_spp_list <- 
  CBF_raw %>% 
  select(species = sciname, 
         genus = Genus, 
         family = Family) %>% 
  distinct() %>% 
  add_count(species, name = "n_byspp")

# Location-level information of CBF data
CBF_loc_list <- 
  CBF_raw %>% 
  select(location = Location, 
         lat = Lat, 
         lon = Long) %>%
  summarise(lat = mean(lat, na.rm = TRUE),
            lon = mean(lon, na.rm = TRUE), 
            .by = location) %>% 
  distinct()

rm(CBF_raw)

# ==============================================================================
#  END
# ==============================================================================