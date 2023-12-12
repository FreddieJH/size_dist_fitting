library(tidyverse)

c(
  "input",
  "input/data",
  "input/data/raw", 
  "input/data/processed", 
  "input/data/data_cleaning", 
  "input/models",
  "output", 
  "output/tables",
  "output/model_fits", 
  "output/model_fits/rls",
  "output/model_fits/rls/species", 
  "output/model_fits/rls/ecoregion", 
  "output/model_fits/rls/gridcell", 
  "output/model_fits/cbf",
  "output/model_fits/cbf/species",
  "output/model_fits/cbf/location",
  "output/figures", 
  "output/figures/data_vis", 
  "output/figures/data_vis/empirical_dists", 
  "output/figures/data_vis/bimodal_dists", 
  "output/figures/ms_figs", 
  "scripts/"
) %>% 
  map(dir.create, showWarnings = FALSE)
