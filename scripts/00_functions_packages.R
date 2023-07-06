# Required packages ============================================================

library(tidyverse)
library(arrow)
library(rstan)
library(tidybayes, include.only = c("spread_draws"))
library(cowplot)
library(rlang)
library(posterior)
library(fitdistrplus, include.only = "fitdist")
library(patchwork)
library(scales)
library(magick)

rstan_options(auto_write = TRUE) # avoid recompilation of stan files

# Create required directories ==================================================

c(
  "input/",
  "input/RLS/",
  "input/CBF/",
  "input/RLS/data",
  "input/CBF/data",
  "input/RLS/models",
  "input/CBF/models",
  "output/",
  "output/RLS/",
  "output/CBF/",
  "output/RLS/model_fits",
  "output/RLS/model_plots",
  "output/RLS/data_figs",
  "output/RLS/data_figs",
  "output/RLS/likelihood_plots",
  "output/RLS/likelihood_values",
  "output/RLS/model_checks",
  "output/RLS/model_param_plots",
  "output/RLS/model_pars",
  "output/CBF/model_fits",
  "output/CBF/model_plots",
  "output/CBF/data_figs",
  "output/CBF/data_figs",
  "output/CBF/likelihood_plots",
  "output/CBF/likelihood_values",
  "output/CBF/model_checks",
  "output/CBF/model_param_plots",
  "output/CBF/model_pars",
  "output/CBF/traceplots"
  ) %>% 
  map(.f = function(dir) if(!dir.exists(dir)) dir.create(dir))


# Functions for binning into RLS bins ==========================================

rls_bin_breaks <- 
c(2.5, 5.0, 7.5,  10.0, 12.5, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 
  50.0, 62.5, 75.0, 87.5, 100.0, 112.5, 125.0, 137.5, 150.0, 
  162.5, 175.0, 187.5, 200.0, 250.0, 300.0, 350.0, 400.0)

rls_bin_table <-
  tibble(size_class = c(0, rls_bin_breaks, 500)) %>% 
  mutate(
    size_indx  = 0:(length(size_class)-1),
    size_min = (size_class + lag(size_class))/2,
    size_max = lead(size_min)
  ) %>% 
  filter(size_class %in% c(rls_bin_breaks, 500))


rls_bin <- function(size) {
  
  rls_bin_table$size_class[.bincode(size, breaks = c(0, rls_bin_table$size_max))]
}

# Data wrangling functions =====================================================

# outputs clean body size bin table given a vector of body sizes
get_bintable <- function(size_vector){
  
  sizebins <- 
    size_vector %>% 
    unique() %>% 
    c(0) %>%
    sort()
  
  tibble(size_class = sizebins) %>% 
    mutate(size_indx = 0:(nrow(.)-1),
           size_min = (size_class + lag(size_class))/2,
           size_max = lead(size_min)) %>% 
    filter(size_class != 0) %>% 
    mutate(size_max = case_when(size_class == max(sizebins) ~ size_class + (size_class-size_min), 
                                TRUE ~ size_max))
}

# function that takes population number, size_class, and n
clean_data <- function(count_table, 
                       sizes = "size_class", 
                       count = "n"){
  
  sizebin_tbl <- 
    count_table %>% 
    pull(!!sizes) %>% 
    get_bintable()
  
  popln_tbl <- 
    count_table %>% 
    select(population) %>% 
    distinct() %>% 
    mutate(population_indx = row_number())
  
  count_table %>% 
    rename(size_class := !!sizes) %>% 
    left_join(popln_tbl, by = join_by(population)) %>% 
    left_join(sizebin_tbl, by = join_by(size_class)) %>% 
    add_count(population_indx, wt = n, name = "population_n") %>% 
    arrange(population_indx, size_indx) %>%
    mutate(p = n/population_n) %>% 
    mutate(row = 1:n()) %>% 
    mutate(min_row = min(row), 
           max_row = max(row),
           .by = population_indx) 
  
}


# make the stan datalist for models 03 and 04
make_standata <- function(cleaned_data){
  
  poplvl_data <- 
    cleaned_data %>% 
    select(population_indx, 
           population_n, 
           min_row, 
           max_row) %>% 
    distinct() %>% 
    arrange(population_indx)
  
  sizebin_tbl <- 
    cleaned_data %>% 
    pull(size_class) %>% 
    get_bintable()
  
  
  mean_sizes <- 
    cleaned_data %>% 
    summarise(mean_size = mean(size_class, wt = n), 
              .by = population_indx) %>% 
    arrange(population_indx)
  
  list(
    I = nrow(cleaned_data),
    S = max(cleaned_data$population_indx),
    B = max(cleaned_data$size_indx),
    l = c(sizebin_tbl$size_min, sizebin_tbl$size_max[nrow(sizebin_tbl)]),
    # N_species = poplvl_data$population_n,
    i_min     = poplvl_data$min_row,
    i_max     = poplvl_data$max_row,
    s = cleaned_data$population_indx,
    b = cleaned_data$size_indx,
    n = cleaned_data$n,
    meansize = mean_sizes$mean_size
  )
}