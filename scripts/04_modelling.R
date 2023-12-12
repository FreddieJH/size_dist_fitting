# About ========================================================================

# This script creates a plot with all the empirical body size distributions 

# Packages =====================================================================

library(tidyverse)
library(scales)
library(cowplot)
library(rstan)
library(posterior)

# Parameters ====================================================================

# how many bins for the integration process?
n_steps <- 100
force_run <- FALSE

# RLS model fitting ============================================================
src = "rls"
for(pop in c("species", "ecoregion", "gridcell")){
  
  population_list <- 
    get(paste0("obsdata_", src, "_", pop)) %>% 
    pull(population) %>% 
    unique() 
  
  res_table <- tibble()
  
  for(i in population_list) {
    
    # extract population data   
    current_data <-
      get(paste0("obsdata_", src, "_", pop)) %>% 
      filter(population == i) %>% 
      arrange(size_class)
    
    p_firstbin <- 
      current_data %>% 
      pull(p_obs) %>% 
      .[1]
    
    if (p_firstbin < 0.5) {
      
      # cat(paste(i, " normal\n", sep = "")) # used for debugging
      inits <-
        current_data %>% 
        uncount(n) %>% 
        summarise(mu = mean(size_class), 
                  sigma = sd(size_class),
                  meanlog = mean(log(size_class)), 
                  sdlog = sd(log(size_class)))
      
      stan_data <- 
        list(
          B = length(unique(current_data$size_class)),
          N = current_data$population_n[1], # total population size
          b_upr = current_data$size_max,
          b_lwr = current_data$size_min,
          n = current_data$n,               # bin abundances (sum = N)
          low_bound = 1.25
        )
      
      
      filename <- 
        paste0("output/model_fits/", 
               src, "/",
               pop, "/",
               i,
               ".rds")
      if(!file.exists(filename)|force_run){
        
        # inits to prevent log(0) probability errors with the normal dist
        fit <- 
          stan(file = "input/models/normal.stan",
               data = stan_data,
               iter = 2000,
               warmup = 1000,
               chains = 3,
               refresh = 0,
               cores = 1, 
               init = list(list(mu = inits$mu, sigma = inits$sigma),
                           list(mu = inits$mu, sigma = inits$sigma),
                           list(mu = inits$mu, sigma = inits$sigma)))
        
        
        saveRDS(fit, 
                file = filename)
        
      }
      # res_table <-
      #   res_table %>%
      #   bind_rows(
      #     fit %>%
      #       as_draws_df() %>%
      #       as_tibble() %>%
      #       mutate(population = i) %>%
      #       mutate(norm_const = calc_nc(mu = mu, sigma = sigma), 
      #              mean = calc_mean2(mu = mu, sigma = sigma, nc = norm_const), 
      #              var = calc_var2(mu = mu, sigma = sigma, nc = norm_const, mean = mean), 
      #              cv = sqrt(var)/mean)
      #   ) 
      
    }
  }
}


# CBF model fitting ============================================================


f <- function(x, mu, sigma, nc) (x*dnorm(x, mean = mu, sd = sigma))/nc
f2 <- function(x, mu, sigma, nc, mean) (((x-mean)^2)*dnorm(x, mean = mu, sd = sigma))/nc
calc_nc <- Vectorize(function(mu, sigma) integrate(dnorm, lower = 1.25, upper = Inf, mean = mu, sd = sigma, abs.tol = 0)$value)
calc_mean2 <- Vectorize(function(mu, sigma, nc) integrate(f, lower = 1.25, upper = Inf, mu = mu, sigma = sigma, nc = nc, abs.tol = 0)$value)
calc_var2 <- Vectorize(function(mu, sigma, nc, mean) integrate(f2, lower = 1.25, upper = Inf, mu = mu, sigma = sigma, nc = nc, mean = mean, abs.tol = 0)$value)

traceplot(fit)