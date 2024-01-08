# About ========================================================================

# This script runs all the models from the analysis, exports the posterior 
# distributions for each parameter for each population and model. This script 
# also calculates the summary statistics for the Coefficient of Variation and 
# mean size for each population. 
# The script loops over the data-types (RLS, CBF), the spatial scales and the 
# two models (normal, lognormal).

# Packages =====================================================================

library(tidyverse)
library(scales)
library(cowplot)
library(rstan)
library(posterior)
library(stringr)
library(arrow)

# Parameters ===================================================================

force_run <- FALSE

# Model fitting ================================================================

# model fits are saved as compressed RDS files in the output/models/fits folder

for(dat in c("cbf", "rls")){
  for(pop in c("species", "ecoregion", "gridcell")){
    for(mod in c("normal", "lognormal")){
      
      if(dat == "cbf" & pop == "ecoregion") next
      pop <- ifelse(dat == "cbf" & pop == "gridcell", "location", pop)
      
      population_list <- 
        get(paste0("obsdata_", dat, "_", pop)) %>% 
        pull(population) %>% 
        unique() 
      
      res_table <- tibble()
      
      for(i in population_list) {
        
        filename <- 
          paste0("output/models/fits/", 
                 dat, "/",
                 pop, "/",
                 mod, "/",
                 i,
                 ".rds")
        
        if(!file.exists(filename)|force_run){
        
        # extract population data   
        current_data <-
          get(paste0("obsdata_", dat, "_", pop)) %>% 
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
              low_bound = ifelse(dat == "rls", 1.25, 0.05)
            )
          
            if(mod == "normal"){
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
            }
            if(mod == "lognormal"){
              # inits to prevent log(0) probability errors with the normal dist
              fit <- 
                stan(file = "input/models/lognormal.stan",
                     data = stan_data,
                     iter = 2000,
                     warmup = 1000,
                     chains = 3,
                     refresh = 0,
                     cores = 1, 
                     init = list(list(meanlog = inits$meanlog, sdlog = inits$sdlog),
                                 list(meanlog = inits$meanlog, sdlog = inits$sdlog),
                                 list(meanlog = inits$meanlog, sdlog = inits$sdlog)))
            }
            
            saveRDS(fit, 
                    file = filename, 
                    compress = TRUE)
            
            print(paste(filename, "saved."))
            
        }
        rm(current_data)
        }
      }
      rm(mod)
    }
    rm(pop)
  }
  rm(dat)
}

# Model diagnostics ============================================================

# diagnostics are tibbles with the Rhat and n_eff values for each of the parameter
# estimates for each population. They are saved as csv files in the 
# output/models/diagnostics folder

for(dat in c("rls", "cbf")){
  for(mod in c("lognormal", "normal")){
    for(pop in c("gridcell", "species", "ecoregion")){
      
      out_conv <- tibble()
      out_errs <- tibble()
      
      diag_filename <- 
        paste0("output/models/diagnostics/", 
               dat, "_", 
               pop, "_",
               mod, "_diags.csv")
      
      if(!file.exists(diag_filename)|force_run){
        
        files <- 
          list.files(paste("output/models/fits",
                           dat, 
                           pop,
                           mod,
                           sep = "/"),
                     full.names = TRUE)
        
        
        for(j in files){
          
          curr_pop <- str_extract(j, "([^\\//]+(?=.rds))")
          
          fit <- readRDS(j)
          
          if(!is.null(summary(fit)$summary[,"Rhat"])){
            
            rhat <- 
              summary(fit)$summary[,"Rhat"] %>% 
              as_tibble(rownames = "pars") %>% 
              mutate(population = curr_pop) %>% 
              rename(rhat = value)
            
            n_eff <- 
              summary(fit)$summary[,"n_eff"]%>% 
              as_tibble(rownames = "pars") %>% 
              mutate(population = curr_pop) %>% 
              rename(n_eff = value)
            
            out_conv <- 
              bind_rows(
                out_conv,
                left_join(rhat, n_eff, 
                          by = join_by(pars, population)))
            
            print(paste0(curr_pop, " added."))
            
          } else {
            
            out_errs <- 
              out_errs %>% 
              bind_rows(tibble(population = curr_pop))
            print(paste0(curr_pop, " had error."))
          } 
          
        }
        
        write_csv(out_conv,
                  diag_filename)
        
        write_csv(out_errs,
                  paste0("output/models/diagnostics/", 
                         dat, "_", 
                         pop, "_",
                         mod, "_errors.csv"))
        
        rm(out_conv, out_errs)
      } else {
        print(paste(diag_filename, "already exists."))
      }
      rm(pop)
    }
    rm(mod)
  }
  rm(dat)
}


# Model traceplots =============================================================

# Traceplots are plotted for every parameter and every population

# j = "output/models/fits/rls/species/normal/Caranx lugubris__Caranx lugubris.rds"
for(dat in c("rls", "cbf")){
  for(mod in c("normal", "lognormal")){
    for(pop in c("species", "ecoregion","gridcell" )){
      
      files <- 
        list.files(paste("output/models/fits",
                         dat, 
                         pop,
                         mod,
                         sep = "/"),
                   full.names = TRUE)
      
      count = 0
      for(j in files){
        count = count + 1
        curr_pop <- str_extract(j, "([^\\//]+(?=.rds))")
        
        trace_filename <- 
          paste0("output/models/traceplots/", 
                 dat, "/", 
                 pop, "/",
                 mod, "/", 
                 curr_pop, 
                 ".png")
        
        if(!file.exists(trace_filename)|force_run){
          fit <- readRDS(j) 
          if(!is.null(summary(fit)$summary[,"Rhat"])){
            
            p <- 
              fit %>% 
              traceplot() 
            
            ggsave(trace_filename, 
                   plot = p)
            
            print(paste0("(", count, ")", curr_pop, " traceplot saved."))
            
          } else {
            print(paste(curr_pop, "stan fit has errors (no traceplot saved)."))
          }
          rm(fit)
        } else {
          print(paste(curr_pop, "traceplot exists."))
        }
      }
      
      
      rm(pop)
    }
    rm(mod)
  }
  rm(dat)
}

# Posteriors ===================================================================

for(dat in c("rls", "cbf")){
  for(pop in c("gridcell", "species", "ecoregion")){
    for(mod in c("lognormal", "normal")){
      
      if(dat == "cbf" & pop == "ecoregion") next
      pop <- ifelse(dat == "cbf" & pop == "gridcell", "location", pop)
      
      files <- 
        list.files(paste("output/models/fits",
                         dat, 
                         pop,
                         mod,
                         sep = "/"),
                   full.names = TRUE)
      
      count = 0
      for(j in files){
        count = count + 1
        
        curr_pop <- str_extract(j, "([^\\//]+(?=.rds))")
        
        posteriors_filename <- 
          paste0("output/models/posteriors/", 
                 dat, "/", 
                 pop, "/",
                 mod, "/", 
                 curr_pop, 
                 ".parquet")
        
        
        
        if(!file.exists(posteriors_filename)|force_run){
          fit <- readRDS(j) 
          if(!is.null(summary(fit)$summary[,"Rhat"])){
            
            fit %>% 
              as_draws_df() %>% 
              as_tibble() %>% 
              write_parquet(posteriors_filename)
            
            print(paste0("(", count,") ", curr_pop, " posteriors saved."))
            
          } else {
            print(paste0("(", count,") ", curr_pop, " has an error."))
          }
          
          
          rm(fit)
        } else {
          # print(paste0("(", count, ") ",curr_pop, " posteriors previously saved."))
        }
        
      }
      
      posteriors_files <- 
        list.files(paste0("output/models/posteriors/", 
                          dat, "/", 
                          pop, "/",
                          mod, "/"), 
                   pattern = ".parquet", 
                   full.names = TRUE)
      
      summary_filename <- 
        paste0("output/models/summary/posteriors/", 
               dat, "_", 
               pop, "_",
               mod,  
               ".parquet")
      
      
      read_mutate <- function(file) {
        read_parquet(file) %>% 
          mutate(population = str_extract(file, "([^\\//]+(?=.parquet))"))
        }
      
      if(!file.exists(summary_filename)|force_run){
        bind_rows(
          lapply(
            posteriors_files, 
            read_mutate
          )
        ) %>% 
          write_parquet(summary_filename)
      }
      
      
      rm(mod, read_mutate)
    }
    rm(pop)
  }
  rm(dat)
}


# CV calculation ===============================================================

# Because we are dealing with truncated distributions (no individuals observed 
# less than the minimum observable body size) we need to correct for the true 
# mean and true variance of the distributions if we want to calculate the true CV

# 
# # calculate normalising constant for truncated normal
# CompSimp_TN_k <- function(I = step_size, l_min = low_bound, l_max = 500.0, 
#                           mu = mu, sigma = sigma) {
#   
#   h <- (l_max - l_min) / I
#   l <- l_min + h*(0:I)
#   f <-  dnorm(l, mean = mu, sd = sigma)
#   j <- 1:((I/2)-1)
#   sum_1 <- sum(f[2*j + 1])
#   j <- 1:(I/2)
#   sum_2 <- sum(f[2*j])
#   
#   return((h/3)*(f[1] + 2*sum_1 + 4*sum_2 + f[I+1]))
# }
# 
# calc_k = Vectorize(CompSimp_TN_k)
# 
# # calculate mean for the truncated normal
# CompSimp_TN_mean <- function(I = step_size, l_min = low_bound, l_max = 500.0, 
#                              mu = mu, sigma = sigma, k) {
#   
#   h <- (l_max - l_min) / I
#   l <- l_min + h*(0:I)
#   f <- l*dnorm(l, mean = mu, sd = sigma) / k
#   j <- 1:((I/2)-1)
#   sum_1 <- sum(f[2*j + 1])
#   j <- 1:(I/2)
#   sum_2 <- sum(f[2*j])
#   
#   return((h/3)*(f[1] + 2*sum_1 + 4*sum_2 + f[I+1]))
# }
# 
# calc_mean <- Vectorize(CompSimp_TN_mean)
# 
# # calculate variance for truncated normal
# CompSimp_TN_var <- function(I = step_size, l_min = low_bound, l_max = 500.0, 
#                             mu = mu, sigma = sigma, k, l_bar) {
#   
#   h <- (l_max - l_min) / I
#   l <- l_min + h*(0:I)
#   f <- ((l - l_bar)^2)*dnorm(l, mean = mu, sd = sigma) / k
#   j <- 1:((I/2)-1)
#   sum_1 <- sum(f[2*j + 1])
#   j <- 1:(I/2)
#   sum_2 <- sum(f[2*j])
#   
#   return((h/3)*(f[1] + 2*sum_1 + 4*sum_2 + f[I+1]))
# }
# 
# calc_var <- Vectorize(CompSimp_TN_var)
# 
# 
# # calculate coefficient of variation for truncated normal
# CompSimp_TN_cv <- function(l_bar, l_var) {
#   
#   return(sqrt(l_var)/l_bar)
# }



f_mean_norm <- function(x, mu, sigma, nc) (x*dnorm(x, mean = mu, sd = sigma))/nc
f_var_norm <- function(x, mu, sigma, nc, mean) (((x-mean)^2)*dnorm(x, mean = mu, sd = sigma))/nc
calc_nc_norm <- Vectorize(function(mu, sigma, lowbound) 1 - integrate(dnorm, lower = -Inf, upper = lowbound, mean = mu, sd = sigma, abs.tol = 0)$value)
calc_mean_norm <- Vectorize(function(mu, sigma, nc, lowbound) integrate(f_mean_norm, lower = lowbound, upper = Inf, mu = mu, sigma = sigma, nc = nc, abs.tol = 0)$value)
calc_var_norm <- Vectorize(function(mu, sigma, nc, mean, lowbound) integrate(f_var_norm, lower = lowbound, upper = Inf, mu = mu, sigma = sigma, nc = nc, mean = mean, abs.tol = 0)$value)
# 
# f_mean_lnorm <- function(x, meanlog, sdlog, nc) (x*dlnorm(x, meanlog = meanlog, sdlog = sdlog))/nc
# f_var_lnorm <- function(x, meanlog, sdlog, nc, mean) (((x-mean)^2)*dlnorm(x, meanlog = meanlog, sdlog = sdlog))/nc
# calc_nc_lnorm <- Vectorize(function(meanlog, sdlog, lowbound) 1 - integrate(dlnorm, lower = -Inf, upper = lowbound, meanlog = meanlog, sdlog = sdlog, abs.tol = 0, stop.on.error = FALSE)$value)
# calc_mean_lnorm <- Vectorize(function(meanlog, sdlog, nc, lowbound) integrate(f_mean_lnorm, lower = lowbound, upper = 1000, meanlog = meanlog, sdlog = sdlog, nc = nc, abs.tol = 0, stop.on.error = FALSE)$value)
# calc_var_lnorm <- Vectorize(function(meanlog, sdlog, nc, mean, lowbound) integrate(f_var_lnorm, lower = lowbound, upper = 1000, meanlog = meanlog, sdlog = sdlog, nc = nc, mean = mean, abs.tol = 0, stop.on.error = FALSE)$value)



for(dat in c("rls", "cbf")){
  for(pop in c("gridcell", "species", "ecoregion")){
    for(mod in c("lognormal", "normal")){
      
      
      if(dat == "cbf" & pop == "ecoregion") next
      pop <- ifelse(dat == "cbf" & pop == "gridcell", "location", pop)
      
      lower_bound <- ifelse(dat == "rls", 1.25, 0.05)
      
      
      input_filename <- 
        paste0("output/models/summary/posteriors/", 
               dat, "_", 
               pop, "_",
               mod,  
               ".parquet")
      
      output_filename <- 
        paste0("output/models/summary/cv/", 
               dat, "_", 
               pop, "_",
               mod,  
               ".parquet")
      
      if(!file.exists(output_filename) | force_run){
        
        if(mod == "normal"){
          input_filename %>% 
            read_parquet() %>% 
            mutate(nc = calc_nc_norm(lowbound = lower_bound, 
                                     mu = mu,
                                     sigma = sigma),
                   mean = calc_mean_norm(lowbound = lower_bound,
                                         mu = mu,
                                         sigma = sigma,
                                         nc = nc),
                   var = calc_var_norm(lowbound = lower_bound,
                                       mu = mu,
                                       sigma = sigma,
                                       nc = nc,
                                       mean = mean),
                   cv = sqrt(var)/mean) %>%
            summarise(mean_05 = quantile(mean, probs = 0.05),
                      mean_10 = quantile(mean, probs = 0.1),
                      mean_50 = quantile(mean, probs = 0.5),
                      mean_90 = quantile(mean, probs = 0.9),
                      mean_95 = quantile(mean, probs = 0.95),
                      cv_05 = quantile(cv, probs = 0.05),
                      cv_10 = quantile(cv, probs = 0.1),
                      cv_50 = quantile(cv, probs = 0.5),
                      cv_90 = quantile(cv, probs = 0.9),
                      cv_95 = quantile(cv, probs = 0.95),
                      logl_50 = quantile(lp__, prob = 0.5), 
                      .by = population) %>%
            write_parquet(output_filename)
          
        }
        if(mod == "lognormal"){
          
          input_filename %>% 
            read_parquet() %>%
            mutate(mean =  exp(meanlog + 0.5*(sdlog)^2),
                   var =  (exp((sdlog)^2)- 1)*exp(2*meanlog + (sdlog)^2),
                   cv = sqrt(var) / mean) %>%
            summarise(mean_05 = quantile(mean, probs = 0.05),
                      mean_10 = quantile(mean, probs = 0.1),
                      mean_50 = quantile(mean, probs = 0.5),
                      mean_90 = quantile(mean, probs = 0.9),
                      mean_95 = quantile(mean, probs = 0.95),
                      cv_05 = quantile(cv, probs = 0.05),
                      cv_10 = quantile(cv, probs = 0.1),
                      cv_50 = quantile(cv, probs = 0.5),
                      cv_90 = quantile(cv, probs = 0.9),
                      cv_95 = quantile(cv, probs = 0.95),
                      logl_50 = quantile(lp__, prob = 0.5), 
                      .by = population) %>%
            write_parquet(output_filename)
        }
      } else {
        print(paste0(output_filename, " CV previously saved."))
      }
      
      
      
      rm(mod)
    }
    rm(pop)
  }
  rm(dat)
}


# Summarise parameters =========================================================

for(dat in c("rls", "cbf")){
  for(pop in c("gridcell", "species", "ecoregion")){
    for(mod in c("lognormal", "normal")){
      
      if(dat == "cbf" & pop == "ecoregion") next
      pop <- ifelse(dat == "cbf" & pop == "gridcell", "location", pop)
      
      lower_bound <- ifelse(dat == "rls", 1.25, 0.05)
      
      input_filename <- 
        paste0("output/models/summary/posteriors/", 
               dat, "_", 
               pop, "_",
               mod,  
               ".parquet")
      
      output_filename <- 
        paste0("output/models/summary/pars/", 
               dat, "_", 
               pop, "_",
               mod,  
               ".parquet")
      
      if(!file.exists(output_filename)|force_run){
      
      if(mod == "normal"){
        input_filename %>% 
          read_parquet() %>% 
          summarise(mu_05 = quantile(mu, probs = 0.05),
                    mu_10 = quantile(mu, probs = 0.1),
                    mu_50 = quantile(mu, probs = 0.5),
                    mu_90 = quantile(mu, probs = 0.9),
                    mu_95 = quantile(mu, probs = 0.95),
                    sigma_05 = quantile(sigma, probs = 0.05),
                    sigma_10 = quantile(sigma, probs = 0.1),
                    sigma_50 = quantile(sigma, probs = 0.5),
                    sigma_90 = quantile(sigma, probs = 0.9),
                    sigma_95 = quantile(sigma, probs = 0.95),
                    logl_50 = quantile(lp__, prob = 0.5), 
                    .by = population) %>%
          write_parquet(output_filename)
        
      }
      if(mod == "lognormal"){
        
        input_filename %>% 
          read_parquet() %>%
          summarise(meanlog_05 = quantile(meanlog, probs = 0.05),
                    meanlog_10 = quantile(meanlog, probs = 0.1),
                    meanlog_50 = quantile(meanlog, probs = 0.5),
                    meanlog_90 = quantile(meanlog, probs = 0.9),
                    meanlog_95 = quantile(meanlog, probs = 0.95),
                    sdlog_05 = quantile(sdlog, probs = 0.05),
                    sdlog_10 = quantile(sdlog, probs = 0.1),
                    sdlog_50 = quantile(sdlog, probs = 0.5),
                    sdlog_90 = quantile(sdlog, probs = 0.9),
                    sdlog_95 = quantile(sdlog, probs = 0.95),
                    logl_50 = quantile(lp__, prob = 0.5), 
                    .by = population) %>%
          write_parquet(output_filename)
      }
      
      } else {
        print(paste(output_filename, "already saved."))
      }
      
      rm(mod)
    }
    rm(pop)
  }
  rm(dat)
}


# End ==========================================================================
