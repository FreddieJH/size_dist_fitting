# ==============================================================================
# BAYESIAN MODELLING
# ==============================================================================

# RLS modelling ----------------------------------------------------------------

# min_count = the minimum allowable individuals per population
# level = spatial scale of the data (choices: "species", "gridcell", "ecoregion")
for(min_count in c(5000, 1000, 500, 200)){
  for(population_level in c("species", "gridcell", "ecoregion")){
    
    if(!(min_count %in% c(500, 200) & population_level == "gridcell")){
      
      stan_model       <- "RLS_mod13"    # which stan model is being used
      min_bins         <- 4           # the minimum allowable bins per population
      stan_iter        <- 1e4         # number of stan iterations
      
      obs_data <-
        paste0("RLS_data_obs_", population_level) %>%
        get() %>%
        filter(n_sizebins >= min_bins,
               population_n >= min_count) %>%
        clean_data()
      
      n_populations <-
        obs_data %>%
        pull(population) %>%
        n_distinct()
      
      stan_data <-
        obs_data %>%
        make_standata()
      
      model_name  <- paste0(stan_model,
                            "_nbin", min_bins,
                            "_n", min_count, "_",
                            population_level,
                            n_populations)
      
      run_model(model_name = model_name,
                stan_model = stan_model,
                overwrite = FALSE,
                iter = stan_iter)
      extract_pars(model_name = model_name, overwrite = FALSE)
      run_traceplots(model_name = model_name, overwrite = FALSE)
      run_ess_check(model_name = model_name, overwrite = FALSE)
      run_model_vis(model_name = model_name, 
                    obs_data = obs_data, 
                    overwrite = FALSE)
      run_ll(model_name = model_name, overwrite = FALSE)
      run_ll_plot(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
      run_mean_vs_ll(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
      run_par_regression(model_name = model_name, overwrite = FALSE)
      
    }
  }
}

# CBF modelling (binned) -------------------------------------------------------

# min_count = the minimum allowable individuals per population
# level = spatial scale of the data (choices: "species", "location")
for(min_count in c(500, 100, 50, 20)){
  for(population_level in c("species", "location")){
    
    stan_model       <- "CBF_mod13"    # which stan model is being used
    min_bins         <- 4           # the minimum allowable bins per population
    stan_iter        <- 1e4         # number of stan iterations
    
    obs_data <-
      paste0("CBF_data_obs_", population_level, "_binned") %>%
      get() %>%
      filter(n_sizebins >= min_bins,
             population_n >= min_count) %>%
      clean_data()
    
    n_populations <-
      obs_data %>%
      pull(population) %>%
      n_distinct()
    
    population_indx_table <-
      obs_data %>%
      select(population, population_indx) %>%
      distinct()
    
    mod_data <-
      obs_data %>%
      left_join(population_indx_table,
                by = join_by(population))
    
    stan_data <-
      obs_data %>%
      make_standata()
    
    model_name  <- paste0(stan_model,
                          "_nbin", min_bins,
                          "_n", min_count, "_",
                          population_level,
                          n_populations)
    
    run_model(model_name = model_name, stan_model = stan_model, overwrite = FALSE, iter = stan_iter)
    extract_pars(model_name = model_name, overwrite = FALSE)
    run_traceplots(model_name = model_name, overwrite = FALSE)
    run_ess_check(model_name = model_name, overwrite = FALSE)
    run_model_vis(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
    run_ll(model_name = model_name, overwrite = FALSE)
    run_ll_plot(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
    run_mean_vs_ll(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
    run_par_regression(model_name = model_name, overwrite = FALSE)
  }
}



# CBF modelling (continuous) -------------------------------------------------------



# min_count = the minimum allowable individuals per population
# level = spatial scale of the data (choices: "species", "location")
for(min_count in c(500, 100, 50, 20, 10)){
  for(population_level in c("location", "species")){
    
    stan_model       <- "CBF_mod02"    # which stan model is being used
    stan_iter        <- 1e4         # number of stan iterations
    
    obs_data_step1 <-
      paste0("CBF_data_obs_", population_level) %>%
      get() %>%
      filter(population_n >= min_count) %>%
      mutate(tl = tl/10)
    
    n_populations <-
      obs_data_step1 %>%
      pull(population) %>%
      n_distinct()
    
    population_indx_table <-
      obs_data_step1 %>%
      select(population) %>%
      distinct() %>%
      mutate(population_indx = row_number())
    
    obs_data <-
      obs_data_step1 %>%
      left_join(population_indx_table,
                by = join_by(population))
    
    stan_data <-
      list(
        n = nrow(obs_data),
        n_population = max(obs_data$population_indx),
        population_id = obs_data$population_indx,
        y = obs_data$sl
      )
    
    model_name  <- paste0(stan_model,
                          "_n", min_count, "_",
                          population_level,
                          n_populations)
    
    run_model(model_name = model_name, stan_model = stan_model, overwrite = FALSE, iter = stan_iter)
    extract_pars(model_name = model_name, overwrite = FALSE)
    run_traceplots(model_name = model_name, overwrite = FALSE)
    run_ess_check(model_name = model_name, overwrite = FALSE)
    run_model_vis_continuous(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
    run_ll_continuous(model_name = model_name, overwrite = FALSE)
    run_ll_plot(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
    run_par_regression(model_name = model_name, overwrite = FALSE)
    run_mean_vs_ll_continuous(model_name = model_name, obs_data = obs_data, overwrite = FALSE)
  }
}

# ==============================================================================
#  END
# ==============================================================================