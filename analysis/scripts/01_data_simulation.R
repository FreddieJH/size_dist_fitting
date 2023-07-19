
if(!file.exists("analysis/input/data/data_simulated.parquet")){
  # set number of simulated species -----------------------------------------
  sim_n_spp <- 20
  sim_n_reps <- 1000
  set.seed(1)
  ## bin breaks ----
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
  
  is_odd <- function(x) ((x/2) %% 1) == 0.5
  
  data_out <- tibble()
  
  for(i in 1:sim_n_spp){
    
    # is odd number?
    if(is_odd(i)) {
      xx <- 
        tibble(num = 1:sim_n_reps) %>% 
        mutate(mean = rlnorm(1, 2.45, 0.78), 
               sd =  0.09687 + mean*0.33557) %>% 
        mutate(size_class_raw = rnorm(n = sim_n_reps, 
                                      mean = mean, 
                                      sd = sd)) %>% 
        mutate(species_indx = i) %>% 
        mutate(species_name = paste0("spp_", species_indx)) %>% 
        filter(size_class_raw > 1.25) %>% 
        mutate(size_class = rls_bin(size_class_raw)) %>% 
        dplyr::select(-num)
    } else {
      xx <- 
        tibble(num = 1:sim_n_reps) %>% 
        mutate(mean = rlnorm(1, 1.5, 0.1), 
               sd =  rlnorm(1, 0.5, 0.1)) %>% 
        mutate(size_class_raw = rlnorm(n = sim_n_reps, 
                                       meanlog = log(mean), 
                                       sdlog = log(sd))) %>% 
        mutate(species_indx = i) %>% 
        mutate(species_name = paste0("spp_", species_indx)) %>% 
        filter(size_class_raw > 1.25) %>% 
        mutate(size_class = rls_bin(size_class_raw)) %>% 
        dplyr::select(-num)
    }
    
    data_out <- bind_rows(data_out, xx)
  }
  
  data_sim <- 
    data_out %>% 
    left_join(rls_bin_table %>% dplyr::select(size_class, size_indx),
              by = join_by(size_class)) 
  
  write_parquet(data_sim, "analysis/input/data/data_simulated.parquet")
} else {
  data_sim <- read_parquet("analysis/input/data/data_simulated.parquet")
}
