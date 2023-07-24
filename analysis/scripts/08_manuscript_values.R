library(readr)
library(dplyr)

plot_data <- read_csv("ms_data/parameter_summary.csv", show_col_types = FALSE)
get_ci <- function(data_type, param, conf_level = 0.95, signif = 1) {
  
  vec <-
    plot_data %>% 
    filter(type == data_type) %>% 
    pull(param) %>% 
    na.omit(vec)
  sd <- sd(vec)
  n <- length(vec)
  se <- sd(vec)/sqrt(length(vec))
  xx <- qt(1 - ((1 - conf_level) / 2), n - 1) * se
  xx %>% 
    signif(signif)
}

get_mean <- function(data_type, param, signif = 2) {
  
  vec <-
    plot_data %>% 
    filter(type == data_type) %>% 
    pull(param) %>% 
    na.omit(vec)
  vec %>% 
    mean() %>% 
    signif(signif)
}

param_estimates <-
tibble(
  CBF_sdlog_mean = get_mean("CBF", "mean_sdlog"),
  CBF_sdlog_ci = get_ci("CBF", "mean_sdlog"),
  CBF_meanlog_mean = get_mean("CBF", "mean_meanlog"),
  CBF_meanlog_ci = get_ci("CBF", "mean_meanlog"),
  CBF_sigma_mean = get_mean("CBF", "mean_sigma"),
  CBF_sigma_ci = get_ci("CBF", "mean_sigma"),
  CBF_mu_mean = get_mean("CBF", "mean_mu"),
  CBF_mu_ci = get_ci("CBF", "mean_mu"),
  RLS_sdlog_mean = get_mean("RLS", "mean_sdlog"),
  RLS_sdlog_ci = get_ci("RLS", "mean_sdlog"),
  RLS_meanlog_mean = get_mean("RLS", "mean_meanlog"),
  RLS_meanlog_ci = get_ci("RLS", "mean_meanlog"),
  RLS_sigma_mean = get_mean("RLS", "mean_sigma"),
  RLS_sigma_ci = get_ci("RLS", "mean_sigma"),
  RLS_mu_mean = get_mean("RLS", "mean_mu"),
  RLS_mu_ci = get_ci("RLS", "mean_mu"),
)

write_csv(x = param_estimates, 
          file = "ms_data/param_estimates.csv")


# Table with the count of norm vs lnorm for each spatial scales
percent_norm_better <- tibble()
for(min_count in c(5000, 1000, 500, 200)){
  for(population_level in c("species", "gridcell", "ecoregion")){
    
    if(!(min_count %in% c(500, 200) & population_level == "gridcell")){
      
      
      RLS_stan_model <- "RLS_mod13"
      RLS_min_bins   <- 4
      RLS_min_count  <- min_count
      RLS_population_level = population_level
      
      
      RLS_obs_data <-
        paste0("RLS_data_obs_", RLS_population_level) %>%
        get() %>%
        filter(n_sizebins >= RLS_min_bins,
               population_n >= RLS_min_count) %>%
        clean_data()
      
      RLS_n_populations <-
        RLS_obs_data %>%
        pull(population) %>%
        n_distinct()
      
      RLS_model_name  <- paste0(RLS_stan_model,
                                "_nbin", RLS_min_bins,
                                "_n", RLS_min_count, "_",
                                RLS_population_level,
                                RLS_n_populations)
      
      norm_better_tbl <- 
        paste0("analysis/output/likelihood_values/", RLS_model_name, ".parquet") %>%
        read_parquet() %>% 
        pivot_wider(values_from = ll, 
                    names_from = dist) %>% 
        mutate(norm_better = N > LN) %>% 
        summarise(prop_norm_better = mean(norm_better),
                  .by = population) %>% 
        mutate(norm_better = prop_norm_better > 0.5)
      
      xx <- 
        norm_better_tbl %>% 
        count(norm_better) %>% 
        pivot_wider(names_from = norm_better, 
                    values_from = n) %>% 
        rename(norm_better = `TRUE`,
               lnorm_better = `FALSE`) %>% 
        mutate(model = RLS_stan_model, 
               min_count = RLS_min_count,
               population = RLS_population_level) %>% 
        mutate(percent_norm_better = (norm_better)/(norm_better + lnorm_better))
      
      percent_norm_better <- 
        percent_norm_better %>% 
        bind_rows(xx)
      
    }
  }
}

percent_norm_better %>% 
  ggplot(aes(x = population, 
             y = percent_norm_better, 
             fill = min_count %>% as.factor())) +
  geom_col(position = "dodge") +
  geom_hline(yintercept =  0.5)

