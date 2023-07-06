
# Run Stan model ===============================================================

run_model <- function(model_name, stan_model, overwrite, iter = 5e3){
  
  if(!file.exists(paste0("output/model_fits/", model_name, ".rds")) | overwrite){
    
    model_fit <-
      stan(file = paste0("input/models/", stan_model, ".stan"),
           data = stan_data,
           iter = iter,
           warmup = iter*0.9,
           chains = 3,
           refresh = 100,
           init_r = 1,
           seed = 1,
           cores = 3
      )
    
    if(!dir.exists("output/model_fits/")) dir.create("output/model_fits/")
    
    write_rds(x = model_fit, 
              file = paste0("output/model_fits/", model_name, ".rds"))
    
    cat(paste(model_name, ": Model run.\n"))
  } 
}


# Extracts Parameters ==========================================================

extract_pars <- function(model_name, overwrite){
  
  if(!file.exists(paste0("output/model_pars/", model_name, ".csv")) | overwrite){
    
    param_table <-
      paste0("output/model_fits/", model_name, ".rds") %>% 
      read_rds() %>% 
      as_draws_df() %>% 
      pivot_longer(cols = -one_of(c('.chain', '.iteration', '.draw')), names_to = "variable") %>% 
      mutate(
        param = ifelse(str_detect(variable, "\\["), 
                       str_extract(variable, ".*(?=\\[)"), 
                       variable),
        population_indx = ifelse(str_detect(variable, "\\["), 
                                 str_extract(variable, "(?<=\\[).*(?=\\])") %>% 
                                   as.numeric(), 
                                 NA) 
      ) %>% 
      filter(!str_detect(variable, "ln_|lp__|logit_")) %>% 
      select(.draw, 
             param, 
             value, population_indx) 
    
    
    param_table_summary <- 
      paste0("output/model_fits/", model_name, ".rds") %>% 
      read_rds() %>% 
      summarise_draws() %>% 
      mutate(
        param = ifelse(str_detect(variable, "\\["), 
                       str_extract(variable, ".*(?=\\[)"), 
                       variable),
        population_indx = ifelse(str_detect(variable, "\\["), 
                                 str_extract(variable, "(?<=\\[).*(?=\\])") %>% 
                                   as.numeric(), 
                                 NA) 
      ) %>% 
      filter(!str_detect(variable, "ln_|lp__|logit_")) %>% 
      select(population_indx, 
             param,
             mean,
             q5, 
             q95, rhat, ess_bulk) %>% 
      pivot_wider(values_from = c(mean, q5, q95, rhat, ess_bulk), 
                  names_from  = param)
    
    write_parquet(x = param_table, 
                  sink = paste0("output/model_pars/", model_name, ".parquet"))
    
    write_parquet(x = param_table_summary, 
                  sink = paste0("output/model_pars/", model_name, "_summary.parquet"))
    
    cat(paste(model_name, ": Parameters extracted.\n")) 
  }
  
}


# Save Traceplots  =============================================================
# Excludes untransformed parameters (log, logit etc.)
# Will limit to 20 parameters per plot
run_traceplots <- function(model_name){
  
  model_par_list <-
    paste0("output/model_pars/", model_name, ".parquet") %>% 
    read_parquet() %>% 
    count(param, population_indx) %>% 
    count(param)
  
  model_fit <- 
    paste0("output/model_fits/", model_name, ".rds") %>% 
    read_rds()
  
  for(i in 1:nrow(model_par_list)){
    
    if(model_par_list$n[i] == 1){
      param_name  <- model_par_list$param[i]
    } else {
      param_name <- paste0(model_par_list$param[i], "[", 1:20, "]")
    }
    
    p <-
      traceplot(model_fit, param_name, inc_warmup = FALSE) 
    
    if(!dir.exists("output/traceplots/")){
      dir.create("output/traceplots/")
    }
    if(!dir.exists(paste0("output/traceplots/", model_name, "/"))){
      dir.create(paste0("output/traceplots/", model_name, "/"))
    }
    
    ggsave(filename = paste0("output/traceplots/", model_name, "/", model_par_list$param[i],".png"), 
           plot = p,
           height = 20, 
           width = 20*1.618, 
           units = "cm", 
           dpi = 96)
    
    cat(paste(model_name, ": Traceplot for", model_par_list$param[i], "saved.\n"))  
  }
  
}


# Save Diagnostic plots ========================================================
# Check for model convergence (using Rhat and ESS)
run_ess_check <- function(model_name){
  
  rhat_values <-
    paste0("output/model_pars/", model_name, "_summary.parquet") %>% 
    read_parquet() %>% 
    select(population_indx, starts_with("rhat_")) %>% 
    pivot_longer(cols = starts_with("rhat_"), names_prefix = "rhat_", names_to = "param", values_to = "rhat")  
  
  ess_values <- 
    paste0("output/model_pars/", model_name, "_summary.parquet") %>% 
    read_parquet() %>% 
    select(population_indx,  starts_with("ess_")) %>% 
    pivot_longer(cols = starts_with("ess_"), names_prefix = "ess_bulk_", names_to = "param", values_to = "ess") 
  
  
  
  p <- 
    left_join(rhat_values, ess_values, by = join_by(population_indx, param)) %>% 
    arrange(desc(rhat)) %>%  
    mutate(lab = paste(param, population_indx, sep = "_")) %>% 
    ggplot(aes(x = rhat, 
               y = ess)) + 
    geom_text(aes(label = lab)) +
    geom_vline(xintercept = 1.1, col = "red", lty = 2) +
    theme_cowplot() +
    ggtitle("Model checks, anything with a high rhat (> 1.1) value may be a problem")
  
  ggsave(filename = paste0("output/model_checks/", model_name, ".png"), 
         plot = p,
         height = 20, 
         width = 20*1.618, 
         units = "cm", 
         dpi = 96)
  
  cat(paste(model_name, ": Model diagnostic plot saved.\n"))
  
}


# Visualise model predictions ==================================================
# Visualise model fit with error ribbons
run_model_vis <- function(model_name, obs_data){
  
  plot <-
    paste0("output/model_pars/", model_name, ".parquet") %>% 
    read_parquet() %>% 
    filter(population_indx %in% 1:20) %>% 
    pivot_wider(names_from = param, values_from = value) %>% 
    full_join(obs_data %>% filter(population_indx %in% 1:20), 
              by = join_by(population_indx),
              multiple = "all") %>% 
    mutate(p_N = pnorm(size_max, mu, sigma) -  pnorm(size_min, mu, sigma), 
           p_LN = plnorm(size_max, meanlog, sdlog) -  plnorm(size_min, meanlog, sdlog)) %>% 
    summarise(mean_p_N = mean(p_N, na.rm = TRUE), 
              median_p_N = median(p_N, na.rm = TRUE),
              mean_p_LN = mean(p_LN, na.rm = TRUE), 
              median_p_LN = median(p_LN, na.rm = TRUE), 
              q5_p_N = quantile(p_N, probs = 0.05, na.rm = TRUE), 
              q5_p_LN = quantile(p_LN, probs = 0.05, na.rm = TRUE), 
              q95_p_N = quantile(p_N, probs = 0.95, na.rm = TRUE), 
              q95_p_LN = quantile(p_LN, probs = 0.95, na.rm = TRUE), 
              .by = c("population_indx", "population", "size_class", "p")) %>% 
    ggplot() +
    aes(x = size_class, 
        y = p) +
    geom_point() +
    geom_line() + 
    geom_point(aes(y = mean_p_N), 
               col = "blue",
               alpha = 0.5) +
    geom_ribbon(aes(ymin = q5_p_N,
                    ymax = q95_p_N),
                fill = "blue",
                alpha = 0.5) +
    geom_line(aes(y = mean_p_N), 
              col = "blue",
              alpha = 0.5) +
    geom_point(aes(y = mean_p_LN), 
               col = "red",
               alpha = 0.5) +
    geom_ribbon(aes(ymin = q5_p_LN,
                    ymax = q95_p_LN),
                fill = "red",
                alpha = 0.5) +
    geom_line(aes(y = mean_p_LN), 
              col = "red",
              alpha = 0.5) +
    facet_wrap(~population, 
               scales = "free")  +
    ggtitle(label = "Normal = blue, Lognormal = red")
  
  ggsave(filename = paste0("output/model_plots/", model_name, ".png"), 
         plot = plot,
         height = 20, 
         width = 20*1.618, 
         units = "cm", 
         dpi = 96)
  
  cat(paste(model_name, ": Model prediction plot saved.\n"))
  
}

# Parameter Regression =========================================================

# Plot regression of mean distribution parameter with variance distribution parameter
run_par_regression <- function(model_name){
  
  model_par_raw <-
    paste0("output/model_pars/", model_name, ".parquet") %>% 
    read_parquet() %>% 
    # filter(population_indx %in% 1:20) %>% 
    summarise(value = mean(value),
              .by = c(population_indx, param)) %>% 
    pivot_wider(values_from = value, names_from = param, names_prefix = "mean_")
  
  
  plot_pars <-  function(data = model_par_raw, 
                         x = "mean_mu", 
                         y = "mean_sigma"){
    
    coefs <- 
      lm(as.formula(paste0(y, "~", x)), data = data) %>% 
      coef()
    
    
    data %>% 
      ggplot() +
      aes(
        x = !!sym(x), 
        y = !!sym(y)
      ) +
      geom_point()  +
      stat_smooth(method = "lm", col = "orange", formula = 'y ~ x') +
      ggtitle(paste0(y, " = ", round(coefs[1], 2), " + ", round(coefs[2], 2), "*", x))
  }
  
  p1 <- plot_pars(x = "mean_mu", y = "mean_sigma")
  p2 <-  plot_pars(x = "mean_mu", y = "mean_cv")
  p3 <-  plot_pars(x = "mean_meanlog", y = "mean_sdlog")
  
  p4 <- p1 + p2 + p3
  
  ggsave(filename = paste0("output/model_param_plots/", model_name, ".png"), 
         plot = p4,
         height = 10, 
         width = 20*1.618, 
         units = "cm", 
         dpi = 96)
  
  cat(paste(model_name, ": Model parameter regression plot saved.\n"))
}


# Calculate loglikelihood (mod_13) =============================================

run_ll <- function(model_name, overwrite){
  if(!file.exists(paste0("output/likelihood_values/", model_name, ".parquet")) | overwrite){
    
    f_table <- 
      tibble(size_indx = 1:stan_data$B) %>% 
      mutate(f = (stan_data$l[size_indx+1] - stan_data$l[size_indx]) / (stan_data$l[max(size_indx)+1] - stan_data$l[1]))
    
    likelihood_table <- 
      paste0("output/model_pars/", model_name, ".parquet") %>% 
      read_parquet() %>% 
      pivot_wider(names_from = param, values_from = value) %>% 
      full_join(obs_data, 
                by = join_by(population_indx),
                multiple = "all") %>% 
      mutate( norm_c_N = 1.0 - pnorm(1.25, mu, sigma), 
              norm_c_LN = 1.0 - plnorm(1.25, mu, sigma)) %>% 
      mutate(
        p_N = (pnorm(size_max, mu, sigma) - pnorm(size_min, mu, sigma)) / norm_c_N, 
        p_LN = (plnorm(size_max, meanlog, sdlog) - plnorm(size_min, meanlog, sdlog)) / norm_c_LN, 
        p_N_adjust = (1.0 - eps_N)*p_N + eps_N*f_table$f[size_min], 
        p_LN_adjust = (1.0 - eps_LN)*p_LN + eps_LN*f_table$f[size_min],
        ll_N = n*log(p_N_adjust),
        ll_LN = n*log(p_LN_adjust)
      ) %>% 
      summarise(ll_N = sum(ll_N, na.rm = TRUE), 
                ll_LN = sum(ll_LN, na.rm = TRUE), 
                .by = c(.draw, population_indx, population)) %>% 
      pivot_longer(cols = starts_with("ll_"), names_to = "dist", values_to = "ll", names_prefix = "ll_") 
    
    if(!dir.exists("output/likelihood_values/")) dir.create("output/likelihood_values/")
    
    write_parquet(likelihood_table, 
                  paste0("output/likelihood_values/", model_name, ".parquet"))
  }
  
}

# Loglikelihood boxplots =======================================================

run_ll_plot <- function(model_name, obs_data){
  
  p <- 
    paste0("output/likelihood_values/", model_name, ".parquet") %>% 
    read_parquet() %>% 
    filter(population_indx %in% 1:20) %>% 
    mutate(dist = case_when(dist == "N" ~ "Normal", 
                            dist == "LN" ~ "Lognormal")) %>% 
    ggplot() +
    aes(x = dist, 
        y = ll, 
        col = dist) +
    geom_boxplot() +
    facet_wrap(~population, scales = "free") +
    scale_color_manual(values=c("Lognormal" = "red", 
                                "Normal" = "blue"))
  
  if(!dir.exists("output/likelihood_plots/")) dir.create("output/likelihood_plots/")
  
  ggsave(filename = paste0("output/likelihood_plots/boxplot_", model_name, ".png"), 
         plot = p,
         height = 20, 
         width = 20*1.618, 
         units = "cm", 
         dpi = 96)
}

# Loglikelihood vs meansize ====================================================

run_mean_vs_ll <- function(model_name, obs_data){
  
  popln_count <- 
    obs_data %>% 
    select(population_indx, 
           population_n) %>% 
    distinct()
  
  p_data <-   
    paste0("output/likelihood_values/", model_name, ".parquet") %>% 
    read_parquet() %>% 
    left_join(popln_count,
              by = join_by(population_indx)) %>% 
    mutate(ll = ll/population_n) %>% 
    pivot_wider(names_from = dist, 
                values_from = ll) %>% 
    mutate(diff = N - LN) %>% 
    summarise(diff = median(diff), 
              ll_norm = median(N), 
              ll_lnorm = median(LN),
              .by = population_indx) %>%
    left_join(summarise(obs_data, 
                        mean_size = mean(size_class, wt = n), 
                        .by = population_indx),
              by = join_by(population_indx)) 
  
  norm_better_n <- 
    p_data %>% 
    mutate(norm_better = diff > 0) %>% 
    count(norm_better) %>% 
    filter(norm_better) %>% 
    pull(n)
  
  p <- 
    p_data %>% 
    ggplot() +
    aes(
      x = mean_size,
      y = diff
    ) + 
    geom_rect(aes(ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.5) +
    geom_rect(aes(ymin = 0, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "blue", alpha = 0.5) +
    geom_point(size = 2) +
    annotate(geom = "text", y = 0, x = Inf, label = paste0("n = ", norm_better_n), 
             hjust = 1.5, vjust = -0.5, col = "white", size = 5) +
    annotate(geom = "text", y = 0, x = Inf, label = paste0("n = ", nrow(p_data)-norm_better_n), 
             hjust = 1.5, vjust = 1.5, col = "white", size = 5) +
    stat_smooth(formula = 'y ~ x', col = "orange", method = "lm") +
    labs(x = "Species mean size", 
         y = "Normalised loglikelighood ratio (norm_LL_N-norm_LL_LN)") +
    theme_cowplot(20)
  
  ggsave(filename = paste0("output/likelihood_plots/sizeVSll_", model_name, ".png"), 
         plot = p,
         height = 20, 
         width = 20*1.618, 
         units = "cm", 
         dpi = 96)
  p_data2 <-   
    paste0("output/likelihood_values/", model_name, ".parquet") %>% 
    read_parquet() %>% 
    pivot_wider(names_from = dist, 
                values_from = ll) %>% 
    mutate(norm_better = N>LN) %>% 
    summarise(prop_norm_better = mean(norm_better), 
              .by = population_indx) %>%
    left_join(summarise(obs_data, 
                        mean_size = mean(size_class, wt = n), 
                        .by = population_indx),
              by = join_by(population_indx)) 
  
  norm_better_n2 <- 
    p_data2 %>% 
    mutate(norm_better = prop_norm_better > 0.5) %>% 
    count(norm_better) %>% 
    filter(norm_better) %>% 
    pull(n)
  
  p2 <- 
    p_data2 %>% 
    ggplot() +
    aes(
      x = mean_size,
      y = prop_norm_better
    ) + 
    geom_rect(aes(ymin = -Inf, ymax = 0.5, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.5) +
    geom_rect(aes(ymin = 0.5, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "blue", alpha = 0.5) +
    geom_point(size = 2) +
    annotate(geom = "text", y = 0.5, x = Inf, label = paste0("n = ", norm_better_n), 
             hjust = 1.5, vjust = -0.5, col = "white", size = 5) +
    annotate(geom = "text", y = 0.5, x = Inf, label = paste0("n = ", nrow(p_data)-norm_better_n), 
             hjust = 1.5, vjust = 1.5, col = "white", size = 5) +
    stat_smooth(formula = 'y ~ x', col = "orange", method = "lm") +
    labs(x = "Species mean size", 
         y = "Normalised loglikelighood ratio (norm_LL_N-norm_LL_LN)") +
    theme_cowplot(20)
  
  ggsave(filename = paste0("output/likelihood_plots/sizeVSll_", model_name, "_method2.png"), 
         plot = p2,
         height = 20, 
         width = 20*1.618, 
         units = "cm", 
         dpi = 96)
  cat(paste(model_name, ": Loglikelihood vs mean body size plots saved.\n"))
  
}