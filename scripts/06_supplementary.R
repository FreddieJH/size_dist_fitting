# About ========================================================================

# This script outputs all the figures for the main paper and the supplementary
# material.

# Packages =====================================================================

library(tidyverse)
library(arrow)
library(ggsci)
library(patchwork)
library(ggpubr)

# Parameters ===================================================================

force_run <- TRUE

# Intra-species CV variation ===================================================

output_filename <- 
  "output/figures/ms_figs/supplementary/intraspecies_cv_variation.png"

if(!file.exists(output_filename) | force_run){
  
  species_highpops <- 
    plotdata_gridcell %>% 
    filter(!bimodal) %>% 
    count(species, dat) %>% 
    arrange(desc(n)) %>% 
    head(20) %>% 
    pull(species) 
  
  
  p <-
    plotdata_gridcell %>% 
    filter(dat == "rls") %>% 
    filter(species %in% species_highpops) %>% 
    left_join(meansizes_species %>% 
                filter(data == "rls") %>% 
                mutate(species = str_extract(population, ".*(?=__)")) %>% 
                select(species, data, species_mean_size = mean_size), 
              by = join_by(species)) %>% 
    ggplot(aes(x = species_mean_size, 
               col = fct_reorder(species, species_mean_size),
               y = cov_pref)) +
    geom_point(pch = 21) +
    scale_x_continuous(label = label_number(suffix = "cm")) +
    theme_cowplot() +
    theme(legend.title = element_blank(), 
          legend.text = element_text(size = 12)) +
    labs(x = "Species mean size", 
         y = "Coefficient of Variation") 
  
  ggsave(filename = output_filename,
         plot = p,
         height = 15,
         width = 15*1.618,
         units = "cm")
  rm(p, species_highpops)
}

# Empirical distributions of populations with high CV ==========================

output_filename <- 
  "output/figures/ms_figs/supplementary/highcov_sdd.png"

if(!file.exists(output_filename) | force_run){
  
  current_plotdata <- 
    plotdata_gridcell %>% 
    filter(!bimodal)
  
  high_cov <- 
    current_plotdata %>% 
    filter(cov_pref > quantile(current_plotdata$cov_pref, 0.99)) %>% 
    pull(population)
  
  s_highcov_ssd <- 
    get(paste0("obsdata_rls_", "gridcell")) %>% 
    filter(population %in% high_cov) %>% 
    mutate(lat = str_extract(population, "(?<=__).\\d+") %>% as.numeric(), 
           lon = str_extract(population, "(?<=\\d_).\\d+") %>% as.numeric(),
           population_name = paste0(species, " \n (", lat, "°, ", lon, "°)")) %>% 
    ggplot() +
    aes(
      x = size_class, 
      y = n
    ) +
    geom_path() +
    geom_point() +
    facet_wrap(~population_name, scales = "free") +
    theme_cowplot(10) +
    labs(
      x = "Body size (cm)", 
      y = "Total abundance"
    )
  
  ggsave(filename = output_filename,
         plot = s_highcov_ssd,
         height = 20,
         width = 32,
         units = "cm")
  
  rm(high_cov, s_highcov_ssd)
  
}


# Distribution preference with mean size  ======================================

output_filename <- 
  "output/figures/ms_figs/supplementary/preference_meansize.png"

if(!file.exists(output_filename) | force_run){
  
  current_plotdata <- 
    plotdata_gridcell %>% 
    filter(!bimodal) %>% 
    left_join(meansizes_gridcell %>% 
                mutate(species = str_extract(population, ".*(?=__)")) %>% 
                select(population, data, pop_mean_size = mean_size), 
              by = join_by(population)) 
  
  logit_mod <-
    current_plotdata %>% 
    mutate(normal_better = as.numeric(normal_better)) %>% 
    glm(normal_better ~ pop_mean_size, data = ., family = "binomial")
  
  logit_plot <- 
    current_plotdata %>% 
    ggplot(aes(
      x = pop_mean_size,
      y = normal_better %>% as.numeric)) +
    geom_point(alpha = 0.6, pch = 21, size = 3) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"),
                col = "darkblue", fill = "darkblue") +
    # geom_line(aes(y = predict(logit_mod, type = "response"))) +
    scale_x_continuous(label = label_number(suffix = "cm")) +
    # scale_x_log10(label = label_number(suffix = "cm")) +
    theme_cowplot(20) +
    labs(
      x = "Mean size", 
      y = "Normal preferred"
    )
  
  ggsave(filename = output_filename,
         plot = logit_plot,
         height = 25,
         width = 35,
         units = "cm")
  
  rm(logit_plot, logit_mod, current_plotdata)
  
}


# CV vs. mean size model =======================================================

output_filename <- 
  "output/figures/ms_figs/supplementary/CV_meansize.png"

if(!file.exists(output_filename) | force_run){
  
  mod_data <- 
    plotdata_gridcell %>% 
    filter(!bimodal) %>% 
    left_join(meansizes_gridcell %>% 
                mutate(species = str_extract(population, ".*(?=__)")) %>% 
                select(population, data, mean_size), 
              by = join_by(population)) 
  
  cv_mod <- lmerTest::lmer(cov_pref ~ mean_size*data + (1|species), data = mod_data) 
  
  summary(cv_mod)
  m1_pred <- effects::effect(term = "mean_size", 
                             mod = cv_mod, 
                             xlevels = list(mean_size = seq(1, 100, 1))) %>% 
    as_tibble()
  
  p <- 
    mod_data %>% 
    ggplot() + 
    aes(x = mean_size) +
    geom_point(aes( y = cov_pref,
                    pch = data,
                    col = better_dist)) +
    geom_ribbon(aes(ymin = (lower), 
                    ymax = (upper)), 
                data = m1_pred, alpha = 0.5) +
    geom_line(aes(y = fit), 
              data = m1_pred) +
    scale_x_log10(label = label_number(suffix = "cm")) +
    scale_shape_manual(values = c("rls" = 21, "cbf" = 24), 
                       labels = c("rls" = "Reef Life Survey (binned)",
                                  "cbf" = "Cryptobenthic (continuous)")) +
    scale_color_manual(values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                                  "lognormal" = rgb(29, 84, 128, maxColorValue=255)),
                       labels = c("normal" = "Normal preferred",
                                  "lognormal" = "Lognormal preferred")) +
    labs(x = "Log mean size", 
         y = "Coefficient of variation") +
    theme_cowplot(20) +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(filename = output_filename,
         plot = p,
         height = 15,
         width = 32,
         units = "cm")
  rm(p)
}


# Population sizes =============================================================


get_npops <- Vectorize(function(src, pop, has_cv, rm_bimodal){
  
  temp_pop <- ifelse(pop == "location", "gridcell", pop)
  cbf_pop  <- ifelse(pop %in% c("ecoregion", "gridcell"), "location", pop)
  
  if(src == "all"){
    
    {if(!has_cv) bind_rows(get(paste("obsdata", "rls", pop, sep = "_")), 
                           get(paste("obsdata", "cbf", cbf_pop, sep = "_"))) else get(paste("plotdata", pop, sep = "_"))} %>% 
      {if(rm_bimodal) filter(., !bimodal) else .} %>% 
      pull(population) %>% 
      unique() %>% 
      length()
  } else {
    
    {if(!has_cv) get(paste("obsdata", src, pop, sep = "_")) else get(paste("plotdata", temp_pop, sep = "_")) %>% filter(dat == src)} %>% 
      {if(rm_bimodal) filter(., !bimodal) else .} %>% 
      pull(population) %>% 
      unique() %>% 
      length()
    
  }
})


if(!file.exists("output/tables/population_sizes.csv") | force_run){
  tibble(src = c("rls", "all")) %>%
    expand_grid(pop = c("species", "ecoregion", "gridcell")) %>% 
    bind_rows(tibble(src = "cbf", 
                     pop = c("species", "location"))) %>% 
    expand_grid(has_cv = c(TRUE, FALSE), 
                rm_bimodal = c(TRUE, FALSE)) %>% 
    filter(!(!has_cv & rm_bimodal)) %>%  
    mutate(popsize = get_npops(src, pop, has_cv, rm_bimodal)) %>% 
    write_csv("output/tables/population_sizes.csv")
}


# CV estimates =================================================================



get_cv <- Vectorize(function(src, pop, rm_bimodal, rm_fished, mod, quant){
  
  pop <- ifelse(pop == "location", "gridcell", pop)
  
  get(paste0("plotdata_", pop)) %>% 
    {if(src != "all") filter(., dat == src) else .} %>% 
    {if(rm_bimodal) filter(., !bimodal) else .} %>% 
    {if(rm_fished) filter(., !fished) else .} %>% 
    rename(cv_50_pref = cov_pref) %>% 
    pull(sym(paste0("cv_50_", mod))) %>% 
    quantile(quant, na.rm = TRUE)
  
})

if(!file.exists("output/tables/cv_quantiles.csv") | force_run){
  
  tibble(src = c("rls", "all")) %>%
    expand_grid(pop = c("species", "ecoregion", "gridcell")) %>% 
    bind_rows(tibble(src = "cbf", 
                     pop = c("species", "location"))) %>% 
    expand_grid(rm_bimodal = c(TRUE, FALSE), 
                rm_fished = c(TRUE, FALSE), 
                mod = c("norm", "lnorm", "pref"), 
                quant = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975)) %>% 
    mutate(cv = get_cv(src, pop, rm_bimodal, rm_fished, mod, quant)) %>% 
    pivot_wider(values_from = cv, 
                names_from = quant) %>% 
    write_csv("output/tables/cv_quantiles.csv")
}


read_csv("output/tables/cv_quantiles.csv") %>% 
  # pivot_longer(`0.025`:`0.975`) %>% 
  ggplot(aes(x = interaction(pop, src), y = `0.5`, col = mod, pch = rm_bimodal)) +
  geom_point()


normal_cov <- plotdata_gridcell %>% filter(better_dist == "normal") %>% pull(cov_pref)
lognormal_cov <- plotdata_gridcell %>% filter(better_dist == "lognormal") %>% pull(cov_pref)
normal_cov_spp <- plotdata_species %>% filter(better_dist == "normal") %>% pull(cov_pref)
lognormal_cov_spp <- plotdata_species %>% filter(better_dist == "lognormal") %>% pull(cov_pref)

# all populations
wilcox.test(
  normal_cov, 
  lognormal_cov
)

wilcox.test(normal_cov, conf.int = TRUE, conf.level = 0.95)
wilcox.test(lognormal_cov, conf.int = TRUE, conf.level = 0.95)

wilcox.test(normal_cov_spp, conf.int = TRUE, conf.level = 0.95)
wilcox.test(lognormal_cov_spp, conf.int = TRUE, conf.level = 0.95)

# unimodal populations
wilcox.test(
  get(paste0("plotdata_", "gridcell"))  %>% 
    filter(!bimodal,
           better_dist == "normal") %>% 
    pull(cov_pref), 
  get(paste0("plotdata_", "gridcell"))  %>% 
    filter(!bimodal,
           better_dist == "lognormal") %>% 
    pull(cov_pref))

# Number of transects vs CV estimate ===========================================


obsdata_rls_gridcell %>% 
  select(population, n_transects) %>% 
  distinct()

read_parquet("output/models/summary/cv/rls_gridcell_normal.parquet") %>% 
  left_join(obsdata_rls_gridcell %>% 
              select(population, n_transects) %>% 
              distinct(), 
            by = join_by(population)) %>% 
  ggplot(aes(cv_50, x = n_transects %>% log())) +
  geom_point() +
  stat_smooth()

# Proportion of normal vs lognormal


plotdata_gridcell %>% 
  filter(!bimodal) %>% 
  count(better_dist, 
        .by = c(dat))


# Sensitivity analysis (mincount) ==============================================

mincount_tbl <- 
  tibble(data = "rls", 
         mincount = c(20, 50, 200, 500)) %>% 
  bind_rows(tibble(
    data = "cbf", 
    mincount = c(2, 5, 20, 50)
  ))


# Data from https://www.fish.gov.au/reports/species
frdc <- read_csv("input/data/data_cleaning/frdc_fished.csv", 
                 show_col_types = FALSE) 

# From Audzijonyte et al (2020) Nat. Evo. Eco.
targeted <- 
  read_csv("input/data/data_cleaning/fishing_intensity.csv", 
           show_col_types = FALSE) %>% 
  filter(fishing_intensity > 1) # only showing targeted species

for(mod in c("lognormal", "normal")){
  
  for(i in 1:nrow(mincount_tbl)){
    
    out2 <- tibble()
    
    dat <- mincount_tbl$data[i]
    pop <- ifelse(dat == "rls", "gridcell", "location")
    mincount <- mincount_tbl$mincount[i]
    
    if(dat == "rls"){
      
      raw_data <- 
        read_parquet("input/data/raw/data_obs_cleaned.parquet") %>% 
        left_join(read_parquet("input/data/raw/survey_list_m1_aus.parquet"), 
                  by = join_by(survey_id)) %>% 
        rename(species = species_name) %>% 
        mutate(ecoregion = str_replace_all(ecoregion, "/", "-"), 
               ecoregion = str_replace_all(ecoregion, "Great Barrier Reef", "GBR"), 
               ecoregion = str_replace_all(ecoregion, "Central", "C"), 
               ecoregion = str_replace_all(ecoregion, "Northern", "N"), 
               ecoregion = str_replace_all(ecoregion, "Southern", "S"), 
               lat_grid = round(latitude), 
               lon_grid = round(longitude), 
               gridcell = paste(lat_grid, lon_grid, sep = "_"), 
               population = paste(species, !!sym(pop), sep = "__")) 
      
      n_transects <- 
        raw_data %>% 
        select(population, survey_date) %>% 
        distinct() %>% 
        count(population, name = "n_transects")
      
      simple_data <- 
        raw_data %>% 
        count(population, size_class, wt = n) 
      
      filtered_data <- 
        simple_data %>% 
        add_count(population, name = "n_sizebins") %>% 
        add_count(population, wt = n, name = "population_n") %>%
        filter(n_sizebins >= 4,
               population_n >= mincount) %>%
        arrange(desc(population_n)) 
      
      clean_data <-
        filtered_data %>% 
        left_join(rls_bin_table, by = join_by(size_class)) %>% 
        left_join(n_transects, by = join_by(population)) %>% 
        mutate(p_obs = n/population_n) %>% 
        mutate(species = str_extract(population, ".*(?=__)")) %>% 
        mutate(genus = str_extract(species, ".*(?=\\s)")) %>% 
        mutate(fished = (species %in% frdc$species)| 
                 (genus %in% frdc$genus)|
                 (species %in% targeted$species)) %>% 
        select(-genus)
      
      assign(paste0("obsdata_rls_", pop, "_", mincount), 
             value = clean_data, 
             envir = .GlobalEnv)
    } 
    if(dat == "cbf"){
      
      raw_data <- 
        read_csv("input/data/raw/cbf_raw.csv", 
                 show_col_types = FALSE) %>%
        mutate(sl_cm = SL/10,
               tl_cm = TL/10) %>% # now in cm, not mm
        select(location = Location,
               species = sciname,
               tl_cm
        ) 
      
      simple_data <- 
        raw_data %>%
        mutate(size_class = tl_cm %>% round(1), 
               size_min = size_class - 0.05, 
               size_max = size_class + 0.05) %>% 
        filter(str_detect(species, "^[A-Z]{1}[a-z]+\\s[a-z]+$")) %>% 
        mutate(population = paste(species, !!sym(pop), sep = "__")) %>% 
        count(population, size_class, size_min, size_max)
      
      filtered_data <- 
        simple_data %>% 
        add_count(population, wt = n,  name = "population_n") %>%
        filter(population_n >= mincount) %>% 
        arrange(desc(population_n)) 
      
      clean_data <-
        filtered_data %>% 
        mutate(p_obs = n/population_n) %>% 
        mutate(species = str_extract(population, ".*(?=__)"),
               genus = str_extract(species, ".*(?=\\s)")) %>% 
        mutate(fished = (species %in% frdc$species)| 
                 (genus %in% frdc$genus)|
                 (species %in% targeted$species)) %>% 
        select(-genus)
      
      assign(paste0("obsdata_cbf_", pop, "_", mincount), 
             value = clean_data, 
             envir = .GlobalEnv)
    }
    
    population_list <- 
      get(paste0("obsdata_", dat, "_", pop, "_", mincount)) %>% 
      pull(population) %>% 
      unique() 
    
    for(ii in population_list){
      
      current_data <-
        get(paste0("obsdata_", dat, "_", pop, "_", mincount)) %>%  
        filter(population == ii) %>% 
        arrange(size_class)
      
      p_firstbin <- 
        current_data %>% 
        pull(p_obs) %>% 
        .[1]
      
      if (p_firstbin < 0.5) {
        
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
        
        if(!is.null(summary(fit)$summary)){
          if(mod == "normal"){
            out <- 
              fit %>% 
              as_draws_df() %>% 
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
                        logl_50 = quantile(lp__, prob = 0.5)) %>% 
              # summarise_draws() %>% 
              mutate(population = ii, 
                     dat = dat, 
                     mod = mod)
          }
          if(mod == "lognormal"){
            out <- 
              fit %>% 
              as_draws_df() %>% 
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
                        logl_50 = quantile(lp__, prob = 0.5)) %>% 
              # summarise_draws() %>% 
              mutate(population = ii, 
                     dat = dat, 
                     mod = mod)
          }
          
          
          out2 <- 
            bind_rows(out2, out)
          
          write_csv(out2, 
                    paste0("output/tables/sensitivity_mincount_", 
                           dat, "_", 
                           mod, "_",
                           mincount, ".csv"))
        }
        gc()
      }
    }
  }
}

