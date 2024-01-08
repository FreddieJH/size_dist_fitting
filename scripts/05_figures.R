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

force_run <- FALSE

# Plotting data ================================================================


for(pop in c("species", 
             "ecoregion", 
             "gridcell")){
  
  cbf_pop <- ifelse(pop %in% c("gridcell", "ecoregion"), "location", pop)
  
  norm_data_rls <- 
    paste0("output/models/summary/pars/rls_", pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "rls")
  
  lnorm_data_rls <- 
    paste0("output/models/summary/pars/rls_", pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "rls")
  
  norm_data_cbf <- 
    paste0("output/models/summary/pars/cbf_", cbf_pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "cbf")
  
  lnorm_data_cbf <- 
    paste0("output/models/summary/pars/cbf_", cbf_pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "cbf")
  
  norm_data_rls_cv <- 
    paste0("output/models/summary/cv/rls_", pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_norm")) %>% 
    mutate(dat = "rls")
  
  lnorm_data_rls_cv <- 
    paste0("output/models/summary/cv/rls_", pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_lnorm")) %>% 
    mutate(dat = "rls")
  
  norm_data_cbf_cv <- 
    paste0("output/models/summary/cv/cbf_", cbf_pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_norm")) %>% 
    mutate(dat = "cbf")
  
  lnorm_data_cbf_cv <- 
    paste0("output/models/summary/cv/cbf_", cbf_pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    select(-logl_50) %>% 
    rename_with(.cols = -population, .fn = function(x) paste0(x, "_lnorm")) %>% 
    mutate(dat = "cbf")
  
  bimodal_pops_rls <- 
    paste0("input/data/data_cleaning/bimodal_rls_", pop, ".csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    pull(population)
  
  bimodal_pops_cbf <- 
    paste0("input/data/data_cleaning/bimodal_cbf_", cbf_pop, ".csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    pull(population)
  
  
  bind_rows(norm_data_rls, 
            norm_data_cbf) %>% 
    rename(logl_50_norm = logl_50) %>% 
    left_join(bind_rows(lnorm_data_rls, 
                        lnorm_data_cbf) %>% 
                rename(logl_50_lnorm = logl_50), 
              by = join_by(population, dat)) %>% 
    left_join(bind_rows(norm_data_rls_cv, 
                        norm_data_cbf_cv),
              by = join_by(population, dat)) %>% 
    left_join(bind_rows(lnorm_data_rls_cv, 
                        lnorm_data_cbf_cv),
              by = join_by(population, dat)) %>% 
    mutate(species = str_extract(population, "^.*(?=__)")) %>% 
    mutate(fished = species %in% fished_species) %>% 
    mutate(bimodal = population %in% c(bimodal_pops_rls, 
                                       bimodal_pops_cbf)) %>% 
    mutate(normal_better = case_when(
      is.na(mu_50) ~ FALSE,
      is.na(meanlog_50) ~ TRUE,
      logl_50_norm > logl_50_lnorm ~ TRUE, 
      logl_50_lnorm > logl_50_norm ~ FALSE
    )) %>% 
    mutate(better_dist = ifelse(normal_better, "normal", "lognormal")) %>% 
    mutate(cov_pref = ifelse(normal_better, cv_50_norm, cv_50_lnorm)) %>% 
    assign(x = paste0("plotdata_", pop),  
           value = ., 
           envir = .GlobalEnv)
  
}


# Scaling body size ============================================================

fig_filename <- "output/figures/ms_figs/scaling_distributions.png"

if(!file.exists(fig_filename)|force_run){
  
  # 10 randomly selected RLS populations at the gridcell level
  rls_sample_populations <- 
    obsdata_rls_gridcell %>% 
    mutate(population_indx = cur_group_id(), 
           .by = population) %>% 
    filter(population_indx < 10)
  
  # natural size vs natural count
  p_scaling_size_1 <- 
    rls_sample_populations %>% 
    left_join(meansizes_gridcell, by = join_by(population)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    mutate(scaled_n = n/population_n) %>%
    ggplot(aes(x = size_class, y = n, color = population)) +
    scale_color_simpsons() +
    scale_x_continuous(labels = label_number(suffix="cm")) +
    geom_line(aes(y = n)) +
    theme(legend.position = "none") +
    theme_cowplot(20) +
    labs(
      x = "Body size",
      y = "N"
    ) +
    theme(legend.position = "none")
  
  # natural size vs scaled count
  p_scaling_size_2 <- 
    rls_sample_populations %>% 
    left_join(meansizes_gridcell, by = join_by(population)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    mutate(scaled_n = n/population_n) %>%
    ggplot(aes(x = size_class, y = scaled_n, col = population)) +
    geom_line(aes(y = scaled_n)) +
    scale_color_simpsons() +
    scale_x_continuous(labels = label_number(suffix="cm")) +
    theme(legend.position = "none") +
    theme_cowplot(20) +
    labs(
      x = "Body size",
      y = "Relative N"
    ) +
    theme(legend.position = "none")
  
  # scaled size vs scaled count
  p_scaling_size_3 <- 
    rls_sample_populations %>% 
    left_join(meansizes_gridcell, by = join_by(population)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    mutate(scaled_n = n/population_n) %>%
    ggplot(aes(x = scaled_size, y = scaled_n, col = population)) +
    geom_line(aes(y = scaled_n)) +
    scale_color_simpsons() +
    scale_x_continuous(label = label_number(suffix = "x")) +
    theme(legend.position = "none") +
    theme_cowplot(20) +
    labs(
      x = "Relative body size",
      y = "Relative N"
    )+
    theme(legend.position = "none")
  
  # scaled size vs scaled count (log-log)
  p_scaling_size_4 <-
    rls_sample_populations %>% 
    left_join(meansizes_gridcell, by = join_by(population)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    mutate(scaled_n = n/population_n) %>%
    ggplot(aes(x = scaled_size, y = scaled_n, col = population)) +
    geom_line(aes(y = scaled_n)) +
    scale_x_log10(label = label_number(suffix = "x")) +
    scale_y_log10() +
    scale_color_simpsons() +
    theme(legend.position = "none") +
    theme_cowplot(20) +
    labs(
      x = "Relative body size (log)",
      y = "Relative N (log)"
    ) +
    theme(legend.position = "none")
  
  # median parameter estimates of normal and lognormal
  all_pars_median <- 
    plotdata_gridcell %>% 
    filter(!bimodal) %>% 
    drop_na() %>% 
    summarise(
      mu = median(mu_50),
      sigma = median(sigma_50),
      meanlog = median(meanlog_50),
      sdlog = median(sdlog_50)
    ) 
  
  rls_bin <- function(size) {
    
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
    
    rls_bin_table$size_class[.bincode(size, breaks = c(0, rls_bin_table$size_max))]
  }
  
  # function to scale body size and abundance
  scale_size_vec <- function(size_vector, lower_limit = 1.25) {
    tibble(size = size_vector) %>% 
      filter(size > lower_limit) %>% 
      mutate(size = rls_bin(size)) %>%
      mutate(mean_size = mean(size)) %>% 
      mutate(scaled_size = size/mean_size) %>% 
      count(scaled_size) %>% 
      mutate(scaled_n = n/sum(n)) %>% 
      select(scaled_size, scaled_n)
  }
  
  # scaling the median parameter values
  plot_lines <- 
    scale_size_vec(rnorm(1e6, all_pars_median$mu, all_pars_median$sigma)) %>% 
    mutate(dist = "normal") %>% 
    bind_rows(
      scale_size_vec(rlnorm(1e6, all_pars_median$meanlog, all_pars_median$sdlog)) %>% 
        mutate(dist = "lognormal") 
    )
  
  # scaled size vs scaled count (log-log, all populations)
  p_scaling_size_5 <- 
    obsdata_rls_gridcell %>% 
    bind_rows(obsdata_cbf_location) %>% 
    left_join(meansizes_gridcell, by = join_by(population)) %>% 
    left_join(plotdata_gridcell %>% 
                select(population, normal_better), 
              by = join_by(population)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    mutate(scaled_n = n/population_n) %>%
    ggplot(aes(x = scaled_size, 
               y = scaled_n)) +
    geom_line(aes(y = scaled_n, group = population), 
              col = "grey70", 
              alpha = 0.1) +
    geom_line(aes(col = dist), 
              linewidth = 2,
              data = plot_lines) +
    scale_x_log10(label = label_number(suffix = "x")) +
    scale_y_log10() +
    scale_color_manual(
      values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                 "lognormal" = rgb(29, 84, 128, maxColorValue=255))) +
    theme_cowplot(20) +
    theme(legend.position = c(0.1,0.1), 
          legend.justification = c(0,0),
          legend.title = element_blank()) +
    labs(
      x = "Relative body size (log)",
      y = "Relative abundance (log)"
    )
  
  # combining the five figures
  p_scaling_size_all <- 
    p_scaling_size_1 + 
    p_scaling_size_2 + 
    p_scaling_size_3 + 
    p_scaling_size_4 - 
    p_scaling_size_5 + 
    plot_layout(ncol=1)  + 
    plot_annotation(tag_levels = 'A')
  
  ggsave(filename = fig_filename,
         plot = p_scaling_size_all,
         height = 35,
         width = 25,
         units = "cm")
  
  rm(rls_sample_populations,
     plot_lines,
     scale_size_vec,
     all_pars_median,
     p_scaling_size_1, 
     p_scaling_size_2,
     p_scaling_size_3, 
     p_scaling_size_4, 
     p_scaling_size_5, 
     fig_filename)
  
} else {
  print(paste0(fig_filename, " figure already saved."))
  rm(fig_filename)
}

# Mean size vs CV ==============================================================

for(inc_bimodal in c(TRUE, FALSE)){
for(inc_fished in c(TRUE, FALSE)){
    for(pop in c("species", "ecoregion", "gridcell")){
      
      fig_filename <- paste0("output/figures/ms_figs/meansize_cv_", 
                             pop, 
                             {if(inc_fished) "_incfished" else "_unfished"}, 
                             {if(inc_bimodal) "_incbimodal" else "_nobimodal"},
                             ".png")
      
      if(!file.exists(fig_filename)|force_run){
        
        plotdata_current <- 
          get(paste0("plotdata_", pop)) %>% 
          {if(inc_fished) . else filter(., !fished)} %>% 
          {if(inc_bimodal) . else filter(., !bimodal)} %>% 
          left_join(get(paste0("meansizes_", pop)) %>% 
                      rename(dat = data), 
                    by = join_by(population, dat, species)) 
        
        
        q_cov <- function(quant) quantile(plotdata_current$cov_pref, quant)
        
        p1 <- 
          plotdata_current %>% 
          ggplot() + 
          aes(
            x = mean_size, 
            y = cov_pref, 
            pch = dat,
            col = better_dist
          ) +
          annotate(geom = "rect",
                   xmin = min(plotdata_current$mean_size), 
                   xmax = max(plotdata_current$mean_size), 
                   ymin = q_cov(0.975),
                   ymax = q_cov(0.025), 
                   fill = "grey90", 
                   col = "transparent") +
          # annotate(geom = "rect",
          #          xmin = min(plotdata_current$mean_size), 
          #          xmax = max(plotdata_current$mean_size), 
          #          ymin = q_cov(0.05),
          #          ymax = q_cov(0.95), 
          #          fill = "grey70", 
          # col = "transparent") +
          annotate(geom = "rect",
                   xmin = min(plotdata_current$mean_size), 
                   xmax = max(plotdata_current$mean_size), 
                   ymin = q_cov(0.1),
                   ymax = q_cov(0.9), 
                   fill = "grey50", 
                   col = "transparent") +
          geom_point(alpha = 1) +
          geom_point(col = "red", pch = 4, 
                     data = plotdata_current %>% filter(bimodal)) +
          # geom_point(col = "red", pch = 4, data = plotdata_current %>% filter(fished)) +
          annotate(geom = "text", 
                   x = 90, 
                   y = q_cov(0.9) - ((q_cov(0.9) - q_cov(0.1))/2),
                   label = "80%", 
                   size = 8) +
          # annotate(geom = "text", 
          #          x = 100, 
          #          y = q_cov(0.95) - ((q_cov(0.95) - q_cov(0.9))/2), 
          #          label = "90%", 
          #          size = 8) +
          annotate(geom = "text", 
                   x = 90, 
                   y = q_cov(0.975) - ((q_cov(0.975) - q_cov(0.9))/2), 
                   label = "95%", 
                   size = 8) +
          
          scale_x_continuous(trans = "log10", 
                             labels = label_number(suffix="cm"), 
                             limits = range(plotdata_current$mean_size)) +
          scale_shape_manual(
            values = c("rls" = 21, "cbf" = 24), 
            labels = c("rls" = "Reef Life Survey (binned)",
                       "cbf" = "Cryptobenthic (continuous)")) +
          scale_color_manual(
            values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                       "lognormal" = rgb(29, 84, 128, maxColorValue=255)),
            labels = c("normal" = "Normal preferred",
                       "lognormal" = "Lognormal preferred")) +
          labs(x = "Log mean size", 
               y = "Coefficient of variation") +
          theme_cowplot(20) +
          theme(legend.position = "none")
        
        assign(x = paste("p1", 
                         pop, 
                         ifelse(inc_fished, "incfished", "unfished"), 
                         sep = "_"), 
               value = p1, 
               envir = .GlobalEnv)
        
        p2 <- 
          plotdata_current %>% 
          ggplot() +
          aes(
            x = cov_pref, 
            fill = better_dist
          ) +
          scale_fill_manual(
            values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                       "lognormal" = rgb(29, 84, 128, maxColorValue=255))) +
          geom_density(alpha = 0.3, col = "black") +
          coord_flip() +
          theme_void(20) +
          theme(legend.position = "none") 
        
        
        assign(x = paste("p2",
                         pop, 
                         ifelse(inc_fished, "incfished", "unfished"), 
                         sep = "_"),
               value = p2, 
               envir = .GlobalEnv)
        
        ggsave(filename = fig_filename,
               plot = p1 + p2 + plot_layout(widths = c(6,1)),
               height = 15,
               width = 15*1.618,
               units = "cm")
        
        rm(p1, 
           p2,
           pop, 
           q_cov, 
           plotdata_current)
      
    

  
  legend_only <- get_legend(
    get(paste0("p1_gridcell_", ifelse(inc_fished, "incfished", "unfished"))) + 
      theme(legend.position = "bottom",
            legend.justification = c(0.5,0.5),
            legend.background = element_rect(color = "black"),
            legend.margin=margin(10,10,10,10),
            legend.title = element_blank()))
  
  plot_ylim <- max(layer_scales(
    get(paste0("p1_gridcell_", 
               ifelse(inc_fished, "incfished", "unfished"))))$y$range$range,
    layer_scales(
      get(paste0("p1_species_",
                 ifelse(inc_fished, "incfished", "unfished"))))$y$range$range)
  
  fig_design <- {"AAAAAABCCCCCCD
      AAAAAABCCCCCCD
      AAAAAABCCCCCCD
      AAAAAABCCCCCCD
      AAAAAABCCCCCCD
      ####EEEEEE####"}
  
  p_cv_meansize <-
    get(paste0("p1_gridcell_", 
               ifelse(inc_fished, "incfished", "unfished"))) + 
    annotate("text", x = 10, 
             y = plot_ylim, 
             label = "Population-level", 
             size = 12) + 
    ylim(0, plot_ylim) +
    get(paste0("p2_gridcell_", 
               ifelse(inc_fished, "incfished", "unfished"))) + 
    xlim(0, plot_ylim) +
    get(paste0("p1_species_", 
               ifelse(inc_fished, "incfished", "unfished"))) + 
    annotate("text", x = 10, 
             y = plot_ylim, 
             label = "Species-level", 
             size = 12) + 
    ylim(0, plot_ylim) + 
    theme(axis.title.y = element_blank()) +
    get(paste0("p2_species_", 
               ifelse(inc_fished, "incfished", "unfished"))) + 
    xlim(0, plot_ylim) + 
    ggpubr::as_ggplot(legend_only) +
    plot_layout(design = fig_design)
  
  ggsave(filename = fig_filename,
         plot = p_cv_meansize,
         height = 15,
         width = 40,
         units = "cm")
  
      } else {
        print(paste0(fig_filename, " figure already saved."))
      }
    }
  
  # rm(inc_bimodal, p_cv_meansize, fig_design, plot_ylim, legend_only)
}
# rm(inc_fished)
}

# Estimated body size ==========================================================

for(inc_fished in c(TRUE, FALSE)){
  
  fig_filename <- paste0("output/figures/ms_figs/var_explained_", 
                         {if(inc_fished) "_incfished" else "_unfished"}, 
                         ".png")
  
  if(!file.exists(fig_filename) | force_run){
    
    if(!file.exists(paste0("output/tables/var_explained_bycv",  
                           ifelse(inc_fished, "_incfished", "_unfished"), 
                           ".csv")) | force_run){
      
      d1 <- 
        obsdata_rls_gridcell %>% 
        select(population, size_class, size_min, size_max, n, p_obs)
      
      d2 <- 
        obsdata_cbf_location %>% 
        mutate(p = n/population_n) %>% 
        select(population, size_class, size_min, size_max, n, p_obs)
      
      joined_data <-  
        meansizes_gridcell %>% 
        right_join(bind_rows(d1, d2), 
                   by = join_by(population)) %>% 
        left_join(plotdata_gridcell) %>% 
        filter(!bimodal) %>% 
        {if(inc_fished) . else filter(., !fished)} %>% 
        select(population, species, size_class, mean_size, size_min, 
               size_max, n, p_obs, mu_50, sigma_50, meanlog_50, sdlog_50, better_dist)
      
      out <- tibble()
      for(c in seq(0, 1.5, by = 0.01)){
        
        estimated_prob <- 
          joined_data %>% 
          mutate(mu = mean_size,
                 sd = mu*c,
                 sdlog = sqrt(log((c^2)+1)),
                 meanlog = log(mean_size) - ((sdlog^2)/2)) %>%
          mutate(p_norm = pnorm(size_max, mean = mu, sd = sd) -  pnorm(size_min, mean = mu, sd = sd),
                 p_lnorm = plnorm(size_max, meanlog = meanlog, sdlog = sdlog) -  plnorm(size_min, meanlog = meanlog, sdlog = sdlog), 
                 p_pref = ifelse(better_dist == "normal", p_norm, p_lnorm)) %>% 
          select(population, species,
                 size_class, better_dist, 
                 p_obs, p_norm, p_lnorm, p_pref)
        
        lmer_norm <- lmerTest::lmer(p_obs ~ p_norm + (1|species), data = estimated_prob) 
        lmer_lnorm <-  lmerTest::lmer(p_obs ~ p_lnorm + (1|species), data = estimated_prob) 
        lmer_pref <-  lmerTest::lmer(p_obs ~ p_pref + (1|species), data = estimated_prob) 
        
        out <- 
          out %>% 
          bind_rows(
            tibble(cv = c, 
                   marginal_r2_normal = MuMIn::r.squaredGLMM(lmer_norm)[1], 
                   conditional_r2_normal = MuMIn::r.squaredGLMM(lmer_norm)[2], 
                   marginal_r2_lognormal = MuMIn::r.squaredGLMM(lmer_lnorm)[1], 
                   conditional_r2_lognormal = MuMIn::r.squaredGLMM(lmer_lnorm)[2], 
                   marginal_r2_pref = MuMIn::r.squaredGLMM(lmer_pref)[1], 
                   conditional_r2_pref = MuMIn::r.squaredGLMM(lmer_pref)[2])
          )
        
        cat(c, "\n")
      }
      
      write_csv(out, paste0("output/tables/var_explained_bycv",  ifelse(inc_fished, "_incfished", "_unfished"), ".csv"))
      rm(c, out, estimated_prob, lmer_norm, lmer_lnorm, lmer_pref, joined_data, d1, d2)
    }
    
    
    p <- 
      paste0("output/tables/var_explained_bycv",  ifelse(inc_fished, "_incfished", "_unfished"), ".csv") %>% 
      read_csv(show_col_types = FALSE) %>% 
      pivot_longer(cols = contains("r2"), values_to = "r2") %>% 
      mutate(dist = str_extract(name, "(?<=r2_)[a-z]+"),
             r2_type = str_extract(name, "[a-z]+(?=_)")) %>% 
      filter(r2_type == "marginal") %>% 
      ggplot(aes(x = cv, 
                 y = r2, 
                 color = dist)) +
      geom_line(linewidth = 2, alpha = 0.8) +
      scale_color_manual(values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                                    "lognormal" = rgb(29, 84, 128, maxColorValue=255), 
                                    "pref" = "pink"),
                         labels = c("normal" = "Normal",
                                    "lognormal" = "Lognormal", 
                                    "pref" = "Preferred")) +
      guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
      scale_y_continuous(label = label_percent(), limits = c(0,0.8)) +
      labs(x = "Assumed Coefficient of Variation", 
           y = "Variance explained") +
      theme_cowplot(15) +
      theme(legend.position = c(0.95,0.95), 
            legend.justification = c(1,1),
            legend.title = element_blank(), 
            plot.background = element_rect(color = "black")) 
    
    all_median_cov <- 
      plotdata_gridcell %>% 
      {if(inc_fished) . else filter(., !fished)} %>%
      pull(cov_pref) %>% 
      median() 
    
    estimated_prob <- 
      meansizes_gridcell %>%  
      mutate(mu = mean_size, 
             sd = mu*all_median_cov, 
             sdlog = sqrt(log((all_median_cov^2)+1)), 
             meanlog = log(mean_size) - ((sdlog^2)/2)) %>% 
      right_join(obsdata_rls_gridcell %>% {if(inc_fished) . else filter(., !fished)}) %>% 
      mutate(p_norm = pnorm(size_max, mean = mu, sd = sd) -  pnorm(size_min, mean = mu, sd = sd),
             plnorm_upper = plnorm(size_max, meanlog = meanlog, sdlog = sdlog), 
             plnorm_lower = plnorm(size_min, meanlog = meanlog, sdlog = sdlog),
             p_lnorm = plnorm(size_max, meanlog = meanlog, sdlog = sdlog) - plnorm(size_min, meanlog = meanlog, sdlog = sdlog)) %>% 
      select(population, species,
             size_class, 
             p_obs, p_norm, p_lnorm)
    
    lmer_norm <- lmerTest::lmer(p_obs ~ p_norm + (1|species), data = estimated_prob) 
    lmer_lnorm <-  lmerTest::lmer(p_obs ~ p_lnorm + (1|species), data = estimated_prob) 
    
    p_obs_vs_exp <-
      estimated_prob %>% 
      mutate(lmer_norm = predict(lmer_norm),
             lmer_lnorm = predict(lmer_lnorm)) %>%
      pivot_longer(cols = contains("norm"), 
                   names_to = "dist", 
                   values_to = "p_est") %>% 
      ggplot(aes(x = p_est, 
                 y = p_obs, 
                 col = dist), pch = 21) +
      geom_point(alpha = 0.1) +
      geom_abline(slope = 1, lty = 2) +
      labs(y = "Observed probability in size bin", 
           x = "Predicted probabiliy in size bin") +
      scale_x_continuous(label = label_percent()) +
      scale_y_continuous(label = label_percent()) +
      scale_color_manual(values = c("p_norm" = rgb(181, 144, 19, maxColorValue=255),
                                    "p_lnorm" = rgb(29, 84, 128, maxColorValue=255)), 
                         labels = c("p_norm" = "Normal",
                                    "p_lnorm" = "Lognormal")) +
      guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
      theme_cowplot(20) +
      theme(legend.position = c(0.05,0.95), 
            legend.justification = c(0,1),
            legend.title = element_blank())
    
    
    p2 <- 
      p_obs_vs_exp + inset_element(
        p, 
        left =  0.6, right = 0.995, 
        bottom =  0.05, top = 0.495)
    
    ggsave(filename = fig_filename,
           plot = p2,
           height = 25,
           width = 35,
           units = "cm")
    
    # ggsave(filename = "output/figures/ms_figs/fig3_main.png",
    #        plot = p_obs_vs_exp,
    #        height = 25,
    #        width = 35,
    #        units = "cm")
    # 
    # ggsave(filename = "output/figures/ms_figs/fig3_inset.png",
    #        plot = p,
    #        height = 25,
    #        width = 35,
    #        units = "cm")
  } else {
    print(paste0(fig_filename, " figure already saved."))
  }
}



