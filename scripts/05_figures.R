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

read_edit <- function(filename) read_csv(filename, col_names = FALSE) %>% 
  mutate(file = filename)

# giometto <- 
#   list.files("giometto", full.names = TRUE) %>% 
#   map(read_edit) %>% 
#   bind_rows()

giometto <-
  read_csv("giometto_raw.csv")


giometto_scaled <- 
giometto %>% 
  mutate(sum_y = sum(X2), .by = file) %>% 
  mutate(scaled_y = X2/sum_y) %>% 
  mutate(mean_x = weighted.mean(X1, w = X2), .by = file) %>% 
  mutate(scaled_x = X1/mean_x) 

giometto_scaled %>% 
  ggplot(aes(scaled_x, scaled_y)) +
  geom_point(aes(col = file)) +
  scale_x_log10() +
  scale_y_log10()


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
      count(scaled_size, size) %>% 
      mutate(scaled_n = n/sum(n)) %>% 
      select(scaled_size, scaled_n, size)
  }
  
  # scaling the median parameter values
  plot_lines <- 
    scale_size_vec(rnorm(1e6, all_pars_median$mu,
                         all_pars_median$sigma)) %>% 
    mutate(dist = "normal") %>% 
    bind_rows(
      scale_size_vec(rlnorm(1e6, all_pars_median$meanlog, 
                            all_pars_median$sdlog)) %>% 
        mutate(dist = "lognormal") 
    ) %>% 
    rename(size_class = size) %>% 
    left_join(rls_bin_table, by = join_by(size_class)) %>% 
    mutate(scaled_n = scaled_n/(size_max-size_min))
  
  plot_dat <- 
    obsdata_rls_gridcell %>% 
    mutate(dat = "rls") %>% 
    bind_rows(obsdata_cbf_location %>% 
                mutate(dat = "cbf")) %>% 
    left_join(meansizes_gridcell, by = join_by(population, species)) %>% 
    left_join(plotdata_gridcell %>% 
                select(population, normal_better), 
              by = join_by(population)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    mutate(scaled_n = case_when(dat == "rls" ~ (n/population_n)/(size_max-size_min), 
                                TRUE ~ (n/population_n)))
  
  # scaled size vs scaled count (log-log, all populations)
  p_scaling_size_5 <-
    plot_dat %>%
    ggplot(aes(x = scaled_size, 
               y = scaled_n)) +
    geom_line(aes(y = scaled_n, group = population), 
              col = "grey70", 
              alpha = 0.1) +
    geom_line(aes(col = dist), 
              linewidth = 2,
              data = plot_lines) +
    geom_line(aes(x = scaled_x, 
                  y = scaled_y, 
                  col = file, 
                  group = file),
              linewidth = 0.5,
              lty = 2,
              data = giometto_scaled) +
    scale_x_log10(label = label_number(suffix = "x")) +
    scale_y_log10() +
    scale_color_manual(
      values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                 "lognormal" = rgb(29, 84, 128, maxColorValue=255),
                 "giometto/Dataset 0.csv" = "grey50"),
      label = c("normal" = "Normal fit",
                "lognormal" = "Lognormal fit",
                "giometto/Dataset 0.csv" = "13 protist species*")) +
    theme_cowplot(20) +
    theme(legend.position = c(0.1,0.1), 
          legend.justification = c(0,0),
          legend.title = element_blank()) +
    labs(
      x = "Relative body size (log)",
      y = "Relative N (log)"
    ) #+
    # guides(col = guide_legend(override.aes = list(lty = c(1,  2))))
  
  plot_dat %>% 
    select(scaled_size, scaled_n, population) %>% 
    bind_rows(plot_lines, )
  
  p_scaling_size_5_simple <-
    plot_lines %>%
      mutate(file = dist) %>% 
      bind_rows(giometto_scaled %>% 
                  mutate(dist = "giometto") %>% 
                  select(dist, file, scaled_size = scaled_x, scaled_n = scaled_y)) %>% 
    mutate(col = case_when(dist == "normal" ~ rgb(181, 144, 19, maxColorValue=255), 
                           dist == "lognormal" ~ rgb(29, 84, 128, maxColorValue=255), 
                           TRUE ~ "grey50"),
           lty = case_when(dist == "normal" ~ "solid", 
                           dist == "lognormal" ~ "solid", 
                           TRUE ~ "dashed"),
           lwd = case_when(dist == "normal" ~ 2, 
                           dist == "lognormal" ~ 2, 
                           TRUE ~ 0.5)) %>% 
    ggplot(aes(x = scaled_size, 
               y = scaled_n, 
               group = file)) +
    geom_line(aes(y = scaled_n, 
                  group = population),
              col = "grey70",
              alpha = 0.1, 
              data = plot_dat) +
    geom_line(aes(col = col, 
                  lty = lty, 
                  linewidth = lwd)) +
    scale_x_log10(label = label_number(suffix = "x")) +
    scale_y_log10() +
      scale_color_identity(labels = c("Lognormal fit", 
                                      "Normal fit", 
                                      "13 protist species*"), guide = "legend") +
      scale_linetype_identity() +
      scale_linewidth_identity() +
    theme_cowplot(20) +
    theme(legend.position = c(0.1,0.1), 
          legend.justification = c(0,0),
          legend.title = element_blank(), 
          plot.margin = margin(10, 20, 10, 10)) +
    labs(
      x = "Relative body size (log)",
      y = "Abundance density (log)"
    ) +
      guides(col = guide_legend(override.aes = list(lty = c(1, 1, 2), 
                                                    linewidth = c(2,2,0.5))))
  
    ggsave(filename = "output/figures/ms_figs/scaling_p5.png",
           plot = p_scaling_size_5_simple,
           height = 15,
           width = 15*1.618,
           units = "cm")
  
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
      # for(pop in c("gridcell", "ecoregion", "species")){
      
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
        
        
        main_plot <- function(.data) {
          .data %>% 
            mutate(level = "gridcell") %>% 
            bind_rows(plotdata_species%>% 
                        mutate(level = "species")) %>% 
            
            {if(inc_fished) . else filter(., !fished)} %>% 
            {if(inc_bimodal) . else filter(., !bimodal)} %>% 
            left_join(meansizes_gridcell %>% 
                        mutate(level = "gridcell") %>% 
                        bind_rows(meansizes_species%>% 
                                    mutate(level = "species")) %>% 
                        rename(dat = data), 
                      by = join_by(population, level, species, dat)) %>% 
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
                     fill = "grey95", 
                     col = "transparent") +
            annotate(geom = "rect",
                     xmin = min(plotdata_current$mean_size), 
                     xmax = max(plotdata_current$mean_size), 
                     ymin = q_cov(0.1),
                     ymax = q_cov(0.9), 
                     fill = "grey80", 
                     col = "transparent") +
            geom_point(alpha = 1, size = 2) +
            geom_point(col = "red", pch = 4, size = 2,
                       data = plotdata_current %>% filter(bimodal)) +
            scale_x_continuous(trans = "log10", 
                               labels = label_number(suffix="cm"), 
                               limits = range(plotdata_current$mean_size)) +
            scale_shape_manual(
              values = c("rls" = 21, "cbf" = 24), 
              labels = c("rls" = "Visual census",
                         "cbf" = "Exhaustive sampling")) +
            scale_color_manual(
              values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                         "lognormal" = rgb(29, 84, 128, maxColorValue=255)),
              labels = c("normal" = "Normal preferred",
                         "lognormal" = "Lognormal preferred")) +
            labs(x = "Mean body size (log)", 
                 y = "Coefficient of variation") +
          
            theme_cowplot(20) +
            theme(legend.position = "none", 
                  axis.title.y = element_blank(),
                  # panel.background = element_rect(colour = "black", size = 1)
                  ) 
        }
        
        side_plot <- function(.data){
          .data %>% 
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
        }
        
        pop_data <- 
          plotdata_gridcell %>% 
          {if(inc_fished) . else filter(., !fished)} %>% 
          {if(inc_bimodal) . else filter(., !bimodal)}
        
        spp_data <- 
          plotdata_species %>% 
          {if(inc_fished) . else filter(., !fished)} %>% 
          {if(inc_bimodal) . else filter(., !bimodal)}
        
        
        ppp <- 
        main_plot(pop_data) + theme(axis.title.x = element_blank(), 
                                    # axis.text.x = element_blank()
                                    ) +
          ylim(0, plot_ylim) +
          annotate("text", x = 50, 
                                y = plot_ylim*0.95,
                                label = "Population-level",
                                size = 10) +
          side_plot(pop_data) + xlim(0, plot_ylim) +
        main_plot(spp_data) + ylim(0, plot_ylim) +
          annotate("text", x = 50, 
                   y = plot_ylim*0.95,
                   label = "Species-level",
                   size = 10) +
          side_plot(spp_data) +xlim(0, plot_ylim) +
          plot_layout(design = {"
      AAAAAAB
      AAAAAAB
      CCCCCCD
      CCCCCCD"
          })
        
        leg <- 
        ggpubr::as_ggplot(
          get_legend(
            main_plot(pop_data) +
              theme(axis.title = element_blank(), 
                    axis.text = element_blank()) +
                  guides(colour = guide_legend(override.aes = list(size=6)),
                         pch = guide_legend(override.aes = list(size=6))) +
                  theme(legend.position = "bottom",
                        legend.justification = c(0.5,0.5),
                        legend.background = element_rect(color = "black"),
                        legend.margin=margin(10,10,10,10),
                        legend.title = element_blank()))) 
        
        cc <- 
          wrap_elements(ppp) +
          labs(tag = "Coeffient of variation") +
          theme(
            plot.tag = element_text(size = 20, angle = 90),
            plot.tag.position = "left"
          )
        
        p_simple <- cc + leg +
          plot_layout(design = {"
      AAAAAAA
      AAAAAAA
      AAAAAAA
      AAAAAAA
            #BBBBB#"
          })
        
        # 
        # 
        # plot_ylim <- max(layer_scales(
        #   main_plot(pop_data))$y$range$range,
        #   layer_scales(
        #     main_plot(pop_data))$y$range$range)
        # 
        # p_simple <-
        # plotdata_gridcell %>% 
        #   mutate(level = "gridcell") %>% 
        #   bind_rows(plotdata_species%>% 
        #               mutate(level = "species")) %>% 
        #   
        #   {if(inc_fished) . else filter(., !fished)} %>% 
        #   {if(inc_bimodal) . else filter(., !bimodal)} %>% 
        #   left_join(meansizes_gridcell %>% 
        #               mutate(level = "gridcell") %>% 
        #               bind_rows(meansizes_species%>% 
        #                           mutate(level = "species")) %>% 
        #               rename(dat = data), 
        #             by = join_by(population, level, species, dat)) %>% 
        #   ggplot() + 
        #   aes(
        #     x = mean_size, 
        #     y = cov_pref, 
        #     pch = dat,
        #     col = better_dist
        #   ) +
        #   annotate(geom = "rect",
        #            xmin = min(plotdata_current$mean_size), 
        #            xmax = max(plotdata_current$mean_size), 
        #            ymin = q_cov(0.975),
        #            ymax = q_cov(0.025), 
        #            fill = "grey95", 
        #            col = "transparent") +
        #   annotate(geom = "rect",
        #            xmin = min(plotdata_current$mean_size), 
        #            xmax = max(plotdata_current$mean_size), 
        #            ymin = q_cov(0.1),
        #            ymax = q_cov(0.9), 
        #            fill = "grey80", 
        #            col = "transparent") +
        #   geom_point(alpha = 1, size = 2) +
        #   geom_point(col = "red", pch = 4, size = 2,
        #              data = plotdata_current %>% filter(bimodal)) +
        #   scale_x_continuous(trans = "log10", 
        #                      labels = label_number(suffix="cm"), 
        #                      limits = range(plotdata_current$mean_size)) +
        #   scale_shape_manual(
        #     values = c("rls" = 21, "cbf" = 24), 
        #     labels = c("rls" = "Visual census",
        #                "cbf" = "Exhaustive sampling")) +
        #   scale_color_manual(
        #     values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
        #                "lognormal" = rgb(29, 84, 128, maxColorValue=255)),
        #     labels = c("normal" = "Normal preferred",
        #                "lognormal" = "Lognormal preferred")) +
        #   labs(x = "Mean body size (log)", 
        #        y = "Coefficient of variation") +
        #   facet_wrap(~level, labeller = as_labeller(c("gridcell" = "Population", 
        #                                                "species" = "Species")), 
        #              strip.position="top", 
        #              nrow = 2) +
        #   theme_cowplot(20) +
        #   theme(legend.position = "bottom", 
        #         panel.background = element_rect(colour = "black", size = 1), 
        #         strip.background = element_rect(fill="white", colour = "black", size = 1), 
        #         legend.justification = c(0.5,0.5),
        #         legend.background = element_rect(color = "black"),
        #         legend.title = element_blank()) +
        #   guides(colour = guide_legend(override.aes = list(size=5)),
        #          pch = guide_legend(override.aes = list(size=5))) 
        # 
        # 
        # plotdata_gridcell %>% 
        #     ggplot() +
        #     aes(
        #       x = cov_pref,
        #       fill = better_dist
        #     ) +
        #     scale_fill_manual(
        #       values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
        #                  "lognormal" = rgb(29, 84, 128, maxColorValue=255))) +
        #     geom_density(alpha = 0.3, col = "black") +
        #     coord_flip() +
        #     theme_void(20) +
        #     theme(legend.position = "none")
        #   
        # 
      #   p1 <- 
      #     plotdata_current %>% 
      #     ggplot() + 
      #     aes(
      #       x = mean_size, 
      #       y = cov_pref, 
      #       pch = dat,
      #       col = better_dist
      #     ) +
      #     annotate(geom = "rect",
      #              xmin = min(plotdata_current$mean_size), 
      #              xmax = max(plotdata_current$mean_size), 
      #              ymin = q_cov(0.975),
      #              ymax = q_cov(0.025), 
      #              fill = "grey95", 
      #              col = "transparent") +
      #     # annotate(geom = "rect",
      #     #          xmin = min(plotdata_current$mean_size), 
      #     #          xmax = max(plotdata_current$mean_size), 
      #     #          ymin = q_cov(0.05),
      #     #          ymax = q_cov(0.95), 
      #     #          fill = "grey70", 
      #     # col = "transparent") +
      #     annotate(geom = "rect",
      #              xmin = min(plotdata_current$mean_size), 
      #              xmax = max(plotdata_current$mean_size), 
      #              ymin = q_cov(0.1),
      #              ymax = q_cov(0.9), 
      #              fill = "grey80", 
      #              col = "transparent") +
      #     geom_point(alpha = 1) +
      #     geom_point(col = "red", pch = 4, 
      #                data = plotdata_current %>% filter(bimodal)) +
      #     # geom_point(col = "red", pch = 4, data = plotdata_current %>% filter(fished)) +
      #     annotate(geom = "text", 
      #              x = 90, 
      #              y = q_cov(0.9) - ((q_cov(0.9) - q_cov(0.1))/2),
      #              label = "80%", 
      #              size = 8) +
      #     # annotate(geom = "text", 
      #     #          x = 100, 
      #     #          y = q_cov(0.95) - ((q_cov(0.95) - q_cov(0.9))/2), 
      #     #          label = "90%", 
      #     #          size = 8) +
      #     annotate(geom = "text", 
      #              x = 90, 
      #              y = q_cov(0.975) - ((q_cov(0.975) - q_cov(0.9))/2), 
      #              label = "95%", 
      #              size = 8) +
      #     
      #     scale_x_continuous(trans = "log10", 
      #                        labels = label_number(suffix="cm"), 
      #                        limits = range(plotdata_current$mean_size)) +
      #     scale_shape_manual(
      #       values = c("rls" = 21, "cbf" = 24), 
      #       labels = c("rls" = "Visual census",
      #                  "cbf" = "Exhaustive sampling")) +
      #     scale_color_manual(
      #       values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
      #                  "lognormal" = rgb(29, 84, 128, maxColorValue=255)),
      #       labels = c("normal" = "Normal preferred",
      #                  "lognormal" = "Lognormal preferred")) +
      #     labs(x = "Log mean size", 
      #          y = "Coefficient of variation") +
      #     theme_cowplot(20) +
      #     theme(legend.position = "none")
      #   
      #   assign(x = paste("p1", 
      #                    pop, 
      #                    ifelse(inc_fished, "incfished", "unfished"), 
      #                    sep = "_"), 
      #          value = p1, 
      #          envir = .GlobalEnv)
      #   
      #   p2 <- 
      #     plotdata_current %>% 
      #     ggplot() +
      #     aes(
      #       x = cov_pref, 
      #       fill = better_dist
      #     ) +
      #     scale_fill_manual(
      #       values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
      #                  "lognormal" = rgb(29, 84, 128, maxColorValue=255))) +
      #     geom_density(alpha = 0.3, col = "black") +
      #     coord_flip() +
      #     theme_void(20) +
      #     theme(legend.position = "none") 
      #   
      #   
      #   assign(x = paste("p2",
      #                    pop, 
      #                    ifelse(inc_fished, "incfished", "unfished"), 
      #                    sep = "_"),
      #          value = p2, 
      #          envir = .GlobalEnv)
      #   
      #   ggsave(filename = fig_filename,
      #          plot = p1 + p2 + plot_layout(widths = c(6,1)),
      #          height = 15,
      #          width = 15*1.618,
      #          units = "cm")
      #   
      #   
      #   legend_only <- get_legend(
      #     get(paste0("p1_gridcell_", ifelse(inc_fished, "incfished", "unfished"))) + 
      #       guides(colour = guide_legend(override.aes = list(size=6)),
      #              pch = guide_legend(override.aes = list(size=6))) +
      #       theme(legend.position = "bottom",
      #             legend.justification = c(0.5,0.5),
      #             legend.background = element_rect(color = "black"),
      #             legend.margin=margin(10,10,10,10),
      #             legend.title = element_blank()))
      #   
        # plot_ylim <- max(layer_scales(
        #   get(paste0("p1_gridcell_",
        #              ifelse(inc_fished, "incfished", "unfished"))))$y$range$range,
        #   layer_scales(
        #     get(paste0("p1_species_",
        #                ifelse(inc_fished, "incfished", "unfished"))))$y$range$range)
      #   
      #   fig_design <- {"AAAAAABCCCCCCD
      # AAAAAABCCCCCCD
      # AAAAAABCCCCCCD
      # AAAAAABCCCCCCD
      # AAAAAABCCCCCCD
      # ####EEEEEE####"}
      #   
      #   fig_design <- {"AAAAAAB
      #     AAAAAAB
      # AAAAAAB
      # CCCCCCD
      # CCCCCCD
      #     CCCCCCD
      # EEEEEEE"}
      #   
      #   p2 <- 
      #     get(paste0("p2_gridcell_", 
      #                ifelse(inc_fished, "incfished", "unfished"))) + 
      #     xlim(0, plot_ylim) +
      #     plot_layout(tag_level = 'new')
      #   
      #   p4 <- 
      #     get(paste0("p2_species_", 
      #                ifelse(inc_fished, "incfished", "unfished"))) + 
      #     xlim(0, plot_ylim) +
      #     plot_layout(tag_level = 'new')
      #   
      #   p5 <- 
      #     ggpubr::as_ggplot(legend_only) +
      #     plot_layout(tag_level = 'new')
      #     
      #   p_cv_meansize <-
      #     get(paste0("p1_gridcell_", 
      #                ifelse(inc_fished, "incfished", "unfished"))) + 
      #     annotate("text", x = 10, 
      #              y = plot_ylim, 
      #              label = "Population-level", 
      #              size = 12) + 
      #     ylim(0, plot_ylim) +
      #     p2 +
      #     get(paste0("p1_species_", 
      #                ifelse(inc_fished, "incfished", "unfished"))) + 
      #     annotate("text", x = 10, 
      #              y = plot_ylim, 
      #              label = "Species-level", 
      #              size = 12) + 
      #     ylim(0, plot_ylim) + 
      #     theme(axis.title.y = element_blank()) +
      #     p4 + 
      #     p5 +
      #     plot_layout(design = fig_design) +
      #     plot_annotation(tag_levels = 'A')
        
        ggsave(filename = fig_filename,
               plot = p_simple,
               height = 20,
               width = 20*1.618,
               units = "cm")
        
      # } else {
      #   print(paste0(fig_filename, " figure already saved."))
      # }
    }
    
    # rm(inc_bimodal, p_cv_meansize, fig_design, plot_ylim, legend_only)
  }
  # rm(inc_fished)
}

# Estimated body size ==========================================================

for(inc_fished in c(FALSE, TRUE)){
  
  fig_filename <- paste0("output/figures/ms_figs/var_explained_", 
                         {if(inc_fished) "_incfished" else "_unfished"}, 
                         ".png")
  
  if(!file.exists(fig_filename) | force_run){
    
    if(!file.exists(paste0("output/tables/var_explained_bycv",  
                           ifelse(inc_fished, "_incfished", "_unfished"), 
                           ".csv")) | force_run){
      
      d1 <- 
        obsdata_rls_gridcell %>% 
        select(population, size_class, size_min, size_max, n, p_obs)%>% 
        mutate(dat = "rls")
      
      d2 <- 
        obsdata_cbf_location %>% 
        mutate(p = n/population_n) %>% 
        select(population, size_class, size_min, size_max, n, p_obs) %>% 
        mutate(dat = "cbf")
      
      joined_data <-  
        meansizes_gridcell %>% 
        right_join(bind_rows(d1, d2), 
                   by = join_by(population)) %>% 
        left_join(plotdata_gridcell) %>% 
        filter(!bimodal) %>% 
        {if(inc_fished) . else filter(., !fished)} %>% 
        select(population, species, size_class, mean_size, size_min, 
               size_max, n, p_obs, mu_50, sigma_50, meanlog_50, sdlog_50, better_dist, dat)
      
      size_classes <- 
        joined_data %>% 
        filter(dat == "rls") %>% 
        pull(size_class) %>% 
        unique() %>% 
        sort()
      
      out <- tibble()
      for(c in seq(0.1, 1.5, by = 0.01)){
        
        estimated_prob <- 
          joined_data %>% 
          mutate(mu = mean_size,
                 sd = mu*c,
                 sdlog = sqrt(log((c^2)+1)),
                 meanlog = log(mean_size) - ((sdlog^2)/2)) %>%
          mutate(p_norm = 
                   pnorm(size_max, mean = mu, sd = sd) -  
                   pnorm(size_min, mean = mu, sd = sd),
                 p_lnorm = 
                   plnorm(size_max, meanlog = meanlog, sdlog = sdlog) -  
                   plnorm(size_min, meanlog = meanlog, sdlog = sdlog), 
                 p_pref = ifelse(better_dist == "normal", p_norm, p_lnorm)) %>% 
          select(population, species,
                 size_class, better_dist, 
                 p_obs, p_norm, p_lnorm, p_pref, n, dat) %>% 
          filter(dat == "rls") %>% 
          mutate(pop_id = cur_group_id(), .by = population)
        
        
        out_norm <- c()
        out_lnorm <- c()
        out_pref <- c()
        for(ii in 1:max(estimated_prob$pop_id)){
          
          curr_dat <- 
            estimated_prob %>% 
            filter(pop_id == ii)
          
          out_norm[ii]  <- ks.test(curr_dat$p_obs, y = curr_dat$p_norm)$statistic %>% as.numeric()
          out_lnorm[ii] <- ks.test(curr_dat$p_obs, y = curr_dat$p_lnorm)$statistic %>% as.numeric()
          out_pref[ii]  <- ks.test(curr_dat$p_obs, y = curr_dat$p_pref)$statistic %>% as.numeric()
          
        }
        
        out <- 
          out %>% 
          bind_rows(
            tibble(norm_ks = unlist(out_norm), 
                   lnorm_ks = unlist(out_lnorm),
                   pref_ks = unlist(out_pref),
                   c = c))
        
        cat(c, "\n")
      }
      
      write_csv(out, paste0("output/tables/var_explained_bycv_ks",  ifelse(inc_fished, "_incfished", "_unfished"), ".csv"))
      rm(c, out, estimated_prob, lmer_norm, lmer_lnorm, lmer_pref, joined_data, d1, d2)
    }
    
    
    p <- 
      paste0("output/tables/var_explained_bycv_ks",  ifelse(inc_fished, "_incfished", "_unfished"), ".csv") %>% 
      read_csv(show_col_types = FALSE) %>% 
      pivot_longer(cols = contains("ks"), 
                   values_to = "ks",
                   names_to = "dist") %>% 
      summarise(mean_ks = mean(ks), 
                median_ks = median(ks), 
                lwr_ks = quantile(ks, 0.025), 
                upr_ks = quantile(ks, 0.875), 
                .by = c(c, dist)) %>% 
      ggplot(aes(x = c, 
                 y = mean_ks, 
                 color = dist)) +
      geom_line(linewidth = 2, alpha = 0.8) +
      scale_color_manual(values = c("norm_ks" = rgb(181, 144, 19, maxColorValue=255),
                                    "lnorm_ks" = rgb(29, 84, 128, maxColorValue=255), 
                                    "pref_ks" = "pink"),
                         labels = c("norm_ks" = "Normal",
                                    "lnorm_ks" = "Lognormal", 
                                    "pref_ks" = "Preferred")) +
      guides(color = guide_legend(override.aes = list(alpha = 1) ) ) +
      scale_y_continuous() +
      labs(x = "Assumed Coefficient of Variation", 
           y = "Kolmogorov-Smirnov statistic") +
      theme_cowplot(20) +
      theme(legend.position = c(0.95,0.05), 
            legend.justification = c(1,0),
            legend.title = element_blank(), 
            plot.background = element_rect(color = "transparent")) 
    
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
      # mutate(lmer_norm = predict(lmer_norm),
      #        lmer_lnorm = predict(lmer_lnorm)) %>%
      pivot_longer(cols = contains("norm"),
                   names_to = "dist",
                   values_to = "p_est") %>%
      ggplot(aes(x = p_est,
                 y = p_obs,
                 col = dist), pch = 21) +
      geom_point(alpha = 0.1) +
      geom_textabline(slope = 1, lty = 2, label = "Predicted = Observed", size = 6) +
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
    
    
    p2 <- p_obs_vs_exp + p + plot_annotation(tag_levels = 'A')
      # p_obs_vs_exp + inset_element(
      #   p,
      #   left =  0.6, right = 0.995,
      #   bottom =  0.05, top = 0.495)
    
    ggsave(filename = fig_filename,
           plot = p2,
           height = 15,
           width = 30,
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



# New figure 3


tab1 <- 
obsdata_rls_species %>% 
  mutate(population_indx = cur_group_id(), 
         .by = species) %>% 
  filter(population_indx < 7) %>% 
  left_join(meansizes_species) %>% 
  select(species, mean_size, size_class, size_min, size_max, p_obs)

tab1 %>% 
  select(species, mean_size) %>% 
  distinct()

tab1  %>% 
  mutate(mu = mean_size,
         sd = mu*0.35,
         sdlog = sqrt(log((0.35^2)+1)),
         meanlog = log(mean_size) - ((sdlog^2)/2)) %>% 
  mutate(p_norm = 
           pnorm(size_max, mean = mu, sd = sd) -  
           pnorm(size_min, mean = mu, sd = sd),
         p_lnorm = 
           plnorm(size_max, meanlog = meanlog, sdlog = sdlog) -  
           plnorm(size_min, meanlog = meanlog, sdlog = sdlog)) %>% 
  ggplot(aes(x = size_class,
             y = p_obs, 
             width = size_max-size_min)) +
  geom_rect(aes(xmin = size_min, xmax = size_max, ymin = 0, ymax = p_obs), alpha =.2, col = "grey50") + 
  facet_wrap(~species) +
  geom_point(aes(y = p_norm), col = rgb(181, 144, 19, maxColorValue=255)) +
  geom_path(aes(y = p_norm), col = rgb(181, 144, 19, maxColorValue=255)) +
  geom_point(aes(y = p_lnorm), col = rgb(29, 84, 128, maxColorValue=255)) +
  geom_path(aes(y = p_lnorm), col = rgb(29, 84, 128, maxColorValue=255)) +
  theme_cowplot(20) +
  labs(x = "Size class (cm)", 
       y = "Proportion in size class") +
  theme(legend.position="none",
        strip.background=element_rect(colour="black",
                                      fill="grey97"))
