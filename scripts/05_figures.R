
force_run <- FALSE

# Plotting data ================================================================

for(pop in c("species", 
             "ecoregion", 
             "gridcell")){
  
  cbf_pop <- ifelse(pop %in% c("gridcell", "ecoregion"), "location", pop)
  
  
  popsize_tbl <- 
    obsdata_rls_gridcell %>% 
    bind_rows(obsdata_rls_species) %>% 
    bind_rows(obsdata_rls_ecoregion) %>% 
    mutate(dat = "rls") %>% 
    bind_rows(obsdata_cbf_location %>% mutate(dat = "cbf")) %>% 
    bind_rows(obsdata_cbf_species%>% mutate(dat = "cbf")) %>% 
    select(population, dat, population_n) %>% 
    distinct()
  
  norm_data_rls <- 
    paste0("output/models/summary/pars/rls_", pop, "_normal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "rls")
  
  lnorm_data_rls <- 
    paste0("output/models/summary/pars/rls_", pop, "_lognormal.parquet") %>% 
    read_parquet() %>% 
    mutate(dat = "rls")
  
  norm_data_cbf <- 
    paste0("output/models/summary/pars/cbf_", 
           cbf_pop, "_normal.parquet") %>% 
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
    left_join(popsize_tbl, 
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
    filter(!species %in% c("Sepia apama", 
                           "Sepioteuthis australis", 
                           "Sepia plangon", 
                           "Aipysurus laevis", 
                           # "Bathytoshia brevicaudata", 
                           "Arctocephalus pusillus")) %>% # non-fish M1 species 
    assign(x = paste0("plotdata_", pop),  
           value = ., 
           envir = .GlobalEnv)
}


# Figure 1 ============================================================


# calculate CV for Giometto's protists:
# CV = sqrt(exp(sigma^2) - 1)

# for(n in 10^(1:8)){
#   lengths <- rnorm(n, mean = 10, sd = 2)
#   cv_lengths <- sd(lengths)/mean(lengths)
#   weights <- weight_from_length(lengths)
#   cv_weights <- sd(weights)/mean(weights)
#   
#   print(cv_weights/cv_lengths)
# }
# lengths <- rnorm(n, mean = 10, sd = 2)
# cv_lengths <- sd(lengths)/mean(lengths)
# weights <- weight_from_length(lengths)
# cv_weights <- sd(weights)/mean(weights)
# 
# cv_weights/cv_lengths
# 
# 
# cv_weights <- sqrt(exp(0.222)-1)
# 

# 
# weight_from_length  <- function(l) (pi/6)*l^3
# length_from_weight <- function(w)  pracma::nthroot(((6*w)/pi), n = 3)
# giometto_scaled <- 
#   read_csv("input/data/raw/giometto_raw.csv",
#            show_col_types = FALSE) %>% 
#   mutate(length_um = length_from_weight(X1)) %>% 
#   mutate(sum_y = sum(X2), .by = file) %>% 
#   mutate(scaled_y = X2/sum_y) %>% 
#   mutate(mean_x = weighted.mean(length_um, w = X2), .by = file) %>% 
#   mutate(scaled_x = length_um/mean_x)
# 
# giometto_scaled <- 
#   read_csv("input/data/raw/giometto_raw.csv",
#            show_col_types = FALSE) %>% 
#   mutate(length_um = length_from_weight(X1)) %>% 
#   mutate(sum_y = sum(X2), .by = file) %>% 
#   mutate(scaled_y = X2/sum_y) %>% 
#   mutate(mean_x = weighted.mean(length_um, w = X2), .by = file) %>% 
#   mutate(scaled_x = length_um/mean_x)
# 
# 
# 
# read_csv("input/data/raw/giometto_lengths.csv",
#          show_col_types = FALSE) %>% 
#   ggplot(aes(x = x, y = y, colour = as.factor(population))) +
#   geom_path() + 
#   scale_y_continuous(limits = c(0,1), 
#                      breaks = seq(0,1, by = 0.2)) +
#   theme(legend.position = "none")

giometto_scaled <- 
  read_csv("input/data/raw/giometto_lengths.csv",
           show_col_types = FALSE) %>% 
  filter(y > 0) %>% 
  mutate(sum_y = sum(y), 
         max_y = max(y),
         mean_x = weighted.mean(x, w = y),
         n = n(),
         .by = population) %>% 
  mutate(scaled_y = y/max_y, 
         scaled_x = x/mean_x) 


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
    mutate(size_class = rls_bin(size)) %>%
    mutate(mean_size = mean(size_class)) %>% 
    mutate(scaled_size = size_class/mean_size) %>% 
    count(scaled_size, size_class) %>% 
    left_join(rls_bin_table) %>% 
    mutate(norm_n = n/(size_max-size_min)) %>% 
    mutate(scaled_n = norm_n/max(norm_n)) %>% 
    select(scaled_size, scaled_n, size = size_class)
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
  left_join(rls_bin_table, by = join_by(size_class))

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
  mutate(normalised_n = case_when(
    dat == "rls" ~ (n/(size_max-size_min)), 
    TRUE ~ n)) %>% 
  mutate(sum_n = sum(n), 
         max_n = max(n),
         mean_x = weighted.mean(size_class, w = n),
         .by = population) %>% 
  mutate(scaled_n = n/max_n, 
         scaled_size = size_class/mean_x) 



# scaled size vs scaled count (log-log, all populations)

pdat <- 
  plot_lines %>%
  mutate(file = dist) %>% 
  bind_rows(giometto_scaled %>% 
              mutate(dist = "giometto", file = as.factor(population)) %>% 
              select(dist, file, scaled_size = scaled_x, 
                     scaled_n = scaled_y)) %>% 
  mutate(col = case_when(dist == "normal" ~ rgb(181, 144, 19, 
                                                maxColorValue=255), 
                         dist == "lognormal" ~ rgb(29, 84, 128, 
                                                   maxColorValue=255), 
                         TRUE ~ "black"),
         lty = case_when(dist == "normal" ~ "solid", 
                         dist == "lognormal" ~ "solid", 
                         TRUE ~ "dashed"),
         lwd = case_when(dist == "normal" ~ 2, 
                         dist == "lognormal" ~ 2, 
                         TRUE ~ 0.5)) 
fig1a <-
  pdat %>% 
  bind_rows(tibble(col = "dummy_var")) %>% 
  ggplot(aes(x = scaled_size, 
             y = scaled_n, 
             group = file)) +
  geom_line(aes(y = scaled_n, 
                group = population, 
                colour = col),
            alpha = 0.1, 
            data = plot_dat %>% 
              mutate(col = case_when(
                normal_better ~ rgb(181, 144, 19, 
                                    maxColorValue=255), 
                !normal_better ~ rgb(29, 84, 128, 
                                     maxColorValue=255), 
                TRUE ~ "grey80"),
                lty = case_when(normal_better ~ "solid", 
                                !normal_better ~ "solid", 
                                TRUE ~ "dashed"),
                lwd = case_when(normal_better ~ 2, 
                                !normal_better ~ 2, 
                                TRUE ~ 0.5)) ) +
  geom_line(col = "black", lwd = 3, data = pdat %>% filter(file %in% c("normal", "lognormal"))) +
  geom_line(aes(col = col,
                lty = lty,
                linewidth = lwd)) +
  scale_x_continuous(label = label_number(suffix = "x")) +
  # scale_y_log10() +
  scale_color_identity(labels = c(
                                "Observed (Normal preferred)", 
                                "Observed (Lognormal preferred)",
                                "Normal fit",
                                "Lognormal fit", 
                                "13 protist species"
                                ), 
                     guide = "legend") +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  theme_cowplot(20) +
  theme(legend.position = c(1,1), 
        legend.justification = c(1,1),
        legend.title = element_blank(), 
        plot.margin = margin(10, 20, 10, 10)) +
  labs(
    x = "Body size relative to mean size",
    y = "Abundance relative to max"
  ) +
  guides(col = guide_legend(override.aes = list(lty = c(1, 1, 1, 1, 2), 
                                                linewidth = c(0.3, 0.3, 3, 3, 0.5), 
                                                colour = c(rgb(181, 144, 19, maxColorValue=255),
                                                           rgb(29, 84, 128, maxColorValue=255),
                                                           rgb(181, 144, 19, maxColorValue=255),
                                                           rgb(29, 84, 128, maxColorValue=255),
                                                           "black"), 
                                                alpha = 1), 
                            label.position = "left"))


out <- c()
for(i in 1:10) {
  set.seed(1)
  out <- c(out, sample(((320*(i-1))+1):((320*(i-1)+1)+320), 1))
}

range_bodysizes <- 
  plot_dat %>% 
  filter(population %in% plotdata_gridcell$population) %>% 
  select(population, mean_size) %>% 
  distinct() %>% 
  arrange(mean_size) %>% 
  rownames_to_column("size_num") %>% 
  filter(size_num %in% out)

set.seed(1)
cbf_samples <- 
  plot_dat %>% 
  filter(population_n > 100) %>% 
  filter(dat=="cbf") %>% 
  select(population, mean_size) %>% 
  distinct() %>% 
  arrange(mean_size) %>% 
  rownames_to_column("size_num") %>% 
  filter(size_num %in% sample(1:nrow(.), 2))

extra <- 
  tibble(species = c(range_bodysizes$population, 
                     cbf_samples$population) %>% 
           str_extract(pattern = ".*(?=__)")) %>% 
  filter(species %in% 
           {read_csv("input/data/processed/spp_shapes.csv", 
                     show_col_types = FALSE) %>% 
               pull(species)})

fig1b <-
  plot_dat %>% 
  filter(population %in% c(range_bodysizes$population, cbf_samples$population)) %>%
  mutate(highlight_spp = species %in% extra$species) %>% 
  ggplot() +
  geom_line(aes(x = scaled_size, 
                y = scaled_n,
                col = dat,
                group = population),
            alpha = 0.5, 
            linewidth = 1) +
  geom_line(aes(x = scaled_size,
                y = scaled_n,
                col = dat,
                group = population),
            data = . %>% filter(highlight_spp),
            linewidth = 1.5) +
  scale_x_continuous(label = label_number(suffix = "x")) +
  scale_color_manual(values = c("rls" = "#56B4E9",
                                "cbf" = "#E69F00"),
                     labels = c("rls" = "In-situ observation",
                                "cbf" = "Specimen collection")) +
  theme_cowplot(20) +
  theme(legend.position = c(1,1), 
        legend.justification = c(1,1),
        legend.title = element_blank()) +
  labs(
    x = "Body size relative to mean size",
    y = "Abundance relative to max"
  ) +
  guides(col = guide_legend(label.position = "left"))

fig1 <- 
  fig1a + 
  fig1b + 
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A')

ggsave(filename = "output/figures/ms_figs/scaling_distributions.png",
       plot = fig1,
       height = 25,
       width = 30,
       units = "cm")

# Figure 2 =====================================================================

for(inc_bimodal in c(TRUE, FALSE)){
  for(inc_fished in c(TRUE, FALSE)){
    
    fig_filename <- paste0("output/figures/ms_figs/meansize_cv_", 
                           {if(inc_fished) "incfished" else "unfished"}, 
                           {if(inc_bimodal) "_incbimodal" else "_nobimodal"},
                           ".png")
    
    
    if(!file.exists(fig_filename)|force_run){
      
      # q_cov <- function(quant) quantile(.data$cov_pref, quant)
      
      main_plot <- function(.data) {
        .data %>% 
          ggplot() + 
          aes(
            x = mean_size, 
            y = cov_pref, 
            pch = col
          ) +
          geom_hline(yintercept = 0.3, lty = 2, col = "black", alpha = 0.5) +
          geom_hline(yintercept = 0.34, lty = 1, col = "#56B4E9", alpha = 0.5) +
          geom_hline(yintercept = 0.37, lty = 1, col = "#CC79A7", alpha = 0.5) +
          geom_hline(yintercept = 0.23, lty = 1, col = "#E69F00", alpha = 0.5) +
          geom_hline(yintercept = 0.33, lty = 1, col = "#000000", alpha = 0.5) +
          geom_point(alpha = 1, 
                     aes(col = col, 
                         size = (population_n %>% log10()))) +
          geom_point(col = "red",
                     pch = 4, 
                     size = 2,
                     data = .data %>% filter(bimodal)) +
          scale_size_continuous(range = c(0.1,5),
                                breaks = c(2,3,5),
                                labels =  c("100", "1K", "100K"),
                                transform = "log10", 
                                name = "Sample size") + 
          scale_x_continuous(transform = "log10", 
                             labels = label_number(suffix="cm"), 
                             limits = range(.data$mean_size)) +
          scale_shape_manual(
            values = c(
              "rls_norm" = 21, 
              "rls_lnorm" = 21, 
              "cbf_norm" = 24, 
              "cbf_lnorm" = 24),
            label = c(
              "rls_norm" = "Visual census (normal preferred)", 
              "rls_lnorm" = "Visual census (lognormal preferred)", 
              "cbf_norm" = "Cryptobenthic (normal preferred)", 
              "cbf_lnorm" = "Cryptobenthic (lognormal preferred)")) +
          scale_color_manual(values = c(
            "rls_norm" = "#56B4E9", 
            "rls_lnorm" = "#CC79A7", 
            "cbf_norm" = "#E69F00", 
            "cbf_lnorm" = "#000000"),
            label = c(
              "rls_norm" = "Visual census (normal preferred)", 
              "rls_lnorm" = "Visual census (lognormal preferred)", 
              "cbf_norm" = "Cryptobenthic (normal preferred)", 
              "cbf_lnorm" = "Cryptobenthic (lognormal preferred)")) +
          labs(x = "Mean body length (log)", 
               y = "Coefficient of variation in body size") +
          
          theme_cowplot(20) +
          theme(legend.position = "none", 
                axis.title.y = element_blank(),
                # panel.background = element_rect(colour = "black", size = 1)
          ) 
      }
      
      side_plot <- function(.data, ylim){
        .data %>% 
          ggplot() +
          aes(
            x = cov_pref,
            # fill = better_dist, 
            fill = col,
          ) +
          scale_fill_manual(values = c(
            "rls_norm" = "#56B4E9", 
            "rls_lnorm" = "#CC79A7", 
            "cbf_norm" = "#E69F00", 
            "cbf_lnorm" = "#000000"),
            label = c(
              "rls_norm" = "Visual census (normal preferred)", 
              "rls_lnorm" = "Visual census (lognormal preferred)", 
              "cbf_norm" = "Cryptobenthic (normal preferred)", 
              "cbf_lnorm" = "Cryptobenthic (lognormal preferred)"))+
          geom_density(alpha = 0.3, col = "black") +
          scale_x_continuous(position = "top", 
                             limits = c(0, ylim)) +
          coord_flip() +
          theme(legend.position = "none", 
                plot.background = element_rect(fill = "white"),
                panel.background = element_rect(fill = "white"),
                axis.line.x = element_blank(),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.y = element_line(), 
                axis.title.y = element_blank(), 
                axis.text.y = element_text(margin = margin(t = 0, 
                                                           r = 0, 
                                                           b = 0, 
                                                           l = 10), 
                                           size = 15), 
                axis.ticks = element_line(size = 1, colour = "black"))
      }
      
      pop_data <- 
        plotdata_gridcell %>% 
        left_join(meansizes_gridcell) %>%
        {if(inc_fished) . else filter(., !fished)} %>% 
        {if(inc_bimodal) . else filter(., !bimodal)} %>% 
        mutate(col = case_when(
          dat=="rls" & normal_better ~ "rls_norm",
          dat=="rls" & !normal_better ~ "rls_lnorm",
          dat=="cbf" & normal_better ~ "cbf_norm",
          dat=="cbf" & !normal_better ~ "cbf_lnorm",
        )) 
      
      
      spp_data <- 
        plotdata_species %>% 
        mutate(data = dat) %>%
        left_join(meansizes_species) %>% 
        {if(inc_fished) . else filter(., !fished)} %>% 
        {if(inc_bimodal) . else filter(., !bimodal)}%>% 
        mutate(col = case_when(
          dat=="rls" & normal_better ~ "rls_norm",
          dat=="rls" & !normal_better ~ "rls_lnorm",
          dat=="cbf" & normal_better ~ "cbf_norm",
          dat=="cbf" & !normal_better ~ "cbf_lnorm",
        ))
      
      plot_ylim <- 
        max(layer_scales(
          main_plot(pop_data))$y$range$range,
          layer_scales(
            main_plot(spp_data))$y$range$range)
      
      
      ppp <- 
        main_plot(pop_data) + 
        theme(axis.title.x = element_blank(), 
              # axis.text.x = element_blank()
        ) +
        ylim(0, plot_ylim) +
        annotate("text",
                 x = 1, 
                 hjust = 0, 
                 y = plot_ylim*0.95,
                 label = "Population-level",
                 size = 10) +
        side_plot(pop_data %>% 
                    filter(col != "cbf_lnorm"), 
                  plot_ylim) + 
        main_plot(spp_data) + ylim(0, plot_ylim) +
        annotate("text", 
                 x = 1, 
                 hjust = 0, 
                 y = plot_ylim*0.95,
                 label = "Species-level",
                 size = 10) +
        side_plot(spp_data %>% filter(col != "cbf_lnorm"), 
                  plot_ylim)  +
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
              theme(axis.title = element_text(), 
                    axis.text = element_text()) +
              guides(colour = guide_legend(override.aes = list(size=5, stroke = 2),
                                           nrow = 2,
                                           byrow = TRUE,
                                           title = element_blank()), 
                     pch = guide_legend(nrow = 2,
                                        byrow = TRUE,
                                        title = element_blank())) +
              theme(legend.position = "bottom", 
                    legend.box = "vertical",
                    legend.justification = c(0.5,0.5),
                    legend.background = element_rect(color = "black"),
                    legend.margin=margin(5,5,5,6),
                    legend.title = element_text())
          )) 
      
      cc <- 
        wrap_elements(ppp) +
        labs(tag = "Coeffient of variation in body size") +
        theme(
          plot.tag = element_text(size = 20, angle = 90),
          plot.tag.position = "left"
        )
      
      p_simple <- 
        cc + 
        leg +
        plot_layout(design = {"
      AAAAAAA
      AAAAAAA
      AAAAAAA
      AAAAAAA
            #BBBBB#"
        })
      # p_simple
      
      ggsave(filename = fig_filename,
             plot = p_simple,
             height = 20,
             width = 20*1.618,
             units = "cm")
      
    }
  }
}

# Figure 3 =====================================================================

fig3_filename <- "output/figures/ms_figs/lbspr_comparison.png"

if(!file.exists(fig3_filename) | force_run){
  
  all_obs_species <- 
    obsdata_rls_species %>%
    select(population, population_n) %>%
    mutate(data = "rls") %>%
    bind_rows(obsdata_cbf_species %>%
                select(population, population_n)%>%
                mutate(data = "cbf")) %>%
    left_join(meansizes_species) %>%
    distinct() %>%
    select(species, mean_size, data, population_n) %>%
    arrange(mean_size)
  
  # f = definitely fished, l = likely fished in aus, u = unlikely fished
  spp_fishing_lvl <- 
    read_csv("input/data/raw/prince_spp_rss.csv") %>% 
    select(species, fished) %>% 
    distinct()
  
  
  # cleaning the prince data (select the closest std_mk to the median mk of the species)
  prince_mean_values <-
    read_csv("input/data/raw/prince_2023_data.csv", 
             show_col_types = FALSE) %>% 
    janitor::clean_names() %>% 
    mutate(std_linf = published_linf_mm/10) %>% 
    select(
      species = species_name, 
      std_linf, 
      std_k = standardised_k, 
      std_mk = standardised_m_k,
      hoenig_m, 
      largest_fish_in_sample_mm
    ) %>% 
    drop_na(std_mk, std_linf) %>% 
    # mutate(species_median_mk = median(std_mk, na.rm = TRUE),
    #        .by = species) %>% 
    # mutate(abs_diff_mk = abs(median(std_mk, na.rm = TRUE)-std_mk)) %>% 
    filter(abs(median(std_mk, na.rm = TRUE)-std_mk) == 
             min(abs(median(std_mk, na.rm = TRUE)-std_mk), 
                 na.rm = TRUE), 
           .by = species) %>% 
    # filter(min_diff == abs_diff_mk) %>% 
    filter(abs((largest_fish_in_sample_mm/10)-std_linf) == 
             min(abs((largest_fish_in_sample_mm/10)-std_linf), 
                 na.rm = TRUE), 
           .by = species) 
  
  # species from prince that have 
  prince_species <- 
    prince_mean_values %>% 
    pull(species) %>% 
    unique()
  
  # all species that appear in prince and obs data
  species_in_common_all <- 
    all_obs_species %>% 
    filter(species %in% prince_species) %>% 
    left_join(spp_fishing_lvl, by = join_by(species)) %>% 
    filter(fished != "f")
  
  run_lbspr <- function(Linf, MK, BinWidth, BinMin, BinMax){
    ## Creating empty LB_pars object
    MyPars <- new("LB_pars")
    ## Fill in the required parameters in LB_pars
    #Fill in with the correct values
    MyPars@Linf <- Linf
    MyPars@MK <- MK
    MyPars@BinWidth <- BinWidth
    MyPars@BinMin <- BinMin #Lower bound of minimum size bin
    MyPars@BinMax <- BinMax #Upper bound of maximum size bin
    #Just fill in with arbitrary value.
    MyPars@L50 <- 0.5*MyPars@Linf #0 < L50 < Linf
    MyPars@L95 <- 0.6*MyPars@Linf #L50 < L95 < Linf
    MyPars@SL50 <- 0.5*MyPars@Linf #0 < L50 < Linf
    MyPars@SL95 <- 0.6*MyPars@Linf #L50 < L95 < Linf
    MyPars@FM <- 0 #F/M >= 0
    ## Run the simulation to construct the size composition and many other outputs
    MySim <- LBSPRsim(MyPars)
    ## Extract the size composition in unfished population
    # L <- MySim@pLPop[,1] #Mid value of each size bin
    # pL <- MySim@pLPop[,2] #Proportion of fish at size bin
    return(list(L = MySim@pLPop[,1], pL = MySim@pLPop[,2]))
  }
  # function(Linf, MK, Lc, alp, nl){
  output_lbspr <- 
    read_csv("input/data/processed/rls_obsdata_species_4_200.csv") %>% 
    filter(species %in% species_in_common_all$species) %>% 
    nest(.by = species, .key = "data_list") %>% 
    left_join(species_in_common_all) %>% 
    left_join(prince_mean_values) %>% 
    mutate(lc = 0.1*std_linf, 
           alp = 0.7) %>%
    drop_na() %>% 
    filter(data == "rls") %>% 
    mutate(output = 
             pmap(
               .l = list(
                 Linf = std_linf, 
                 MK = std_mk, 
                 BinWidth = 2, 
                 BinMin = 0, 
                 BinMax = 400), 
               .f = run_lbspr)) %>% 
    mutate(L = map(output, function(x) x$L), 
           pL = map(output, function(x) x$pL)) %>% 
    select(-c(data_list, output)) %>% 
    unnest(cols = c(L, pL)) %>% 
    select(species, L, lbspr = pL) %>% 
    mutate(size_class =  rls_bin_table$size_class[.bincode(L, breaks = rls_bin_table$size_min)]) %>% 
    summarise(n = sum(lbspr), .by = c(size_class, species)) %>% 
    drop_na() %>% 
    mutate(total_n = sum(n), .by = c(species)) %>% 
    mutate(p_lbspr = n/total_n) %>% 
    select(species, size_class, p_lbspr) %>% 
    left_join(rls_bin_table) 
  
  
  
  
  
  output_lbspr_setmk <- 
    read_csv("input/data/processed/rls_obsdata_species_4_200.csv") %>% 
    filter(species %in% species_in_common_all$species) %>% 
    nest(.by = species, .key = "data_list") %>% 
    left_join(species_in_common_all) %>% 
    left_join(prince_mean_values) %>% 
    mutate(lc = 0.1*std_linf, 
           alp = 0.7) %>%
    drop_na() %>% 
    filter(data == "rls") %>% 
    mutate(output = 
             pmap(
               .l = list(
                 Linf = std_linf, 
                 MK = 1.5, 
                 BinWidth = 2, 
                 BinMin = 0, 
                 BinMax = 400), 
               .f = run_lbspr)) %>% 
    mutate(L = map(output, function(x) x$L), 
           pL = map(output, function(x) x$pL)) %>% 
    select(-c(data_list, output)) %>% 
    unnest(cols = c(L, pL)) %>% 
    select(species, L, lbspr = pL) %>% 
    mutate(size_class =  rls_bin_table$size_class[.bincode(L, breaks = rls_bin_table$size_min)]) %>% 
    summarise(n = sum(lbspr), .by = c(size_class, species)) %>% 
    drop_na() %>% 
    mutate(total_n = sum(n), .by = c(species)) %>% 
    mutate(p_lbspr_mk15 = n/total_n) %>% 
    select(species, size_class, p_lbspr_mk15) %>% 
    left_join(rls_bin_table) 
  
  
  output_data <- 
    output_lbspr %>% 
    left_join(output_lbspr_setmk) %>% 
    left_join(read_csv("input/data/processed/rls_obsdata_species_4_200.csv") %>% 
                select(species, p_obs, size_min, size_max)) %>% left_join(prince_mean_values) %>% 
    left_join(species_in_common_all %>% filter(data == "rls") %>% select(species, mean_size) %>% distinct()) %>% 
    mutate(p_obs = replace_na(p_obs, 0))
  
  
  set.seed(1)
  sample_output <-
    output_data %>% 
    arrange(std_mk) %>% 
    count(species) %>% 
    select(-n) %>% 
    slice_sample(n = 7) %>% 
    left_join(output_data, by = "species") 
  
  mk_text <- 
    sample_output %>% 
    filter(p_obs != 0) %>% 
    summarise(max_x = max(size_class, na.rm = TRUE), 
              max_y = max(c(p_obs, p_lbspr), na.rm = TRUE), 
              std_mk = round(unique(std_mk), 2), 
              std_linf = round(unique(std_linf),0), 
              .by = species) 
  
  lm1 <- lm(log(mean_size) ~ log(std_linf), data = output_data)
  
  lm1_fit <-
    output_data %>% 
    add_predictions(lm1) %>% 
    ggplot(aes(x = std_linf, 
               y = mean_size)) +
    geom_point() +
    geom_line(aes(y = exp(pred)), col = "red") +
    theme_cowplot() +
    scale_x_continuous(label = label_number(suffix = "cm")) +
    scale_y_continuous(label = label_number(suffix = "cm")) +
    labs(x = "Asymptotic length (literature)", 
         y = "Mean length (observed)")
  
  ggsave(filename = "output/figures/ms_figs/linf_lmu.png",
         plot = lm1_fit,
         height = 10,
         width = 10*1.618,
         units = "cm")
  
  fig3_data <- 
    output_data %>% 
    left_join(mk_text) %>% 
    add_predictions(lm1, var = "pred_mean_size") %>% 
    mutate(mu = exp(pred_mean_size),
           sd = mu*0.34,
           sdlog = sqrt(log((0.34^2)+1)),
           meanlog = log(exp(pred_mean_size)) - ((sdlog^2)/2)) %>% 
    mutate(p_norm = 
             pnorm(size_max, mean = mu, sd = sd) -  
             pnorm(size_min, mean = mu, sd = sd),
           p_lnorm = 
             plnorm(size_max, meanlog = meanlog, sdlog = sdlog) -  
             plnorm(size_min, meanlog = meanlog, sdlog = sdlog)) %>% 
    select(species,
           size_class, 
           size_min, 
           size_max, 
           lbspr = p_lbspr, 
           lbspr_mk15 = p_lbspr_mk15,
           rls_obs = p_obs, 
           p_norm, 
           p_lnorm, 
           std_mk) %>% 
    filter(!(lbspr == 0 & rls_obs == 0)) %>% 
    pivot_longer(cols = c(lbspr, rls_obs, p_norm, p_lnorm, lbspr_mk15)) %>% 
    mutate(species = fct_reorder(species, std_mk)) 
  
  
  make_fig3 <- function(data, spp_select){
    dat1 <- 
      data %>% 
      filter(species == spp_select) 
    
    dat1 %>% 
      filter(name != "p_lnorm") %>% 
      filter(name != "rls_obs") %>% 
      # mutate(name = case_when(name == "lbspr" ~ "LBSPR",
      #                         name == "p_norm" ~ "Normal")) %>% 
      ggplot(aes(x = size_class, 
                 y = value)) +
      geom_rect(aes(xmin = size_min, 
                    xmax = size_max, 
                    ymin = 0, 
                    ymax = value, 
                    fill = size_10), 
                col = "black",
                data = dat1 %>% 
                  filter(name == "rls_obs") %>% 
                  mutate(size_10 = case_when(size_class >= 10 ~ "large", 
                                             TRUE ~ "small"))) +
      scale_fill_manual(values = c("large" = "grey70", 
                                   "small" = "grey90"),
                        labels = c("large" = "Observed (Length \U2265 10cm)", 
                                   "small" = "Observed (Length \U003C 10cm)")) +
      scale_color_manual(values = c("p_norm" = rgb(181, 144, 19, maxColorValue=255), 
                                    "lbspr" = rgb(148, 31, 36, maxColorValue=255), 
                                    "lbspr_mk15" = rgb(70, 172, 200, maxColorValue=255)),
                         labels = c("lbspr" = "LBSPR (estimate)", 
                                    "lbspr_mk15" = "LBSPR (M/k = 1.5, estimate)",
                                    "p_norm" = "Normal (estimate)")) +
      geom_vline(xintercept = 8.75, lty = 2) +
      geom_path(aes(col = name), linewidth = 2, alpha = 0.8) +
      geom_text(aes(x = Inf, y = Inf, label = paste("M/K =", std_mk, "\n", 
                                                    "L\U221E =", std_linf, "cm")),
                hjust = 1,
                vjust = 1.2,
                data = mk_text %>% filter(species == spp_select)) +
      scale_x_continuous(limits = c(0, NA), label = label_number()) +
      scale_y_continuous(label = label_percent()) +
      facet_wrap(~fct_reorder(species, std_mk), scales = "free") +
      labs(x = "Length class (cm)", 
           y = "Proportion") +
      theme_cowplot() +
      theme(legend.position = "none", 
            strip.background = element_rect(fill = "transparent"), 
            strip.text = element_text(face = "italic"), 
            axis.title.y = element_blank(), 
            axis.title.x = element_blank())
  }
  
  ks_data <- 
    fig3_data %>% 
    pivot_wider(names_from = name, values_from = value)
  
  ks_fig_data <- 
    ks_data %>% 
    summarise(ks_norm = ks.test(x = rls_obs, 
                                y = p_norm)$statistic,
              ks_lnorm = ks.test(x = rls_obs, 
                                 y = p_lnorm)$statistic,
              ks_lbspr = ks.test(x = rls_obs, 
                                 y = lbspr)$statistic,
              ks_lbspr_mk15 = ks.test(x = rls_obs, 
                                      y = lbspr_mk15)$statistic,
              .by = species) %>% 
    pivot_longer(cols = contains("ks_"))
  
  ks_fig_data_rmsmall <- 
    ks_data %>% 
    filter(size_class >= 10) %>% 
    summarise(ks_norm = ks.test(x = p_norm, 
                                y = rls_obs)$statistic,
              ks_lnorm = ks.test(x = p_lnorm, 
                                 y = rls_obs)$statistic,
              ks_lbspr = ks.test(x = lbspr, 
                                 y = rls_obs)$statistic,
              ks_lbspr_mk15 = ks.test(x = rls_obs, 
                                      y = lbspr_mk15)$statistic,
              .by = species) %>% 
    pivot_longer(cols = contains("ks_"))
  
  fig3_part2 <-
    ks_fig_data_rmsmall %>%
    rename(rmsmall = value) %>%
    left_join(ks_fig_data %>%
                rename(withsmall = value)) %>%
    pivot_longer(cols = contains("small"),
                 names_to = "inc_small") %>%
    mutate(inc_small = case_when(inc_small == "rmsmall" ~ "Exluding sizes < 10cm",
                                 inc_small == "withsmall" ~ "All size classes")) %>%
    # mutate(name = case_when(name == "ks_lbspr" ~ "LBSPR",
    #                         name == "ks_lnorm" ~ "Lognormal",
    #                         name == "ks_lbspr_mk15" ~ "LBSPR (M/k = 1.5)",
    #                         name == "ks_norm" ~ "Normal")) %>%
    filter(name != "ks_lnorm") %>%
    ggplot(aes(y = value, x = name)) +
    geom_violin()+
    geom_boxplot(width = 0.1) +
    facet_wrap(~inc_small, nrow=2) +
    scale_x_discrete(labels = c("ks_lbspr" = "LBSPR",
                                "ks_lbspr_mk15" = "LBSPR \n (M/k = 1.5)",
                                "ks_norm" = "Normal")) +
    labs(y = "Dissimilarity to observed") +
    theme_cowplot(20) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(colour = "grey80"),
          strip.background = element_rect(fill = "transparent", colour = "black"),
          axis.title.x = element_blank())
  , 
  axis.text.x = element_text(angle = 10, hjust = 1))

for(i in 1:length(unique(sample_output$species))){
  
  if(i==1) {
    out <- make_fig3(fig3_data, spp_select = unique(sample_output$species)[i]) + theme(legend.position = "none")
  } else if(i %in% c(6,7)) {
    out <- out + make_fig3(fig3_data, spp_select = unique(sample_output$species)[i]) + theme(axis.title.x = element_text())
  } else if(i == max(length(unique(sample_output$species)))) {
    out <- out + make_fig3(fig3_data, spp_select = unique(sample_output$species)[i]) + theme(legend.position = "bottom")
  } else {
    out <- out + make_fig3(fig3_data, spp_select = unique(sample_output$species)[i]) + theme(legend.position = "none")
  }
  
}

fig3 <- 
  out + 
  fig3_part2 + 
  {make_fig3(fig3_data, 
             spp_select = unique(sample_output$species)[i]) + 
      theme(legend.position = "bottom", legend.title = element_blank(), legend.justification = 0.5) } %>% 
  get_legend() %>% as_ggplot() +
  plot_layout(design = 
                "abc
              deh
              fgh
              iii
              ", heights = c(4,4,4,1)) +
  plot_annotation(tag_levels = "A")



ggsave(filename = fig3_filename,
       plot = fig3,
       height = 20,
       width = 20*1.618,
       units = "cm")


difference_plots_all_40cm <-
  output_data %>% 
  mutate(obs_minus_exp = p_obs-p_lbspr) %>% 
  ggplot(aes(x = size_class, 
             y = obs_minus_exp)) +
  geom_hline(yintercept = 0, col = "red") +
  geom_path(alpha = 0.1, aes(group = species)) + 
  geom_point(alpha = 0.1) +
  stat_smooth() +
  scale_x_continuous(breaks = c(output_data$size_class), 
                     limits = c(0, 40)) +
  scale_y_continuous(label = label_number())+
  labs(x = "Visual census length classes (cm)", 
       y = "Observed minus LBSPR") +
  annotate(geom = "text", x = Inf, 
           y = -Inf, 
           label = "Observed < LBSPR", 
           hjust = 1.1, 
           vjust = -0.5, 
           size = 8) +
  annotate(geom = "text", 
           x = Inf, 
           y = Inf, 
           label = "Observed > LBSPR", 
           hjust = 1.1, 
           vjust = 1.5, 
           size = 8) +
  theme_cowplot()


ggsave(filename = "output/figures/ms_figs/obs_minnus_lbspr.png",
       plot = difference_plots_all_40cm,
       height = 15,
       width = 15*1.618,
       units = "cm")

}





s# Figure 4 =====================================================================

set.seed(1)
random_select_species <- 
  plotdata_species %>% 
  filter(dat=="rls") %>% 
  mutate(ll_ratio = logl_50_norm-logl_50_lnorm) %>% 
  drop_na(ll_ratio) %>%
  mutate(pref = case_when(ll_ratio < quantile(ll_ratio, 0.025) ~ "lnorm_pref",
                          ll_ratio > quantile(ll_ratio, 0.975) ~ "norm_pref", 
                          ll_ratio < 2 & ll_ratio > -2 ~ "no_pref")) %>% 
  drop_na(pref) %>% 
  select(species, pref, ll_ratio) %>% 
  slice_sample(n=1, by = pref)

random_select_species_data <-
  obsdata_rls_species %>%
  filter(species %in% random_select_species$species) %>%
  left_join(meansizes_species) %>%
  select(species, mean_size, size_class, size_min, size_max, p_obs)


make_plot <- function(data, title){
  data %>% 
    ggplot(aes(x = size_class,
               y = p_obs, 
               width = size_max-size_min)) +
    facet_wrap(~species, ncol = 2, scales = "free") +
    geom_vline(aes(xintercept = mean_size), lty = 2, linewidth = 1.5) +
    geom_rect(aes(xmin = size_min, xmax = size_max, ymin = 0, ymax = p_obs), 
              alpha = 0.2, 
              col = "black",
              fill = "grey70") + 
    # geom_point(aes(col = name, y = value, alpha = alp)) +
    geom_path(aes(col = name, y = value, alpha = alp, lty=lty), linewidth = 2) +
    scale_alpha_identity()+
    scale_linetype_identity() +
    scale_colour_manual(values = c("p_norm" = rgb(181, 144, 19, maxColorValue=255), 
                                   "p_lnorm" = rgb(29, 84, 128, maxColorValue=255)))+
    theme_cowplot(20) +
    scale_y_continuous(label = label_percent())+
    scale_x_continuous(label = label_number()) +
    labs(x = "Length class (cm)", 
         y = "Proportion in length class") +
    theme(legend.position = "none",
          strip.background = element_rect(
            colour = "black",
            fill = "white"), 
          strip.text = element_text(face = "italic"), 
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(label = title)
}

fig4_data <-
  expand_grid(species = unique(random_select_species$species), 
  ) %>% 
  left_join(random_select_species_data %>% 
              select(species, mean_size, size_min) %>% 
              distinct(), 
            by = join_by(species)) %>% 
  left_join(rls_bin_table, 
            by = join_by(size_min)) %>% 
  left_join(random_select_species_data %>% 
              select(species, p_obs, size_class), 
            by = join_by(species, size_class)) %>% 
  left_join(random_select_species, 
            by = join_by(species)) %>% 
  arrange(pref) %>% 
  mutate(mu = mean_size,
         sd = mu*0.34,
         sdlog = sqrt(log((0.34^2)+1)),
         meanlog = log(mean_size) - ((sdlog^2)/2)) %>% 
  mutate(p_norm = 
           pnorm(size_max, mean = mu, sd = sd) -  
           pnorm(size_min, mean = mu, sd = sd),
         p_lnorm = 
           plnorm(size_max, meanlog = meanlog, sdlog = sdlog) -  
           plnorm(size_min, meanlog = meanlog, sdlog = sdlog)) %>% 
  pivot_longer(cols = c(p_norm, p_lnorm)) %>% 
  mutate(alp = ifelse(str_remove(pref, "_pref")==str_remove(name, "p_"), 1, 0.8),
         lty = ifelse(str_remove(pref, "_pref")==str_remove(name, "p_"), "solid", "73")) 

fig4a <- 
  fig4_data %>% 
  filter(pref == "lnorm_pref") %>% 
  make_plot("Strong preference for lognormal")

fig4b <- 
  fig4_data %>% 
  filter(pref == "no_pref") %>% 
  make_plot("No preference for normal or lognormal")

fig4c <- 
  fig4_data %>% 
  filter(pref == "norm_pref") %>% 
  make_plot("Strong preference for normal")

for(inc_fished in c(FALSE, TRUE)){
  
  fig4_filename <- paste0("output/figures/ms_figs/var_explained_", 
                          {if(inc_fished) "incfished" else "unfished"}, 
                          ".png")
  
  if(!file.exists(fig4_filename) | force_run){
    
    if(!file.exists(paste0("output/tables/var_explained_bycv",  
                           ifelse(inc_fished, "_incfished", "_unfished"), 
                           ".csv")) | FALSE){
      
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
      
      write_csv(out, 
                paste0("output/tables/var_explained_bycv_ks",  
                       ifelse(inc_fished, "_incfished", "_unfished"), 
                       ".csv"))
      rm(c, out, estimated_prob, lmer_norm, lmer_lnorm, lmer_pref, joined_data, d1, d2)
    }
    
    
    ks_curve <- 
      paste0("output/tables/var_explained_bycv_ks",  
             ifelse(inc_fished, "_incfished", "_unfished"), 
             ".csv") %>% 
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
           y = "Dissimilarity (K-S statistic)") +
      theme_cowplot(20) +
      theme(legend.position = c(0.95,0.05), 
            legend.justification = c(1,0),
            legend.title = element_blank(), 
            plot.background = element_rect(color = "transparent")) 
    
    fig4 <-  
      fig4a + 
      fig4b + 
      fig4c + 
      ks_curve +
      plot_annotation(tag_levels = 'A')
    
    ggsave(filename = fig4_filename,
           plot = fig4,
           height = 25,
           width = 25*1.618,
           units = "cm")
    
  }
}