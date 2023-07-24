# ==============================================================================
#  FIGURES FOR THE MANUSCRIPT
# ==============================================================================

CBF_stan_model <- "CBF_mod02"
CBF_min_count  <- 10
CBF_population_level = "location"

CBF2_stan_model <- "CBF_mod13"
CBF2_min_bins  <- 4
CBF2_min_count  <- 10
CBF2_population_level = "location"

RLS_stan_model <- "RLS_mod13"
RLS_min_bins   <- 4
RLS_min_count  <- 500
RLS_population_level = "ecoregion"


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


CBF_obs_data_step1 <-
  paste0("CBF_data_obs_", CBF_population_level) %>%
  get() %>%
  filter(population_n >= CBF_min_count)

CBF_n_populations <-
  CBF_obs_data_step1 %>%
  pull(population) %>%
  n_distinct()

CBF_population_indx_table <-
  CBF_obs_data_step1 %>%
  select(population) %>%
  distinct() %>%
  mutate(population_indx = row_number())

CBF_obs_data <-
  CBF_obs_data_step1 %>%
  left_join(CBF_population_indx_table,
            by = join_by(population))



CBF_model_name  <- paste0(CBF_stan_model,
                          "_n", CBF_min_count, "_",
                          CBF_population_level,
                          CBF_n_populations)


CBF2_obs_data_step1 <-
  paste0("CBF_data_obs_", CBF2_population_level, "_binned") %>%
  get() %>%
  filter(n_sizebins >= CBF2_min_bins, 
         population_n >= CBF2_min_count)

CBF2_n_populations <-
  CBF2_obs_data_step1 %>%
  pull(population) %>%
  n_distinct()

CBF2_population_indx_table <-
  CBF2_obs_data_step1 %>%
  select(population) %>%
  distinct() %>%
  mutate(population_indx = row_number())

CBF2_obs_data <-
  CBF2_obs_data_step1 %>%
  left_join(CBF2_population_indx_table,
            by = join_by(population))



CBF2_model_name  <- paste0(CBF2_stan_model,
                           "_nbin", CBF2_min_bins,
                           "_n", CBF2_min_count, "_",
                           CBF2_population_level,
                           CBF2_n_populations)

RLS_model_par_raw <-
  paste0("analysis/output/model_pars/", RLS_model_name, ".parquet") %>%
  read_parquet() %>%
  summarise(value = mean(value),
            .by = c(population_indx, param)) %>%
  pivot_wider(values_from = value, names_from = param, names_prefix = "mean_")

CBF_model_par_raw <-
  paste0("analysis/output/model_pars/", CBF_model_name, ".parquet") %>%
  read_parquet() %>%
  summarise(value = mean(value),
            .by = c(population_indx, param)) %>%
  pivot_wider(values_from = value, names_from = param, names_prefix = "mean_")

CBF2_model_par_raw <-
  paste0("analysis/output/model_pars/", CBF2_model_name, ".parquet") %>%
  read_parquet() %>%
  summarise(value = mean(value),
            .by = c(population_indx, param)) %>%
  pivot_wider(values_from = value, names_from = param, names_prefix = "mean_")

RLS_dat1 <-
  paste0("analysis/output/likelihood_values/", RLS_model_name, ".parquet") %>%
  read_parquet()

CBF_dat1 <-
  paste0("analysis/output/likelihood_values/", CBF_model_name, ".parquet") %>%
  read_parquet()


RLS_point_size <-
  RLS_dat1 %>%
  left_join(RLS_obs_data %>% select(population_indx, population_n) %>% distinct(),
            by = join_by(population_indx)) %>%
  mutate(ll = ll/population_n) %>%
  pivot_wider(names_from = dist,
              values_from = ll) %>%
  mutate(diff = N - LN) %>%
  summarise(diff = median(diff),
            ll_norm = median(N),
            ll_lnorm = median(LN),
            .by = population_indx)

RLS_norm_better_n2 <-
  RLS_dat1 %>%
  pivot_wider(names_from = dist,
              values_from = ll) %>%
  mutate(norm_better = N>LN) %>%
  summarise(prop_norm_better = mean(norm_better),
            .by = population_indx) %>%
  mutate(norm_better = prop_norm_better > 0.5)


CBF_point_size <-
  CBF_dat1 %>%
  left_join(CBF_obs_data %>% select(population_indx, population_n) %>% distinct(),
            by = join_by(population_indx)) %>%
  mutate(ll = ll/population_n) %>%
  pivot_wider(names_from = dist,
              values_from = ll) %>%
  mutate(diff = N - LN) %>%
  summarise(diff = median(diff),
            ll_norm = median(N),
            ll_lnorm = median(LN),
            .by = population_indx)

CBF_norm_better_n2 <-
  CBF_dat1 %>%
  pivot_wider(names_from = dist,
              values_from = ll) %>%
  mutate(norm_better = N>LN) %>%
  summarise(prop_norm_better = mean(norm_better),
            .by = population_indx) %>%
  mutate(norm_better = prop_norm_better > 0.5)


plot_data <-
  RLS_model_par_raw %>%
  left_join(RLS_point_size, by = join_by(population_indx)) %>%
  left_join(RLS_norm_better_n2, by = join_by(population_indx)) %>%
  mutate(type = "RLS") %>%
  bind_rows(CBF_model_par_raw %>%
              left_join(CBF_point_size, by = join_by(population_indx)) %>%
              left_join(CBF_norm_better_n2, by = join_by(population_indx)) %>%
              mutate(type = "CBF"))

p1 <-
  plot_data %>%
  mutate(better_dist = ifelse(norm_better, "NORMAL PREFERED", "LOG-NORMAL PREFERED"),
         data_type_full = ifelse(type == "CBF", "Cryptobenthic (continuous)", "Reef Life Survey (binned)")) %>%
  ggplot() +
  aes(
    x = mean_meanlog,
    y = mean_sdlog,
    col = data_type_full,
    fill= data_type_full,
  ) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = 'y ~ x', alpha = 0.1) +
  scale_color_npg() +
  labs(
    x = "Lognormal meanlog",
    y = "Lognormal sdlog"
  ) +
  facet_wrap(.~better_dist, ncol = 1) +
  theme_cowplot(20) +
  theme(legend.position = "none",
        # legend.justification = c(1,1.2),
        legend.background = element_rect(colour = "transparent"),
        legend.title = element_blank())

p2 <-
  plot_data %>%
  mutate(better_dist = ifelse(norm_better, "NORMAL PREFERED", "LOG-NORMAL PREFERED"),
         data_type_full = ifelse(type == "CBF", "Cryptobenthic (continuous)", "Reef Life Survey (binned)")) %>%
  ggplot() +
  aes(
    x = mean_mu,
    y = mean_sigma,
    col = data_type_full,
    fill= data_type_full,
  ) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = 'y ~ x', alpha = 0.1) +
  scale_color_npg() +
  labs(
    x = "Normal \u03BC",
    y = "Normal \u03C3"
  ) +
  facet_wrap(.~better_dist, ncol = 1, scales = "free") +
  theme_cowplot(20) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1.2),
        legend.background = element_rect(colour = "transparent"),
        legend.title = element_blank())

p3 <-
  plot_data %>%
  mutate(better_dist = ifelse(norm_better, "NORMAL PREFERED", "LOG-NORMAL PREFERED"),
         data_type_full = ifelse(type == "CBF", "Cryptobenthic (continuous)", "Reef Life Survey (binned)")) %>%
  ggplot() +
  aes(
    x = mean_mu,
    y = mean_cv,
    col = data_type_full,
    fill= data_type_full,
  ) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = 'y ~ x', alpha = 0.1) +
  scale_color_npg() +
  labs(
    x = "Normal \u03BC",
    y = "Normal coeficient of Variation (\u03C3 / \u03BC)"
  ) +
  facet_wrap(.~better_dist, ncol = 1, scales = "free") +
  theme_cowplot(20) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1.2),
        legend.background = element_rect(colour = "transparent"),
        legend.title = element_blank())


p4 <- p1 + p2 + plot_annotation(tag_levels = 'A')
p5 <- p1 + p3 + plot_annotation(tag_levels = 'A')

ggsave(filename = paste0("ms_figures/param_regression.png"),
       plot = p4,
       # height = 10,
       # width = 10*1.618,
       units = "cm",
       dpi = 300)

ggsave(filename = paste0("ms_figures/param_regression_cov.png"),
       plot = p5,
       # height = 10,
       # width = 10*1.618,
       units = "cm",
       dpi = 300)

write_csv(x = plot_data, file = "ms_data/parameter_summary.csv")



p1 <-
  plot_data %>%
  mutate(better_dist = ifelse(norm_better, "NORMAL PREFERED", "LOG-NORMAL PREFERED"),
         data_type_full = ifelse(type == "CBF", "Cryptobenthic (continuous)", "Reef Life Survey (binned)")) %>%
  ggplot() +
  aes(
    x = mean_mu,
    y = mean_sigma,
    col = data_type_full,
    fill= data_type_full,
  ) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", formula = 'y ~ x', alpha = 0.1) +
  scale_color_npg() +
  labs(
    x = "Lognormal mu",
    y = "Lognormal sigma"
  ) +
  theme_cowplot(20) +
  theme(legend.position = "none",
        # legend.justification = c(1,1.2),
        legend.background = element_rect(colour = "transparent"),
        legend.title = element_blank())

# Figure 2 ---------------------------------------------------------------------

# overlay of species size distributions 


population_meansize <-
  RLS_data_obs_ecoregion %>% 
  uncount(n) %>% 
  summarise(mean_size_class = mean(size_class),
            max_size_class = max(size_class),
            .by = population) 

summary_tbl <- 
  read_parquet(paste0("analysis/output/model_pars/", RLS_model_name, ".parquet")) %>% 
  summarise(mean = mean(value), 
            median = median(value), 
            q5 = quantile(value, prob = 0.05), 
            q95 = quantile(value, prob = 0.95), 
            .by = param)

pop_indx_tbl <- 
  paste0("RLS_data_obs_", RLS_population_level) %>%
  get() %>%
  filter(n_sizebins >= RLS_min_bins,
         population_n >= RLS_min_count) %>%
  clean_data() %>% 
  select(population, population_indx) %>% 
  distinct()

xx <- 
  read_parquet(paste0("analysis/output/model_pars/", RLS_model_name, ".parquet")) %>%
  left_join(pop_indx_tbl) %>% 
  left_join(population_meansize)


xx2 <- tibble()
for(i in 1:max(pop_indx_tbl$population_indx)){
  
  xx3 <-
    xx %>% 
    filter(population_indx == i) %>%
    expand_grid(
      size = 1:max(population_meansize$max_size_class)
    ) %>% 
    filter(size < max_size_class) %>% 
    pivot_wider(names_from = param, 
                values_from = value) %>% 
    mutate(p_N = dnorm(size, mu, sigma), 
           p_LN = dlnorm(size, meanlog, sdlog)) %>% 
    summarise(mean_p_N = mean(p_N, na.rm = TRUE), 
              median_p_N = median(p_N, na.rm = TRUE),
              mean_p_LN = mean(p_LN, na.rm = TRUE), 
              median_p_LN = median(p_LN, na.rm = TRUE), 
              q5_p_N = quantile(p_N, probs = 0.05, na.rm = TRUE), 
              q5_p_LN = quantile(p_LN, probs = 0.05, na.rm = TRUE), 
              q95_p_N = quantile(p_N, probs = 0.95, na.rm = TRUE), 
              q95_p_LN = quantile(p_LN, probs = 0.95, na.rm = TRUE), 
              .by = c("population_indx", "population", "size")) 
  
  
  xx2 <- bind_rows(xx2, xx3)
  
}

write_csv(xx2, "ms_data/est_size_bypop.csv")



# xx <-
#   read_parquet(paste0("analysis/output/model_pars/", RLS_model_name, ".parquet")) %>%
#   left_join(pop_indx_tbl) %>% 
#   left_join(population_meansize) %>% 
#   expand_grid(
#     size = 1:max(population_meansize$max_size_class)
#   ) %>% 
#   filter(size < max_size_class)

scale_size <- function(size_vector, lower_limit = 1.25) {
  tibble(size = size_vector) %>% 
    filter(size > lower_limit) %>% 
    mutate(size = rls_bin(size)) %>% 
    mutate(mean_size = mean(size)) %>% 
    mutate(scaled_size = size/mean_size) %>% 
    count(scaled_size) %>% 
    mutate(scaled_n = n/sum(n)) %>% 
    select(scaled_size, scaled_n)
}

est_lines <- 
  bind_rows(
    scale_size(rnorm(1e6, 
                     summary_tbl %>% filter(param == "mu") %>% pull(mean), 
                     summary_tbl %>% filter(param == "sigma") %>% pull(mean))) %>% 
      mutate(dist = "norm", quant = "mean"),
    scale_size(rnorm(1e6, 
                     summary_tbl %>% filter(param == "mu") %>% pull(q5), 
                     summary_tbl %>% filter(param == "sigma") %>% pull(q5))) %>% 
      mutate(dist = "norm", quant = "q5"),
    scale_size(rnorm(1e6, 
                     summary_tbl %>% filter(param == "mu") %>% pull(q95), 
                     summary_tbl %>% filter(param == "sigma") %>% pull(q95))) %>% 
      mutate(dist = "norm", quant = "q95"),
    scale_size(rlnorm(1e6, 
                      summary_tbl %>% filter(param == "meanlog") %>% pull(mean), 
                      summary_tbl %>% filter(param == "sdlog") %>% pull(mean))) %>% 
      mutate(dist = "lnorm", quant = "mean"),
    scale_size(rlnorm(1e6, 
                      summary_tbl %>% filter(param == "meanlog") %>% pull(q5), 
                      summary_tbl %>% filter(param == "sdlog") %>% pull(q5))) %>% 
      mutate(dist = "lnorm", quant = "q5"),
    scale_size(rlnorm(1e6, 
                      summary_tbl %>% filter(param == "meanlog") %>% pull(q95), 
                      summary_tbl %>% filter(param == "sdlog") %>% pull(q95))) %>% 
      mutate(dist = "lnorm", quant = "q95")
  )


norm_better_tbl <- 
  paste0("analysis/output/likelihood_values/", RLS_model_name, ".parquet") %>%
  read_parquet() %>% 
  pivot_wider(values_from = ll, 
              names_from = dist) %>% 
  mutate(norm_better = N > LN) %>% 
  summarise(prop_norm_better = mean(norm_better),
            .by = population) %>% 
  mutate(norm_better = prop_norm_better > 0.5)


normalised_ll_ratio <-
  paste0("analysis/output/likelihood_values/", RLS_model_name, ".parquet") %>%
  read_parquet() %>%
  left_join(RLS_obs_data %>% select(population_indx, population_n) %>% distinct(),
            by = join_by(population_indx)) %>%
  mutate(ll = ll/population_n) %>%
  pivot_wider(names_from = dist,
              values_from = ll) %>%
  mutate(diff = N - LN) %>%
  summarise(diff = median(diff),
            ll_norm = median(N),
            ll_lnorm = median(LN),
            .by = population)

 

# p6 <- 
RLS_data_obs_ecoregion %>% 
  filter(population %in% norm_better_tbl$population) %>% 
  left_join(population_meansize) %>% 
  left_join(norm_better_tbl) %>% 
  left_join(normalised_ll_ratio) %>% 
  mutate(scaled_size = size_class/mean_size_class) %>% 
  mutate(scaled_n = n/population_n) %>%
  mutate(better_dist = ifelse(norm_better, "NORMAL PREFERED", "LOG-NORMAL PREFERED")) %>% 
  ggplot(aes(x = size_class)) +
  geom_line(aes(y = scaled_n, group = population), alpha = 0.3, col = "grey") +
  # geom_line(aes(y = scaled_n, group = population), col = "black", linewidth = 2, data = est_lines %>% filter(dist == "norm", quant == "mean")) +
  # geom_ribbon(aes(x = size, ymin = q5_p_N, ymax = q95_p_N), data = xx2) +
  # geom_line(col = "black", linetype = 2, data = est_lines %>% filter(dist == "norm", quant == "q5")) +
  # geom_line(col = "black", linetype = 2, data = est_lines %>% filter(dist == "norm", quant == "q95")) +
  # geom_line(col = "blue", linewidth = 2, data = est_lines %>% filter(dist == "lnorm", quant == "mean")) +
  # geom_line(col = "blue", linetype = 2, data = est_lines %>% filter(dist == "lnorm", quant == "q5")) +
  # geom_line(col = "blue", linetype = 2, data = est_lines %>% filter(dist == "lnorm", quant == "q95")) +
  # scale_x_log10() +
  # scale_y_log10() +
  facet_grid(~better_dist) +
  theme(legend.position = "none") +
  theme_cowplot(20) +
  labs(
    x = "Scaled body size (log)",
    y = "Scaled abundance (log)"
  )

ggsave(filename = paste0("ms_figures/scaled_bodysize.png"),
       plot = p6,
       # height = 10,
       # width = 10*1.618,
       units = "cm",
       dpi = 300)

traits <- read_csv("analysis/input/data/traits_2020.csv") %>% 
  rename(species_name = CURRENT_SPECIES_NAME) %>% 
  janitor::clean_names()

RLS_data_obs_ecoregion %>% 
  filter(population %in% norm_better_tbl$population) %>% 
  left_join(population_meansize) %>% 
  left_join(norm_better_tbl) %>% 
  left_join(normalised_ll_ratio) %>% 
  left_join(traits) %>% 
  ggplot(aes(x = thermal_guild, 
             y = diff)) +
  geom_point(alpha = 0.1) +
  stat_smooth()


load("analysis/input/data/fish_data.RData")
load("analysis/input/data/temporal_fish_data.RData")

# Figure 3 ---------------------------------------------------------------------
# Sdlog and CV estimates at various spatial scales

out <- tibble()
for(poplvl in c("species", "gridcell", "ecoregion")){
  RLS_stan_model <- "RLS_mod13"
  RLS_min_bins   <- 4
  RLS_min_count  <- 5000
  RLS_population_level = poplvl
  
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
  
  out <- 
    read_parquet(paste0("analysis/output/model_pars/", RLS_model_name, "_summary.parquet")) %>% 
    mutate(poplvl = poplvl, data = "RLS") %>% 
    bind_rows(out)
}

for(poplvl in c("species", "location")){
  
  CBF_stan_model <- "CBF_mod02"
  CBF_min_count  <- 10
  CBF_population_level = poplvl
  
  
  CBF_obs_data_step1 <-
    paste0("CBF_data_obs_", CBF_population_level) %>%
    get() %>%
    filter(population_n >= CBF_min_count)
  
  CBF_n_populations <-
    CBF_obs_data_step1 %>%
    pull(population) %>%
    n_distinct()
  
  CBF_population_indx_table <-
    CBF_obs_data_step1 %>%
    select(population) %>%
    distinct() %>%
    mutate(population_indx = row_number())
  
  CBF_obs_data <-
    CBF_obs_data_step1 %>%
    left_join(CBF_population_indx_table,
              by = join_by(population))
  
  
  
  CBF_model_name  <- paste0(CBF_stan_model,
                            "_n", CBF_min_count, "_",
                            CBF_population_level,
                            CBF_n_populations)
  
  out <- 
    read_parquet(paste0("analysis/output/model_pars/", CBF_model_name, "_summary.parquet")) %>% 
    mutate(poplvl = poplvl, data = "CBF") %>% 
    bind_rows(out)
}

out %>% 
  ggplot(aes(x = poplvl, 
             y  = mean_sdlog, 
             col = data)) +
  geom_boxplot()


out %>% 
  ggplot(aes(x = poplvl, 
             y  = mean_cv, 
             col = data)) +
  geom_boxplot()

check_indx <- 
  out %>% 
  filter(mean_cv > 2) %>% pull(population_indx)


RLS_stan_model <- "RLS_mod13"
RLS_min_bins   <- 4
RLS_min_count  <- 5000
RLS_population_level = "gridcell"

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

RLS_obs_data %>% 
  filter(population_indx %in% check_indx) %>% 
  ggplot(aes(x = size_class, y = p)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~population)


paste0("analysis/output/model_pars/", RLS_model_name, ".parquet") %>% 
  read_parquet() %>%
  filter(population_indx %in% check_indx) %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  full_join(RLS_obs_data %>% filter(population_indx %in% check_indx), 
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

# ==============================================================================
#  END
# ==============================================================================