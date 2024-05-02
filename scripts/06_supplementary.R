

# ED Figure 1 ==================================================================

set.seed(1)
selected_pops_unimod <-
  obsdata_rls_gridcell %>% 
  mutate(lat = str_extract(population, "(?<=__).*(?=_)") %>% 
           as.numeric(), 
         lon = str_extract(population, "(?<=__-?\\d{1,2}_).*$") %>% 
           as.numeric()) %>% 
  mutate(pop_label = paste0(species, " (", lat, "°, ", lon, "°)")) %>% 
  mutate(is_bimodal = population %in% bimodal_pops_rls_gridcell) %>% 
  filter(!is_bimodal) %>% 
  mutate(pop_id = cur_group_id(), .by = population) %>% 
  filter(pop_id %in% sample(x = 1:max(pop_id), 
                            size = 4)) 

selected_pops_bimod <-
  obsdata_rls_gridcell %>% 
  mutate(lat = str_extract(population, "(?<=__).*(?=_)") %>% 
           as.numeric(), 
         lon = str_extract(population, "(?<=__-?\\d{1,2}_).*$") %>% 
           as.numeric()) %>% 
  mutate(pop_label = paste0(species, " (", lat, "°, ", lon, "°)")) %>% 
  mutate(is_bimodal = population %in% bimodal_pops_rls_gridcell) %>% 
  filter(is_bimodal) %>% 
  mutate(pop_id = cur_group_id(), .by = population) %>% 
  filter(pop_id %in% sample(x = 1:max(pop_id), 
                            size = 4)) 

bimod_plot <- function(data, nrow){
  ggplot(data) +
    geom_rect(aes(
      xmin = size_min,
      xmax = size_max,
      ymax = n, 
      ymin = 0
    ), fill = "grey70", 
    col = "black") +
    facet_wrap(~pop_label, 
               nrow = nrow,
               scales = "free_y") + 
    scale_x_continuous(labels = label_number(suffix = "cm")) +
    labs(x = "Body length class", 
         y = "Count") +
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "white", colour = "black"), 
          plot.title = element_text(hjust = 0.5)) 
}

p1 <- bimod_plot(selected_pops_unimod, 4) + ggtitle("Unimodal examples (n = 4 of 3,169)")
p2 <- bimod_plot(selected_pops_bimod, 4) + ggtitle("Bimodal examples (n = 4 of 59)")

rls_bimod_dists <- 
  p1 + p2 +
  plot_layout(ncol= 2) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "output/figures/ms_figs/bimodal_rls.png", 
       plot = rls_bimod_dists, 
       height = 10, 
       width = 10)


# ED Figure 2 ==================================================================

set.seed(1)
selected_pops_unimod_cbf <-
  obsdata_cbf_location %>% 
  filter(population_n > 100) %>% 
  mutate(pop_label = paste0(species, " (", str_extract(population, "(?<=__).*"), ")")) %>% 
  mutate(is_bimodal = population %in% bimodal_pops_cbf_location) %>% 
  filter(!is_bimodal) %>% 
  mutate(pop_id = cur_group_id(), .by = population) %>% 
  filter(pop_id %in% sample(x = 1:max(pop_id), 
                            size = 3)) %>% 
  mutate(str_replace_all(pop_label, "cflatifasciata", "latifasciata"))
set.seed(1)
selected_pops_bimod_cbf <-
  obsdata_cbf_location %>% 
  mutate(pop_label = paste0(species, " (", str_extract(population, "(?<=__).*"), ")")) %>% 
  mutate(is_bimodal = population %in% bimodal_pops_cbf_location) %>% 
  filter(is_bimodal) %>% 
  mutate(pop_id = cur_group_id(), .by = population) %>% 
  filter(pop_id %in% sample(x = 1:max(pop_id), 
                            size = 3)) %>% 
  mutate(str_replace_all(pop_label, "cflatifasciata", "latifasciata"))

p1 <- bimod_plot(selected_pops_unimod_cbf, 3) + ggtitle("Unimodal examples (n = 3 of 133)")
p2 <- bimod_plot(selected_pops_bimod_cbf, 3) + ggtitle("Bimodal distributions (n = 3 of 3)")

cbf_bimod_dists <- 
  p1 + p2 +
  plot_layout(ncol= 2) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "output/figures/ms_figs/bimodal_cbf.png", 
       plot = cbf_bimod_dists, 
       height = 10, 
       width = 10)



# ED Figure 3 ==================================================================

pref_meansize_filename <- "output/figures/ms_figs/preference_meansize.png"

if(!file.exists(pref_meansize_filename) | force_run){
  
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
    geom_point(alpha = 0.3, pch = 3, size = 3) +
    geom_smooth(method = "glm",
                method.args = list(family = "binomial"),
                col = "darkblue", fill = "darkblue") +
    # geom_line(aes(y = predict(logit_mod, type = "response"))) +
    scale_x_continuous(label = label_number(suffix = "cm")) +
    # scale_x_log10(label = label_number(suffix = "cm")) +
    theme_cowplot(20) +
    labs(
      x = "Mean body length", 
      y = "Normal preferred?"
    )
  
  ggsave(filename = pref_meansize_filename,
         plot = logit_plot,
         height = 15,
         width = 15*1.618,
         units = "cm")
}

# ED Figure 4 ==================================================================


pop_data <- 
  plotdata_ecoregion %>% 
  left_join(meansizes_ecoregion) %>%
  filter(!bimodal) %>% 
  mutate(col = case_when(
    dat=="rls" & normal_better ~ "rls_norm",
    dat=="rls" & !normal_better ~ "rls_lnorm",
    dat=="cbf" & normal_better ~ "cbf_norm",
    dat=="cbf" & !normal_better ~ "cbf_lnorm",
  )) 



plot_ylim <- 
  max(layer_scales(
    main_plot(pop_data))$y$range$range)


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
           label = "Ecoregion-level",
           size = 10) +
  side_plot(spp_data %>% filter(col != "cbf_lnorm"), 
            plot_ylim)  +
  plot_layout(design = {"
      AAAAAAB
      AAAAAAB"
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
              legend.title = element_text()))) 

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
            #BBBBB#"
  })

ggsave(filename = "output/figures/ms_figs/meansize_cv_ecoregion.png",
       plot = p_simple,
       height = 20,
       width = 20*1.618,
       units = "cm")



# ED Figure 5 ==================================================================

# This figure is generated in the 05_figures.R script

# ED Figure 6 ==================================================================

# How does minimum mean count influence the results of the study?

mincount_tbl <- 
  tibble(data = "rls", 
         mincount = c(20, 50, 200, 500)) %>% 
  bind_rows(tibble(
    data = "cbf", 
    mincount = c(2, 5, 20, 50)
  ))

for(mod in c("lognormal", "normal")){
  
  for(i in 1:nrow(mincount_tbl)){
    
    out2 <- tibble()
    
    dat <- mincount_tbl$data[i]
    pop <- ifelse(dat == "rls", "gridcell", "location")
    mincount <- mincount_tbl$mincount[i]
    
    filename <- 
      paste0("output/tables/sensitivity_mincount_", 
             dat, "_", 
             mod, "_",
             mincount, ".csv")
    
    if(!file.exists(filename)){
      
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
            
            write_csv(out2, filename)
          } 
        }
        gc()
      }
    } else {
      print(paste(filename, "exists already."))
    }
  }
}

files <- 
  list.files("output/tables/", 
             pattern = "sensitivity", 
             full.names = TRUE) 

plotdat0 <- 
  map(files, read_csv, show_col_types = FALSE) %>% 
  bind_rows(.id = "id") %>% 
  mutate(species = str_extract(population, ".*(?=__)")) %>% 
  filter(!species %in% c("Sepia apama", 
                         "Sepioteuthis australis", 
                         "Sepia plangon", 
                         "Aipysurus laevis", 
                         # "Bathytoshia brevicaudata", # this is classed as elasmobranch
                         "Arctocephalus pusillus")) %>% # non-fish M1 species
  mutate(id = as.numeric(id)) %>% 
  left_join(tibble(filename = files) %>% mutate(id = row_number()), 
            by = join_by(id)) %>% 
  mutate(mincount = str_extract(filename, "\\d*(?=.csv)") %>% as.numeric()) %>% 
  summarise(mid = median(cv_50), 
            sd = sd(cv_50), 
            lwr = quantile(cv_50, probs = 0.05),
            upr = quantile(cv_50, probs = 0.95), 
            n = n(),
            .by = c(dat, mod, mincount)) %>% 
  mutate(se = sd/n) %>% 
  mutate(datatype = ifelse(dat == "rls", "In-situ observation", 
                           "Cryptobenthic specimen collection")) 
p0 <-
  plotdat0 %>% 
  mutate(mincount = as.factor(mincount)) %>% 
  mutate(n = max(n), .by = c(mincount, datatype)) %>% 
  ggplot(aes(x = mincount, y = mid, col = mod)) + 
  geom_linerange(aes(ymin = lwr, ymax = upr), 
                 position = position_dodge(0.2), 
                 linewidth = 2) +
  geom_point(position = position_dodge(0.2), 
             size = 4) +
  geom_text(aes(label = paste("n =",n)), y = 0.09, show.legend = FALSE, col = "black") +
  facet_wrap(datatype~., scales = "free_x") +
  scale_color_manual(values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
                                "lognormal" = rgb(29, 84, 128, maxColorValue=255)), 
                     labels = c("normal" = "Normal",
                                "lognormal" = "Lognormal")) +
  labs(x = "Minimum population size (# individuals)", 
       y = "Median CV estimate") +
  theme_cowplot(20) +
  theme(legend.position = c(1,1), 
        legend.justification = c(1,1.1),
        legend.title = element_blank(), 
        plot.margin = margin(10, 50, 10, 10), 
        strip.background = element_rect(fill = "white", colour = "black"))

ggsave(filename = "output/figures/ms_figs/mincount_sensitivity.png",
       plot = p0,
       height = 15,
       width = 25,
       units = "cm")

# ED Figure 7 ==================================================================

# This figure is generated in the 05_figures.R script

# ED Figure 8 ==================================================================

pdat <- 
  plotdata_species %>% 
  mutate(data = dat) %>%
  left_join(meansizes_species, by = join_by(species, data)) %>% 
  filter(!bimodal) 

p_out <- 
  pdat %>% 
  ggplot() + 
  aes(
    x = mean_size, 
    y = cov_pref, 
    pch = dat,
    # col = better_dist
  ) +
  geom_point(alpha = 1, size = 2, col = "grey80") +
  geom_point(alpha = 1, size = 2, col="black",
             data = pdat %>% 
               filter(species %in% c("Scorpis aequipinnis", 
                                     "Cheilodactylus spectabilis"))) +
  geom_point(alpha = 0.4, size = 2, col="red", 
             data = pdat %>% 
               filter(str_detect(species, "Eviota"))) +
  geom_point(alpha = 1, size = 2, col="red",
             data = pdat %>% 
               filter(str_detect(species, "Eviota queenslandica"))) +
  scale_x_continuous(trans = "log10", 
                     labels = label_number(suffix="cm")) +
  scale_shape_manual(
    values = c("rls" = 21, "cbf" = 24), 
    labels = c("rls" = "Visual census",
               "cbf" = "Cryptobenthic (exhaustive)")) +
  scale_color_manual(
    values = c("normal" = rgb(181, 144, 19, maxColorValue=255),
               "lognormal" = rgb(29, 84, 128, maxColorValue=255)),
    labels = c("normal" = "Normal preferred",
               "lognormal" = "Lognormal preferred")) +
  labs(x = "Mean body size (log)", 
       y = "Coefficient of variation in body size") +
  
  theme_cowplot(20) +
  theme(legend.position = "none"
  ) 

ggsave(filename = "output/figures/ms_figs/cv_meansize_spphighlight.png",
       plot = p_out,
       height = 20,
       width = 20*1.618,
       units = "cm")


# ED Figure 9 ==================================================================

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
        legend.text = element_text(size = 12, face = "italic")) +
  labs(x = "Species mean body length", 
       y = "Coefficient of Variation") 

ggsave(filename = "output/figures/ms_figs/intraspecies_cv_variation.png",
       plot = p,
       height = 15,
       width = 15*1.618,
       units = "cm")


# ED Figure 10 =================================================================

set.seed(1)
xx <- 
  plotdata_species %>% 
  filter(dat=="rls") %>% 
  mutate(ll_ratio = logl_50_norm-logl_50_lnorm) %>% 
  drop_na(ll_ratio) %>%
  mutate(pref = case_when(ll_ratio < quantile(ll_ratio, 0.025) ~ "lnorm_pref",
                          ll_ratio > quantile(ll_ratio, 0.975) ~ "norm_pref", 
                          ll_ratio < 2 & ll_ratio > -2 ~ "no_pref")) %>% 
  drop_na(pref) %>% 
  select(species, pref, ll_ratio) %>% 
  slice_sample(n=4, by = pref)

tab1 <-
  obsdata_rls_species %>%
  filter(species %in% xx$species) %>%
  left_join(meansizes_species) %>%
  select(species, mean_size, size_class, size_min, size_max, p_obs)


make_plot <- function(data, title){
  data %>% 
    ggplot(aes(x = size_class,
               y = p_obs, 
               width = size_max-size_min)) +
    facet_wrap(~species, ncol = 2, scales = "free") +
    geom_vline(aes(xintercept = mean_size), lty = 2) +
    geom_rect(aes(xmin = size_min, xmax = size_max, ymin = 0, ymax = p_obs), 
              alpha = 0.2, 
              col = "black",
              fill = "grey70") + 
    # geom_point(aes(col = name, y = value, alpha = alp)) +
    geom_path(aes(col = name, y = value, alpha = alp, lty=lty), 
              linewidth = 2) +
    scale_alpha_identity()+
    scale_linetype_identity() +
    scale_colour_manual(values = c("p_norm" = rgb(181, 144, 19, maxColorValue=255), 
                                   "p_lnorm" = rgb(29, 84, 128, maxColorValue=255)))+
    theme_cowplot(20) +
    scale_y_continuous(label = label_percent())+
    scale_x_continuous(label = label_number()) +
    labs(x = "Size class (cm)", 
         y = "Proportion in size class") +
    theme(legend.position="none",
          strip.background=element_rect(colour="black",
                                        fill="grey97"), 
          strip.text = element_text(face = "italic", size = 12), 
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(label = title)
}

p3_alldat <-
  expand_grid(species = unique(xx$species), 
              # size_min = sort(unique(tab1$size_min))
  ) %>% 
  left_join(tab1 %>% select(species, mean_size, size_min) %>% distinct()) %>% 
  left_join(rls_bin_table) %>% 
  left_join(tab1 %>% select(species, p_obs, size_class)) %>% 
  left_join(xx) %>% 
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
         lty = ifelse(str_remove(pref, "_pref")==str_remove(name, "p_"), 1, 2)) 

p3_1 <- 
  p3_alldat %>% 
  filter(pref == "lnorm_pref") %>% 
  make_plot("Strong preference for lognormal")
p3_2 <- 
  p3_alldat %>% 
  filter(pref == "no_pref") %>% 
  make_plot("No preference for normal or lognormal")
p3_3 <- 
  p3_alldat %>% 
  filter(pref == "norm_pref") %>% 
  make_plot("Strong preference for normal")

for(inc_fished in c(FALSE, TRUE)){
  
  fig3_filename_supple <- paste0("output/figures/ms_figs/var_explained_", 
                                 {if(inc_fished) "incfished" else "unfished"}, 
                                 "_extra.png")
  
  if(!file.exists(fig3_filename_supple) | force_run){
    
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
    
    p5 <-  
      p3_1 + 
      p3_2 + 
      p3_3 + 
      ks_curve +
      plot_annotation(tag_levels = 'A')
    
    ggsave(filename = fig3_filename_supple,
           plot = p5,
           height = 25,
           width = 25*1.618,
           units = "cm")
    
  }
}

# ED Figure 11 =================================================================

binned_effect <- 
  bind_rows(
    read_parquet("output/models/summary/cv/cbf_binned_location_normal.parquet") %>%
      select(population, cv_50, logl_50) %>% 
      mutate(dist = "norm", 
             binning = "binned"), 
    read_parquet("output/models/summary/cv/cbf_binned_location_lognormal.parquet") %>%
      select(population, cv_50, logl_50) %>% 
      mutate(dist = "lnorm", 
             binning = "binned"), 
    read_parquet("output/models/summary/cv/cbf_location_normal.parquet") %>%
      select(population, cv_50, logl_50) %>% 
      mutate(dist = "norm", 
             binning = "cont"), 
    read_parquet("output/models/summary/cv/cbf_location_lognormal.parquet") %>%
      select(population, cv_50, logl_50) %>% 
      mutate(dist = "lnorm", 
             binning = "cont")
  ) %>% 
  pivot_wider(names_from = dist, values_from = c(cv_50, logl_50)) %>% 
  mutate(cv_pref = case_when(is.na(logl_50_norm) ~ cv_50_lnorm, 
                             is.na(logl_50_lnorm) ~ cv_50_norm, 
                             logl_50_norm > logl_50_lnorm ~ cv_50_norm,
                             logl_50_norm < logl_50_lnorm ~ cv_50_lnorm)) %>% 
  select(population, binning, cv_pref) %>% 
  pivot_wider(values_from = cv_pref, names_from = binning) %>% 
  left_join(meansizes_gridcell)

p_binning <- 
  binned_effect %>%
  ggplot() +
  aes(
    x = mean_size,
    xend = mean_size,
    y = cont,
    yend = binned
  ) +
  geom_hline(yintercept = median(binned_effect$binned, na.rm = TRUE), lty = 1) +
  geom_hline(yintercept = median(binned_effect$cont, na.rm = TRUE), lty = "35") +
  geom_segment(arrow = arrow(length = unit(0.3, "cm"))) +
  scale_x_log10(label = label_number(suffix = "cm")) +
  theme_cowplot(20) +
  labs(x = "Mean body length (log)", 
       y = "CV value (before and after binning)")


ggsave(filename = "output/figures/ms_figs/binning_effect.png",
       plot = p_binning,
       height = 15,
       width = 15*1.618,
       units = "cm")



# ED Figure 12 =================================================================

# This figure is generated in the 05_figures.R script

# ED Figure 13 =================================================================

# This figure is generated in the 05_figures.R script

# ED Figure 14 =================================================================
set.seed(1)
fig_ed14_data <- 
  obsdata_rls_gridcell %>% 
  mutate(species_id = cur_group_id(),
         .by = population) %>% 
  filter(population_n > 500) %>% 
  filter(species_id %in% sample(1:max(species_id), 6)) %>% 
  mutate(lat = str_extract(population, "(?<=__).*(?=_)") %>% 
           as.numeric(), 
         lon = str_extract(population, "(?<=__-?\\d{1,2}_).*$") %>% 
           as.numeric()) %>% 
  mutate(pop_label = paste0(species, "\n(", lat, "°, ", lon, "°)")) 

fig_ed14_p1 <- 
fig_ed14_data %>% 
  ggplot(aes(x = size_class, y = n)) +
  geom_point() +
  geom_path() +
  facet_wrap(~pop_label, 
             scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "white", colour = "black"))+ 
  labs(x = "Length class (cm)", 
       y = "Count") 

fig_ed14_p2 <- 
fig_ed14_data %>% 
  mutate(mass_min = (size_min^3)*0.01, 
         mass_max = (size_max^3)*0.01, 
         mass_mid = (size_class^3)*0.01) %>% 
  mutate(mass_width = mass_max-mass_min) %>% 
  mutate(norm_n = n/mass_width) %>% 
  ggplot(aes(x = mass_mid, y = norm_n)) +
  geom_point() +
  geom_path() +
  scale_x_continuous(labels = label_number(scale_cut = cut_si("g"))) +
  facet_wrap(~pop_label, scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "white", colour = "black"))+ 
  labs(x = "Weight class", 
       y = "Count")

fig_ed14 <- 
  fig_ed14_p1 + 
  fig_ed14_p2 +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "output/figures/ms_figs/biomass_dists.png",
       plot = fig_ed14,
       height = 20,
       width = 20,
       units = "cm")


# Other information ============================================================


plotdata_gridcell %>% 
  pull(species) %>% 
  unique() %>% 
  length()

plotdata_species %>% 
  pull(species) %>% 
  unique() %>% 
  length()

obsdata_cbf_species %>% 
  pull(species) %>% 
  c(obsdata_rls_species %>% 
      pull(species)) %>% 
  unique() %>% 
  length()



# median CV values -------------------------------------------------------------
plotdata_gridcell %>% 
  filter(!bimodal) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1))

plotdata_gridcell %>% 
  drop_na(cov_pref) %>% 
  filter(!bimodal) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05) %>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1),
            .by = dat)

plotdata_gridcell %>% 
  filter(!bimodal) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05) %>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1),
            .by = better_dist)

plotdata_ecoregion %>% 
  filter(!bimodal) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1))

plotdata_gridcell %>% 
  filter(!bimodal) %>%
  filter(!fished) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1))

plotdata_species %>% 
  filter(!bimodal) %>%
  filter(!fished) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1))

plotdata_gridcell %>%
  filter(!bimodal) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1))

plotdata_species %>% 
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1))

plotdata_gridcell %>%
  filter(!bimodal) %>%
  summarise(cv_05 = quantile(cov_pref, probs = 0.05)%>% round(2),
            cv_10 = quantile(cov_pref, probs = 0.1)%>% round(2),
            cv_50 = quantile(cov_pref, probs = 0.5)%>% round(2),
            cv_90 = quantile(cov_pref, probs = 0.9)%>% round(2),
            cv_95 = quantile(cov_pref, probs = 0.95)%>% round(2), 
            cv_sd = sd(cov_pref)%>% round(2),
            cv_n = n(),
            se = (cv_sd/sqrt(cv_n))%>% signif(digits = 1), 
            .by = c(dat, normal_better))

# Range of max sizes -----------------------------------------------------------

obsdata_rls_gridcell %>% 
  bind_rows(obsdata_cbf_location) %>% 
  summarise(max_size = max(size_class), 
            .by = population) %>% 
  pull(max_size) %>% 
  range()

# Population sizes -------------------------------------------------------------
get_npops <- Vectorize(function(src, pop, has_cv, rm_bimodal){
  
  temp_pop <- ifelse(pop == "location", "gridcell", pop)
  cbf_pop  <- ifelse(pop %in% c("ecoregion", "gridcell"), "location", pop)
  
  if(src == "all"){
    
    
    {if(!has_cv) bind_rows(get(paste("obsdata", "rls", pop, sep = "_")), 
                           get(paste("obsdata", "cbf", cbf_pop, sep = "_"))) else get(paste("plotdata", pop, sep = "_"))} %>% 
      {if(rm_bimodal) filter(., !bimodal) else .} %>% 
      mutate(species = str_extract(population, ".*(?=__)")) %>% 
      filter(!species %in% c("Sepia apama", 
                             "Sepioteuthis australis", 
                             "Sepia plangon", 
                             "Aipysurus laevis", 
                             # "Bathytoshia brevicaudata", # this is classed as elasmobranch
                             "Arctocephalus pusillus")) %>% 
      pull(population) %>% 
      unique() %>% 
      length()
    
  } else {
    
    {if(!has_cv) get(paste("obsdata", src, pop, sep = "_")) else get(paste("plotdata", temp_pop, sep = "_")) %>% filter(dat == src)} %>% 
      mutate(species = str_extract(population, ".*(?=__)")) %>% 
      filter(!species %in% c("Sepia apama", 
                             "Sepioteuthis australis", 
                             "Sepia plangon", 
                             "Aipysurus laevis", 
                             # "Bathytoshia brevicaudata", # this is classed as elasmobranch
                             "Arctocephalus pusillus")) %>% 
      {if(rm_bimodal) filter(., !bimodal) else .} %>% 
      pull(population) %>% 
      unique() %>% 
      length()
    
  }
})



bimodal_pops_rls %>% length()

if(!file.exists("output/tables/population_sizes.csv") | force_run){
  tibble(src = c("rls", "all")) %>%
    expand_grid(pop = c("species", "ecoregion", "gridcell")) %>% 
    bind_rows(tibble(src = "cbf", 
                     pop = c("species", "location"))) %>% 
    expand_grid(has_cv = c(TRUE, FALSE)) %>% 
    expand_grid(rm_bimodal = c(TRUE, FALSE)) %>% 
    filter(!(!has_cv & rm_bimodal)) %>%
    mutate(popsize = get_npops(src, pop, has_cv, rm_bimodal)) %>% 
    write_csv("output/tables/population_sizes.csv")
}

obsdata_rls_gridcell$population %>% unique() %>% length()
obsdata_cbf_location$population %>% unique() %>% length()
obsdata_rls_ecoregion$population %>% unique() %>% length()
obsdata_rls_species$population %>% unique() %>% length()
obsdata_cbf_species$population %>% unique() %>% length()


# intra-species CV estimates

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
        legend.text = element_text(size = 12, face = "italic")) +
  labs(x = "Species mean size", 
       y = "Coefficient of Variation") 

ggsave(filename = "output/figures/ms_figs/intraspecies_cv_variation.png",
       plot = p,
       height = 15,
       width = 15*1.618,
       units = "cm")


# Convergence numbers ========================================================

err_files <- 
  list.files(path = "output/models/diagnostics", 
             pattern = "_errors", 
             full.names = TRUE) 

err_files %>% 
  map(read_csv, show_col_types = FALSE) %>% 
  bind_rows(.id = "id") %>% 
  mutate(id = as.numeric(id)) %>% 
  left_join(tibble(name = str_extract(err_files, "(?<=diagnostics/)(.*)(?=_errors\\.csv)"), 
                   id = 1:length(err_files))) %>% 
  separate(col = name, into = c("dat", "pop", "dist")) %>% 
  count(dat, pop, dist)

plotdata_gridcell %>% filter(is.na(mu_50))
plotdata_gridcell %>% filter(is.na(meanlog_50))
plotdata_gridcell 

# how many populations were removed due to too many indivs in first size bin?
obsdata_rls_gridcell %>% 
  arrange(population, size_class) %>% 
  mutate(min = min(size_class), .by = population) %>% 
  filter(min==size_class&p_obs > 0.5)

obsdata_cbf_location %>% 
  arrange(population, size_class) %>% 
  mutate(min = min(size_class), .by = population) %>% 
  filter(min==size_class&p_obs > 0.5)

summary_data <- 
  bind_rows(obsdata_rls_gridcell %>% mutate(dat = "rls"), 
            obsdata_cbf_location%>% mutate(dat = "cbf")) %>% 
  mutate(bimodal = population %in% c(bimodal_pops_rls, bimodal_pops_cbf)) %>% 
  select(population, dat, fished, bimodal) %>% 
  distinct() %>% 
  left_join(plotdata_gridcell %>% 
              mutate(norm_fit = is.numeric(mu_50), lnorm_fit = is.numeric(meanlog_50)) %>% 
              select(population, norm_fit, lnorm_fit)) 

count(summary_data, dat)

count(summary_data, fished)
count(summary_data, fished, dat)
count(summary_data, bimodal, dat)

count(summary_data, norm_fit)



# Binning CRF data =============================================================




obsdata_rls_gridcell %>% 
  mutate(species_id = cur_group_id(),
         .by = population) %>% 
  filter(species_id < 10) %>% 
  ggplot(aes(x = size_class, y = n)) +
  geom_point() +
  geom_path() +
  facet_wrap(~population, scales = "free")


obsdata_rls_gridcell %>% 
  mutate(mass_min = (size_min^3)*0.01, 
         mass_max = (size_max^3)*0.01, 
         mass_mid = (size_class^3)*0.01) %>% 
  mutate(mass_width = mass_max-mass_min) %>% 
  mutate(norm_n = n/mass_width) %>% 
  mutate(species_id = cur_group_id(), 
         .by = population) %>% 
  filter(species_id < 20) %>% 
  ggplot(aes(x = mass_mid, y = norm_n)) +
  geom_point() +
  geom_path() +
  scale_y_continuous(labels = label_number()) +
  facet_wrap(~population, scales = "free") +
  labs(x = "Body mass bin", 
       y = "Abundance density")
