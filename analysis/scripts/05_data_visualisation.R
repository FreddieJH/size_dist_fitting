# ==============================================================================
#  DATA VISUALISATION
# ==============================================================================

# RLS data visualisation -------------------------------------------------------

if(!file.exists("analysis/output/data_figs/ssd_obs_20spp.png")){
  
  p <-
    RLS_data_obs_sample %>% 
    count(species_name, size_class) %>% 
    ggplot() +
    geom_point(mapping = aes(x = size_class, y = n)) +
    geom_line(mapping = aes(x = size_class, y = n)) +
    facet_wrap( ~ species_name, scales = "free") +
    labs(x = "Fish length (cm)",  
         y = "Fish abundance") +
    theme_bw() +
    ggtitle("Body size distribution of 20 randomly sampled species")
  
  save_plot(plot = p, 
            filename = "analysis/output/data_figs/ssd_obs_20spp.png", 
            base_height = 8)
  
} 

# CBF data visualisation -------------------------------------------------------




# ==============================================================================
#  END
# ==============================================================================