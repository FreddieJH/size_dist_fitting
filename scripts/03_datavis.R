# About ========================================================================

# This script creates a plot with all the empirical body size distributions 

# Parameters ===================================================================

force_run <- FALSE

# Packages =====================================================================

library(tidyverse)
library(scales)
library(cowplot)

# Empirical distributions ======================================================

for(src in c("rls", "cbf")){
  for(pop in c("species", "ecoregion", "gridcell")){
    
    if(src == "cbf" & pop == "ecoregion") next
    pop <- ifelse(src == "cbf" & pop == "gridcell", "location", pop)
    
    current_obsdata <- 
      get(paste0("obsdata_", src, "_", pop)) %>% 
      arrange(population)
    
    dists <- current_obsdata$population %>% unique()
    plots_per_page <- ifelse(src == "cbf", 50, 100)
    
    for(page in 1:ceiling(length(dists)/plots_per_page)){
      
      first_dist <- ((page-1)*plots_per_page)+1
      second_dist <- first_dist + (plots_per_page-1)
      
      if(second_dist > length(dists)) {
        second_dist <- length(dists)
      }
      
      filename <- paste0("output/figures/data_vis/empirical_dists/", 
                         src, "/", 
                         pop, "/", 
                         first_dist, "_", 
                         second_dist, ".png")
      
      if(!file.exists(filename)|force_run){
        
        p <- 
          current_obsdata %>% 
          filter(population %in% dists[first_dist:second_dist]) %>% 
          ggplot() +
          aes(
            x = size_class, 
            y = n
          ) +
          geom_point() +
          geom_path() +
          facet_wrap(~population, scales = "free") +
          scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))
        
        ggsave(filename = filename, 
               plot = p, 
               height = 15, 
               width = 25)
        rm(p)
      }
    }
    rm(current_obsdata, plots_per_page, first_dist, second_dist, dists, page, pop)
  }
  rm(src)
}


# Bimodal distributions ========================================================

# Of the empirical distributions, plotted from the code above, bimodal 
# distributions were identified

for(src in c("rls", "cbf")){
  for(pop in c("species", "ecoregion", "gridcell")){
    
    if(src == "cbf" & pop == "ecoregion") next
    pop <- ifelse(src == "cbf" & pop == "gridcell", "location", pop)
    
    read_csv(file = paste0("input/data/data_cleaning/bimodal_", src, "_", pop, ".csv"), 
             show_col_types = FALSE) %>% 
      pull(population) %>% 
      assign(paste0("bimodal_pops_", src, "_", pop), value = ., envir = .GlobalEnv)
    
    filename <- paste0("output/figures/data_vis/bimodal_dists/", 
                       src, "_", 
                       pop, 
                       ".png")
    
    if(!file.exists(filename)|force_run){
      
      p <- 
        get(paste0("obsdata_", src, "_", pop)) %>% 
        filter(population %in% get(paste0("bimodal_pops_", src, "_", pop))) %>% 
        ggplot() +
        aes(
          x = size_class, 
          y = n
        ) +
        geom_point() +
        geom_path() +
        facet_wrap(~population, scales = "free") +
        theme_cowplot() +
        scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()))
      
      ggsave(filename = filename, 
             plot = p, 
             height = 15, 
             width = 25)
      
      rm(p)
    }
    rm(pop)
  }
  rm(src)
}

# End ==========================================================================