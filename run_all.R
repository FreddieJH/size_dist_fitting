library(tidyverse)
library(scales)
library(cowplot)
library(rstan)
library(posterior)
library(stringr)
library(arrow)
library(ggsci)
library(patchwork)
library(ggpubr)


source("scripts/00_setup.R")
source("scripts/01_import.R")
source("scripts/02_wrangling.R")
source("scripts/03_datavis.R")
source("scripts/04_modelling.R")
source("scripts/05_figures.R")
source("scripts/06_supplementary.R")



spp <- obsdata_rls_species %>% 
  pull(species)

xx <- rfishbase::popchar(species_list = spp) %>% 
  select(species = Species, Lmax)


obsdata_rls_species %>% 
  left_join(xx %>% summarise(lmax = median(Lmax, na.rm = T), .by = species)) %>% 
  drop_na() %>% 
  left_join(meansizes_species %>% filter(data == "rls")) %>% 
  ggplot(aes(lmax, mean_size)) +
  geom_point(alpha = 0.1) +
  # scale_x_log10() +
  # scale_y_log10() +
  geom_abline(slope = 1) +
  geom_line(aes(y = 3.347359+lmax*0.4))


obsdata_rls_species %>% 
  left_join(xx %>% summarise(lmax = median(Lmax, na.rm = T), .by = species)) %>% 
  drop_na() %>% 
  left_join(meansizes_species %>% filter(data == "rls")) %>% 
  lm(mean_size~0+lmax, data = .) %>% 
  summary()
