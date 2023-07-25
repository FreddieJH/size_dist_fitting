# test script
# trying to create a new spatial scale that is more appropriate than the gridcell level


RLS_data_survs %>% 
  ggplot(aes(x = longitude, 
             y = latitude)) +
  geom_point(pch = 21) +
  geom_point(size = 0.05, col = "red") +
  geom_vline(xintercept = seq(100.5, 170.5, by = 1)) +
  geom_hline(yintercept = seq(-45.5, -8.5, by = 1)) +
  theme_cowplot()

RLS_data_survs %>% 
  ggplot(aes(x = longitude, 
             y = latitude)) +
  geom_point(pch = 21) +
  geom_point(size = 0.05, col = "red") +
  geom_vline(xintercept = seq(100.5, 170.5, by = 5)) +
  geom_hline(yintercept = seq(-45.5, -8.5, by = 5)) +
  theme_cowplot()
RLS_data_survs %>% 
  ggplot(aes(x = longitude, 
             y = latitude)) +
  geom_point(pch = 21) +
  geom_point(size = 0.05, col = "red") +
  geom_vline(xintercept = seq(100.5, 170.5, by = .5)) +
  geom_hline(yintercept = seq(-45.5, -8.5, by = .5)) +
  theme_cowplot()
  
# gridcell actually seems okay, maybe we could do smaller gridcells but have a larger cutoff

