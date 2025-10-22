library(tidyverse)

files = list.files('Y:/repos/evapo_scripts/Analysis/Landsat_vs_Guzinski/', full.names = T)

block = map_dfr(files, read_csv)

block %>% 
  ggplot(aes(x=Guzinski, y=Landsat)) + 
  geom_point() + 
  facet_wrap(~Date) + 
  geom_abline(intercept=0, slope=1, color = 'red') +
  coord_fixed()
  

