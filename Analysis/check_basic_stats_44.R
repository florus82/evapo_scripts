library(tidyverse)

dat = read.csv('Y:/repos/evapo_scripts/Analysis/basic_stats_vrt_44.csv')

p <- dat %>%
  group_by(folder, comp, day, mask) %>% 
  pivot_longer(
    cols = c(maxLST, minLST, X25thLST, medianLST, X75thLST),
    names_to = 'stat',
    values_to = 'value'
  ) %>%  
  ggplot(aes(x = day, y = value, color=stat)) + 
  geom_point()
  # facet_grid(~ dat$comp + dat$day, scales = "free") + 
  # theme_minimal()'
  


print(p)
