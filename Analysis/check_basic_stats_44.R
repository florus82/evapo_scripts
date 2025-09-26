library(tidyverse)

dat = read.csv('Y:/repos/evapo_scripts/Analysis/basic_stats_vrt_44.csv')

p <- dat %>%
  mutate(mask_fac = factor(mask, levels = c('S2notMasked_withoutLSTmask',
                                            'S2notMasked_withLSTmask',
                                            'S2Masked_withoutLSTmask',
                                            'S2Masked_withLSTmask'),
                              labels = c('no masking',
                                         'S3 masked', 'S2 masked', 'S2+S3 masked'))) %>% 
  
  mutate(comp_fac = factor(comp, levels = c('minVZA', 'maxLST'),
                           labels = c('Composite: minimum viewing angle',
                                      'Composite: maximum LST'))) %>% 
  
  mutate(folder_fac = factor(folder, levels = c('S2only', 'all_preds'),
                         labels = c('only S2 as predictor', 'S2, slope, aspect, solar incidence as predictors'))) %>% 
  ggplot(aes(x = factor(day), fill = mask_fac)) + 
  geom_boxplot(
    stat = 'identity',
    aes(lower = X25thLST, upper = X75thLST, middle = medianLST, ymin = minLST, ymax = maxLST)
  ) + 
  facet_wrap(comp_fac ~ folder_fac) +
  ggtitle('Effect of different input datasets on resulting sharpened LST values') + 
  labs(x='Date in 2019', y = 'Temperature in Kelvin', fill='Preprocessing') + 
  theme_minimal() + 
  theme_gray() +
  theme(axis.title.x = element_text(margin = margin(t = 10), face = 'bold', size=14),
        axis.title.y = element_text(margin = margin(r = 10), face = 'bold', size=14),
        plot.title = element_text(margin = margin(b=10), hjust = 0.5, face = 'bold', size = 16),
        strip.text = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold'),
        legend.background = element_rect(fill='lightgrey', color = 'black'))
  


print(p)




