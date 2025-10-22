library(tidyverse)

dat = read.csv('Y:/repos/evapo_scripts/Analysis/extract_evap0.csv')

dat_sub = dat %>% 
  filter(name_var == 'LST',
         LST_Mask == 'no',
         S2_Mask == 'yes',
         ET_var == 'Canopy'#,
         #date == '23-07-2019'
         ) %>% 
  mutate(across(starts_with('v'), as.numeric))

  ggplot(dat_sub, aes(x = factor(date), fill = Mask_Name)) + 
    geom_boxplot(
      stat = 'identity',
      aes(lower = vp25, upper = vp75, middle = vp50, ymin = vmin, ymax = vmax)
    ) + 
    facet_wrap(Tile + preds ~ composite, scales = 'free') +
    #ggtitle(label = ) + 
    #labs(x='Date in 2019', y = 'Temperature in Kelvin', fill='Preprocessing') + 
    theme_minimal() + 
    theme_gray()
    # theme(axis.title.x = element_text(margin = margin(t = 10), face = 'bold', size=14),
    #       axis.title.y = element_text(margin = margin(r = 10), face = 'bold', size=14),
    #       plot.title = element_text(margin = margin(b=10), hjust = 0.5, face = 'bold', size = 16),
    #       strip.text = element_text(face = 'bold', size = 12),
    #       legend.title = element_text(face = 'bold'),
    #       legend.background = element_rect(fill='lightgrey', color = 'black'))