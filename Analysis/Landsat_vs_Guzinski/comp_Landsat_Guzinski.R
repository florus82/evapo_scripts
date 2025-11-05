library(tidyverse)
library(lubridate)
# files = list.files('Y:/repos/evapo_scripts/Analysis/Landsat_vs_Guzinski/', full.names = T)
# block = map_dfr(files, read_csv)


# block %>% 
#   ggplot(aes(x=Guzinski, y=Landsat)) + 
#   geom_point() + 
#   facet_wrap(~Date) + 
#   geom_abline(intercept=0, slope=1, color = 'red') +
#   coord_fixed()

block = read.csv('Y:/repos/evapo_scripts/Analysis/Landsat_vs_Guzinski/daily_extracts.csv')

block$Date = as.Date(parse_date_time(block$Date, orders = 'Y_m_d'))
block$Month = format(block$Date, "%Y-%m")

block %>% 
  filter(Guzinski < 6,
         Month != '2019-08') %>% 
  ggplot(aes(x=Guzinski, y=Landsat)) + 
  # geom_point() + 
  geom_hex() + 
  scale_fill_viridis_c(option = 'plasma') +
  facet_wrap(Month~comp) +
  geom_abline(intercept=0, slope=1, color = 'red') +
  coord_fixed()
  


  

