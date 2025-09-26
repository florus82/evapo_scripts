library(tidyverse)

dat = read.csv('Z:/et/Landsat/LST_44/extract.csv')

dat %>% 
  sample_n(100000) %>% 
  ggplot(aes(x=Sentinel2, y=Landsat)) + 
  geom_point()
