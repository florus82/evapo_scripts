library(tidyverse)

# used fields
path  <- 'Z:/et/Auxiliary/landcover/csv'

files = list.files(path,'.csv', full.names = T)
conti = map_dfr(files, read_csv)



aa = conti %>% 
  filter(if_any(c(No_cropland, No_forest, No_grassland, No_impervious, Dry), ~ .x != 76 * 76)) 
  #group_by(row_col) %>% # goes across all landcovertypes
  # need some disolve here per row_col, so that all columns are filled
  
View(aa[1,])
  

