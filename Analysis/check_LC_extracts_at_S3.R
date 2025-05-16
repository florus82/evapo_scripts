library(tidyverse)

# used fields
path  <- 'Z:/et/Auxiliary/landcover/csv'

files = list.files(path,'.csv', full.names = T)

conti = map_dfr(files, read_csv)



aa = conti %>% 
  group_by(row_col) %>% # goes across all landcovertypes
  summarise(across(where(is.numeric), ~sum(.x, na.rm = TRUE), .names = '{.col}_sum')) %>% 
  filter(`0_sum` != 76*76*5)

col_names = c("1110_sum", "1120_sum", "1130_sum", "1140_sum", "1150_sum", "1210_sum", "1220_sum",
  "1310_sum", "1320_sum", "1410_sum", "1420_sum", "1430_sum", "1440_sum", "2100_sum", 
  "2200_sum", "2310_sum", "2320_sum", "3100_sum", "3200_sum", "1_sum", "2_sum", "3_sum",
  "4_sum", "253_sum", "254_sum")


