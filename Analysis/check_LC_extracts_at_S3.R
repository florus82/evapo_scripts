library(tidyverse)

# used fields
path  <- 'Z:/et/Auxiliary/landcover/csv'

# files = list.files(path,'.csv', full.names = T)
# conti = map_dfr(files, read_csv)
# 
# dat = conti %>%
#   group_by(row_col) %>%
#   select(-LC) %>% 
#   summarise(across(
#     .cols = where(is.numeric),
#     .fns = ~ mean(.x, na.rm = TRUE)
#   ), .groups = 'drop') %>% 
#   filter(if_any(c(No_cropland, No_forest, No_grassland, No_impervious, Dry), ~ .x != 76 * 76)) 
# 
# write.csv(dat, paste0(path, '/filtered_merged.csv'), row.names = FALSE)

dat = read.csv(paste0(path, '/filtered_merged.csv'))
dat$Crop_prop = round((((76*76) - dat$No_cropland) / (76*76))*100)
dat$Forest_prop = round((((76*76) - dat$No_forest) / (76*76))*100)
dat$Grassland_prop = round((((76*76) - dat$No_grassland) / (76*76))*100)
dat$Impervious_prop = round((((76*76) - dat$No_impervious) / (76*76))*100)
dat$Water_prop = round((((76*76) - dat$Dry) / (76*76))*100) # do here only permanent water bodies!!!!

ind100 = dat$row_col[which(dat$Crop_prop==100)]
ind90 = dat$row_col[which(dat$Crop_prop>=90)]
ind80 = dat$row_col[which(dat$Crop_prop >= 80)]
ind70 = dat$row_col[which(dat$Crop_prop >= 70)]
ind60 = dat$row_col[which(dat$Crop_prop >= 60)]
ind50 = dat$row_col[which(dat$Crop_prop >= 50)]


# Kürzere Vektoren mit NA auffüllen
ind100 <- c(ind100, rep(NA, length(ind50) - length(ind100)))
ind90 <- c(ind90, rep(NA, length(ind50) - length(ind90)))
ind80 <- c(ind80, rep(NA, length(ind50) - length(ind80)))
ind70 <- c(ind70, rep(NA, length(ind50) - length(ind70)))
ind60 <- c(ind60, rep(NA, length(ind50) - length(ind60)))


# In ein data.frame umwandeln
df <- data.frame(
  'Thresh100' = ind100,
  'Thresh90' = ind90,
  'Thresh80' = ind80,
  'Thresh70' = ind70,
  'Thresh60' = ind60,
  'Thresh50' = ind50,
  stringsAsFactors = FALSE
)

# Als CSV exportieren
write.csv(df, paste0(path, '/row_cols.csv'), row.names = FALSE)

















dat$tot_sum_with_water = dat$Crop_prop + dat$Forest_prop + dat$Grassland_prop + dat$Impervious_prop + dat$Water_prop
dat$tot_sum_without_water = dat$Crop_prop + dat$Forest_prop + dat$Grassland_prop + dat$Impervious_prop

dat %>%
  pivot_longer(cols = c(tot_sum_with_water, tot_sum_without_water),
               names_to = "Variable", values_to = "Value") %>% 
  ggplot(aes(x = Value, fill = Variable)) +
  geom_histogram(binwidth = 5, alpha = 0.6, position = "identity", color = "black") +
  facet_wrap(~ Variable) +
  labs(title = "Histograms of Landcover Variables", x = "Value", y = "Count") +
  theme_minimal()


sum(dat$tot_sum_without_water > 100) / nrow(dat)
sum(dat$tot_sum_with_water > 100) / nrow(dat)



pure_crop = dat[dat$Crop_prop == 100,]
sum(pure_crop$Grassland_prop > 0 ) /nrow(pure_crop)
sum(pure_crop$Forest_prop > 0) / nrow(pure_crop)
sum(pure_crop$Impervious_prop > 0) / nrow(pure_crop)
sum(pure_crop$Water_prop > 0) / nrow(pure_crop)
# there seems to be quite a bit of water within pure crop pixel
sum(pure_crop$Permanent_water > 0) / nrow(pure_crop)
sum(pure_crop$Temporary_water > 0) / nrow(pure_crop)
sum(pure_crop$Permanent_wet > 0) / nrow(pure_crop)
sum(pure_crop$Temporary_wet > 0) / nrow(pure_crop)
sum(pure_crop$Sea_water > 0) / nrow(pure_crop)
sum(pure_crop$unclassifiable > 0) / nrow(pure_crop)
# the issue are 'temporary wet areas'. These are often areas close to rivers (bay + inundation) and close to the coast


dat$Water_prop_non_perm = round((dat$Permanent_water + dat$Permanent_wet + dat$unclassifiable + dat$Sea_water) / (76*76)*100)

dat %>%
  filter(Crop_prop > 80 & Crop_prop < 90) %>% 
  pivot_longer(cols = c(Crop_prop, Forest_prop, Grassland_prop, Impervious_prop, Water_prop_non_perm),
             names_to = "Variable", values_to = "Value") %>% 
  ggplot(aes(x = row_col, y = Value, color = Variable)) +
  #geom_histogram(binwidth = 5, alpha = 0.6, position = "identity", color = "black") +
  geom_point() + 
  theme_minimal()





