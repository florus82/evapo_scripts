library(tidyverse)
library(lubridate)
library(ggpmisc)

# load csvs
landsat = read.csv('Y:/repos/evapo_scripts/Analysis/Flux_vs_Landsat/Landsat_extract.csv')
zalf = read.csv('X:/Evapotranspiration/flux_data/Raw/zalf.ID_5881_KC_LYSIMETER_raw.csv', header=T)

# clean up
colnames(zalf)[1] = 'MeasureDate'
landsat$date_standard = as.Date(landsat$date)
zalf$date_standard = as.Date(ymd_hms(zalf$MeasureDate))


# Filter and plot
zalf_sub = zalf %>%
  filter(date_standard >= min(landsat$date_standard) & 
           date_standard <= max(landsat$date_standard)) %>%
  pivot_longer(cols = starts_with('ETa'), 
               names_to = 'metric', 
               values_to = 'value') %>% 
  select(c(date_standard, metric, value))


landsat_sub = landsat %>% 
  rename(value = val) %>% 
  mutate(metric = 'Landsat')


zalf_sub %>% 
  bind_rows(landsat_sub) %>% 
  filter(value > 0, !is.na(value)) %>% 
  ggplot(aes(x = date_standard, y = value, color = metric)) +
  geom_line() +
  theme_light() +
  labs(x = "Date", y = "ETa value", title = "ETa over time")



# Calculate RMSE per facet
rmse_df <- zalf_sub %>%
  left_join(landsat_sub, by = "date_standard") %>%
  filter(value.x > 0, value.y > 0) %>%
  group_by(metric.x, metric.y) %>%
  summarise(
    rmse = sqrt(mean((value.y - value.x)^2)),
    x_pos = min(value.x),   # positions for the text
    y_pos = max(value.y)
  )


zalf_sub %>% 
  left_join(landsat_sub, by = 'date_standard') %>% 
  filter(value.x > 0,
         value.y > 0) %>% 
  ggplot(aes(x=value.x, y=value.y)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F, color = 'blue') +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = '~~~')),
               formula = y ~ x, parse = T) +
  geom_text(data = rmse_df,
            aes(x = x_pos, y = y_pos, label = paste("RMSE =", round(rmse, 2))),
            inherit.aes = FALSE, hjust = 0) +
  facet_wrap(metric.x ~ metric.y) +
  geom_abline(intercept=0, slope=1, color = 'red') + 
  coord_fixed() + 
  theme_light() +
  labs(x = "Zalf", y = "Guzinski")
