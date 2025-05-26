library(tidyverse)

dat = read.csv('Z:/et/Landsat/composites/Brandenburg/field_stats_3ids.csv')

View(dat)
dat$datetime <- as.Date(ISOdate(dat$Year, dat$month, 1))

dat_long = pivot_longer(
  dat,
  cols = c(band_mean, band_std, band_min, band_max),
  names_to = "key",
  values_to = "value"
)


ggplot(dat_long, aes(x=datetime, y = value, color =key)) +
  geom_point() + 
  facet_wrap(~zone_id, scales = 'free_y', ncol = 1) + 
  theme_minimal()
  

ggplot(dat, aes(x = datetime, y = band_mean, color=band_mean))+#, color = as.factor(zone_id))) +
  geom_point(shape = 21, size = 3, show.legend = FALSE) +
  geom_errorbar(aes(ymin = band_mean - band_std, ymax = band_mean + band_std), width = 0.0, show.legend = FALSE) +
  facet_wrap(~zone_id, scales = 'free_y', ncol = 1) + 
  theme_minimal() + 
  labs(x='Time of observation',
       color = "Field ID",
       y = "Mean ± SD",
       title = "Summary statistics of your selected fields") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
  scale_color_distiller(palette = "GnBu", direction = 1) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust=0.5, face='bold')
    )
