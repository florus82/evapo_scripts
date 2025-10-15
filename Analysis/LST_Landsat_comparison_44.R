library(tidyverse)
library(broom)
library(ggpubr)

# files = list.files('Z:/et/Landsat/LST_44/extracts/', full.names = T, pattern = 'no_agromask')
# dat_conti = list()
# 
# for (i in seq_along(files)){
#   
#   dat <- read.csv(files[i])
# 
#   file_id <- basename(files[i])
#   
# 
#   dat_conti[[i]] <- dat %>%
#     mutate(source_file = file_id) %>% 
#     sample_n(min(10000, nrow(dat)))   
#   
# }
# 
# 
# dat_all <- bind_rows(dat_conti)
# write.csv(dat_all, 'Z:/et/Landsat/LST_44/extracts/sampled_data_no_agro.csv', row.names = FALSE)
dat_all = read.csv('Z:/et/Landsat/LST_44/extracts/sampled_data_no_agro.csv')

dat1 = dat_all[dat_all$source_file == 'extract_allpred_maxLST_July_26_S2notMasked_withoutLSTmask_warped_no_agromask.csv',]




dat_sub = dat_all %>% 
  group_by(source_file) %>% 
  filter(str_detect(source_file, 'maxLST')) %>% 
  mutate(sourci = factor(source_file,
                         levels = c('extract_S2only_maxLST_July_26_S2notMasked_withoutLSTmask_warped_no_agromask.csv', 
                         'extract_S2only_maxLST_July_26_S2notMasked_withLSTmask_warped_no_agromask.csv',
                         'extract_allpred_maxLST_July_26_S2notMasked_withoutLSTmask_warped_no_agromask.csv', 
                         'extract_allpred_maxLST_July_26_S2notMasked_withLSTmask_warped_no_agromask.csv'),
                         labels = c('only S2 as predictor, no masking', 
                                    'only S2 as predictor, S3 masking', 
                                    'S2, slope, aspect, solar incidence as predictors, no masking',
                                    'S2, slope, aspect, solar incidence as predictors, S3 masking'))) %>% 
  sample_n(10000)
 





ggplot(dat_sub, aes(x=Sentinel2, y=Landsat)) +
geom_hex(bins = 20) + 
scale_fill_viridis_c(option = 'plasma') + 
#geom_point() + 
geom_abline(slope=1, intercept = 0, color = 'red', linetype = 'dashed') +
#geom_smooth(method = 'lm', color='black') + 
stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 295, label.y = 325) + 
stat_regline_equation(label.x = 295, label.y = 323) +
facet_wrap(~sourci) +
xlim(275,325) + 
ylim(275,325) + 
ggtitle('Sharpened Sentinel LST vs Landsat LST on 2019/07/26 (n=10000)',) + 
labs(x='Sharpened Sentinel LST', y = 'Landsat LST') + 
theme_minimal() + 
coord_equal() + 
theme_grey() + 
theme(axis.title.x = element_text(margin = margin(t = 10), face = 'bold', size=14),
      axis.title.y = element_text(margin = margin(r = 10), face = 'bold', size=14),
      plot.title = element_text(margin = margin(b=10), hjust = 0.5, face = 'bold', size = 16),
      strip.text = element_text(face = 'bold', size = 10),
      legend.title = element_text(face = 'bold'),
      legend.background = element_rect(fill='lightgrey', color = 'black'))


for (i in unique(dat_sub$sourci)){
  m = lm(dat_sub$Landsat[dat_sub$sourci == i] ~ dat_sub$Sentinel2[dat_sub$sourci == i])
  print(sqrt(mean((m$residuals)^2)))
}


m = lm(dat_sub$Landsat[dat_sub$sourci == i] ~ dat_sub$Sentinel2[dat_sub$sourci == i])

