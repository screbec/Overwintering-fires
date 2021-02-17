### script for temporal drivers analysis ###
## including:
## - script for figure 2
## - script for ED figures 6 & 7

library(tidyverse)
library(egg) #ggplot layouts (tag_facet)

### Load data ------------------------------------------------------------------
# output from overwintering algorithm
data = read_csv('overwintering_out.csv')

# annual burned area
area_yr_reg = data %>% select(Year, ShortNm, FireID, area) %>% distinct() %>% 
  group_by(Year, ShortNm) %>% summarise(area = sum(area, na.rm = TRUE))

# mjjas maximum temperature
temp = read.csv('temp_max_extremes.csv',)

# regional average of burn depth (used in scatter plots)
depth_reg = read_csv('spatial.csv') %>% 
  left_join(data) %>% filter(!is.na(depth_mean), !is.na(Year)) %>% 
  mutate(depth_mean = depth_mean*n_depth) %>%
  group_by(Year, ShortNm) %>% summarise(depth_mean = sum(depth_mean)/sum(n_depth))
# this mean depth calculation is a bit complicated because it includes data
# but it is more consistent, because some fire scars on the border of NT and AK dont have their ignition in NT/AK
# and we remove these by linking the fireids to the ignition points

# burned area since 1975
ba_trend = read_csv('burned_area_trend.csv')

# combine all 
ignyr = data %>% filter(ign == 'overwinter') %>% 
  group_by(Year, ShortNm) %>% summarise(ignNo = n()) %>% ungroup() %>% 
  full_join(area_yr_reg) %>%
  full_join(depth_reg) %>%
  mutate(ignNo = ifelse(is.na(ignNo) & Year > 2001 & Year < 2019 , 0, ignNo)) %>%
  group_by(ShortNm) %>% arrange(Year) %>% mutate(ignNo_lag = lead(ignNo))  %>% full_join(temp)

# extended time series for time series/trend figure
ba_temp = ba_trend %>% pivot_longer(cols = 1:2, names_to = 'ShortNm', values_to = 'ba') %>% 
  full_join(ignyr)

### Subplots (Figure 2a-d) ------------------------------
## Figure 2 consist of three graph parts stitched together in illustrator
## and two tables with correlations

## part 1: time series plot ignitions ~ burned area
p = ggplot(ignyr) + theme_classic() + facet_wrap(~ShortNm) +
  # add area bars
  geom_area(aes(Year, area/10000), col = 'white', alpha = 0.5, fill = 'gray60') +
  # add overwintering ignition bars
  geom_col(aes(Year, ignNo_lag*mult), fill = 'firebrick', width = 0.3) +
  # add second axis
  scale_y_continuous(sec.axis = sec_axis(~ ./0.5, name = '# Fires that overwintered'), labels = comma) +
  # configure axes
  xlim(c(2000.5, 2017.5)) +
  labs(x = '', y = 'Annual burned area (Mha)') +
  theme(axis.line.y.right = element_line(color = "Firebrick"), axis.ticks.y.right = element_line(color = "Firebrick"),
        axis.text.y.right = element_text(color = "Firebrick"), axis.title.y.right = element_text(color = "Firebrick"), 
        axis.line.y.left = element_line(color = "gray40"), axis.ticks.y.left = element_line(color = "gray40"),
        axis.text.y.left = element_text(color = "gray40"), axis.title.y.left = element_text(color = "gray40"),
        text=element_text(size=7),
        aspect.ratio = 4/6) 
tag_facet(p, open = '', close = '', tag_pool = c('c','d'))

## part2: time series MJJAS temperature
p = ggplot(filter(ignyr, param == 'mean')) + theme_classic() + facet_wrap(~ShortNm) +
  geom_line(aes(Year, temp_perc), col = 'tomato') +
  labs(x='', y = 'MJJAS Temperature (K)') + ylim(c(286, 291)) + xlim(c(2001, 2017)) +
  theme(axis.line.y.left = element_line(color = "tomato"), axis.ticks.y.left = element_line(color = "tomato"),
        axis.text.y.left = element_text(color = "tomato"), axis.title.y.left = element_text(color = "tomato"),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        text=element_text(size=7),
        aspect.ratio = 4/6)
tag_facet(p, open = '', close = '', tag_pool = c('',''))

## part 3: burned area and temperature trend
# ggplot has no default for plotting a second y-axis
# MJJAS temperature values thus have to be scaled to range of burned area values (temp_scale)
# and the axis scaled back to temperature values
ba_temp2 = ba_temp %>% filter(is.na(param)| param == 'mean', Year < 2020) %>%
  mutate(ba = ba/10000)  # convert from km2 to Mha
ba_temp2 = mutate(ba_temp2, temp_scale = (temp_perc-min(ba_temp2$temp_perc, na.rm = T))/
                    (max(ba_temp2$temp_perc, na.rm = T) - min(ba_temp2$temp_perc, na.rm = T)) * max(ba_temp2$ba, na.rm = T))

# create annotation dataframe with slopa and p-values of trends
lm_t_AK = lm(temp_perc ~ Year, data = filter(ba_temp2, ShortNm == 'AK'))
lm_t_NT = lm(temp_perc ~ Year, data = filter(ba_temp2, ShortNm == 'NT'))
lm_ba_AK = lm(ba ~ Year, data = filter(ba_temp2, ShortNm == 'AK'))
# there is no burned area trend in Northwest Territories
tablestats = data.frame(slope = c(paste('+', round(summary(lm_ba_AK)$coefficients[2,1],4)*10000, '~km^2/yr'), ''),
                        p_reg = c(paste('p =', round(summary(lm_ba_AK)$coefficients[2,4],2)), ''), 
                        slopet = c(paste('+', round(summary(lm_t_AK)$coefficients[2,1],2), '~K/yr'), 
                                   paste('+', round(summary(lm_t_NT)$coefficients[2,1], 2), '~K/yr')),
                        p_regt = c(paste('p =', round(summary(lm_t_AK)$coefficients[2,4],2)), 
                                   paste('p =', round(summary(lm_t_NT)$coefficients[2,4], 2))),
                        x_slope = rep(2020, 2), y_slopet = c(1.75, 2.55), y_slope = c(0.65, 0.7))

p = ggplot(ba_temp2) + theme_classic() + facet_wrap(~ShortNm) + xlim(c(1975, 2030)) +
  # burned area time series
  geom_line(aes(x = Year, y = ba), col = 'gray40', linetype = 1, size = 0.4) +
  # trend line for burned area
  geom_smooth(data = filter(ba_temp2, ShortNm == 'AK'), aes(x = Year, y = ba),
              col = 'gray40', fill = 'gray40', size = 0.4, method = lm, linetype = 2) +
  # temperature time series
  geom_line(aes(x = Year, y = temp_scale), col = 'tomato', linetype = 1, size = 0.4) +
  # trend lines for temperature
  geom_smooth(aes(x = Year, y = temp_scale), 
              col = 'tomato', fill = 'tomato', size = 0.4, method = 'lm', linetype = 2) +
  # add shading
  geom_rect(aes(xmin = -Inf, xmax = 2001, ymin = -Inf, ymax = Inf), alpha = 0.01, fill = 'gray70') +
  geom_rect(aes(xmin = 2018, xmax = Inf, ymin = -Inf, ymax = Inf), alpha = 0.01, fill = 'gray70') + 
  # add annotations for slope and their p-values
  geom_text(data = anno, aes(x = x_slope, y = y_slope, label = slope), parse = TRUE, size = 1.7, hjust = 'left', col = "gray40") +
  geom_text(data = anno, aes(x = x_slope, y = y_slope, label = p_reg), size = 1.7, hjust = 'left', col = "gray40", nudge_y = -0.22) +
  geom_text(data = anno, aes(x = x_slope, y = y_slopet, label = slopet), parse = TRUE, size = 1.7, hjust = 'left', col = "tomato") +
  geom_text(data = anno, aes(x = x_slope, y = y_slopet, label = p_regt), size = 1.7, hjust = 'left', col = "tomato", nudge_y = -0.22) +
  # second y-axis scaled back to display actual temperature values
  scale_y_continuous(sec.axis = sec_axis(~ .*(max(ba_temp2$temp_perc, na.rm = T) - min(ba_temp2$temp_perc, na.rm = T))/
                                           max(ba_temp2$ba, na.rm = T) + min(ba_temp2$temp_perc, na.rm = T),
                                         name = 'MJJAS Temperature (K)'), labels = comma) +
  # axis labels & colours
  labs(y = 'Annual burned area (Mha)', x = '') +
  theme(axis.line.y.left = element_line(color = "gray40"), axis.ticks.y.left = element_line(color = "gray40"),
        axis.text.y.left = element_text(color = "gray40"), axis.title.y.left = element_text(color = "gray40"), 
        axis.line.y.right = element_line(color = "tomato"), axis.ticks.y.right = element_line(color = "tomato"),
        axis.text.y.right = element_text(color = "tomato"), axis.title.y.right = element_text(color = "tomato"),
        axis.text.x=element_text(color = c(rep('black', 5), 'transparent')),
        axis.ticks.x = element_line(color = c(rep('black', 5), 'transparent')),
        text=element_text(size=7),
        aspect.ratio = 4/6)
tag_facet(p, open = '', close = '', tag_pool = c('a','b'))


### Tables (Figure 2e-f): correlations between temporal drivers ---------------------------------
# computes the data for tables e and f of figure 2. Tables are created in Word

# burned area and MJJAS max temperature
group_by(ba_temp2, ShortNm) %>% 
  summarise(cor = round(cor(ba, temp_perc, 'complete.obs', method = 'spearman'),2), 
            pcor = cor.test(ba, temp_perc, method = 'spearman')$p.value)

# ignition and burned area
group_by(ignyr, ShortNm) %>% 
  summarise(cor = round(cor(ignNo_lag, area, 'complete.obs', method = 'spearman'),2), 
            pcor = cor.test(ignNo_lag, area, method = 'spearman')$p.value)

# ignition and MJJAS max temperature
group_by(ignyr, ShortNm) %>% 
  summarise(cor = round(cor(ignNo_lag, temp_perc, 'complete.obs', method = 'spearman'),2), 
            pcor = cor.test(ignNo_lag, temp_perc, method = 'spearman')$p.value)

## Extended Data Figures 5 & 6: Scatterplots of correlations -----------------------------

# function for plotting scatterplots with correlations & p-values as annotations
scatter_cor = function(data, varX, varY, xanno, yanno, labx, laby, figure_tags, scale_varX = 1, scale_varY = 1, xlim = NA){
  annos = data.frame(data %>% group_by(ShortNm) %>% 
                       summarise(cor = round(cor(get(varX), get(varY), 'complete.obs', method = 'spearman'),2), 
                                 pcor = cor.test(get(varX), get(varY), method = 'spearman')$p.value),
                     x_cor = xanno, y_cor = yanno) %>%
    mutate(pcor = ifelse(pcor < 0.001, 'p < 0.001', ifelse (pcor < 0.01, 'p < 0.01', paste('p =', signif(pcor, 1)))))
  p = ggplot(data) + theme_bw() +
    geom_point(aes(x = get(varX)/scale_varX, y = get(varY)/scale_varY, col = ShortNm)) + facet_wrap(~ShortNm) +
    geom_text(data = annos, aes(x = x_cor, y = y_cor, label = paste("rho == ", cor)), parse = T, size = 2.5, hjust = 'left') +
    geom_text(data = annos, aes(x = x_cor, y = y_cor, label = pcor), size = 2.5, hjust = 'left') +
    scale_color_manual(values = c('chartreuse4', 'chartreuse3'), guide = F) +
    labs(x = labx, y = laby) +
    theme(text=element_text(size=7),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 5/6)
  if (!is.na(xlim)){
    p = p + xlim(xlim)
  }
  print(tag_facet(p, open = '', close = '', tag_pool = figure_tags))
}

# ED Figure 6a-b: burned area vs mean mjjas tmax
scatter_cor(filter(ba_temp, param == 'mean'), varX = 'temp_perc', varY = 'ba', 
            xanno = rep(285.2, 2), yanno = rep(3.4, 2),
            labx = 'MJJAS Tmean (K)', laby = 'Annual burned area (Mha)',
            figure_tags = c('a','b'), scale_varY = 10000)
# ED Figure 6c-d: ign vs area
scatter_cor(filter(ba_temp, param == 'mean', Year > 2001), varX = 'ba', varY = 'ignNo_lag', 
            xanno = rep(0.3, 2), yanno = rep(7, 2),
            labx = 'Annual burned area (Mha)', laby = '# Overwintering flare-ups',
            figure_tags = c('c','d'), scale_varX = 10000)
# ED Figure 6e-f: ign vs mean mjjas tmax
scatter_cor(filter(ba_temp, param == 'mean'), varX = 'temp_perc', varY = 'ignNo_lag', 
            xanno = rep(286.3, 2), yanno = rep(7, 2),
            labx = 'MJJAS Tmean (K)', laby = '# Overwintering flare-ups',
            figure_tags = c('e','f'), xlim = c(286, 290.5))

# ED Figure 7a-b: burned area vs hot days90 mjjas tmax
scatter_cor(filter(ba_temp, param == 'hot_days90'), varX = 'temp_perc', varY = 'ba', 
            xanno = rep(10, 2), yanno = rep(3.4, 2), 
            labx = '# Hot days (T90)', laby = 'Annual burned area (Mha)',
            figure_tags = c('a','b'), scale_varY = 10000)
# ED Figure 7c-d: Overwintering fires vs hot days90 mjjas tmax
scatter_cor(filter(ba_temp, param == 'hot_days90'), varX = 'temp_perc', varY = 'ignNo_lag', 
            xanno = rep(10.5, 2), yanno = rep(7, 2), 
            labx = '# Hot days (T90)', laby = '# Overwintering flare-ups',
            figure_tags = c('c','d'))
# ED Figure 7e-f: Burn deptj vs mjjas tmax 90th percentile
scatter_cor(filter(ba_temp, param == '90'), varX = 'temp_perc', varY = 'depth_mean', 
            xanno = rep(292.6, 2), yanno = rep(20, 2), 
            labx = 'MJJAS Tmax90 (K)', laby = 'Average burn depth (cm)',
            figure_tags = c('e','f'))
