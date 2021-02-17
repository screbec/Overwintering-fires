### script for spatial drivers analysis ###
## including:
## - script for figure 3
## - script for inputs to Table 1

library(tidyverse)
library(egg) #ggplot layouts (tag_facet)


### load data -----------------------------------------------------------------------------------------

# environmental variables including burn depth
spat_vars_mean =  read_csv('spatial.csv') %>% 
  # take each mean value * number of pixels with that value
  mutate(depth_mean = depth_mean*n_depth, # burn depth
         tc = tc * n_tc,     # tree cover
         sp_bs = sp_bs*n_ts, # black spruce cover
         sp_ws = sp_ws*n_ts, # white spruce cover
         sp_dc = sp_dc*n_ts, # deciduous cover
         sp_p = sp_p*n_ts,   # jack pine cover
         Elev = Elev*n, slope = slope*n, # elevation and slope
         soil30 = soil30*n_6) %>% # soil carbon in upper 30 cm
  # summarise for each fire scar (fire perimeters can be fragmented due to satellite detection)
  group_by(FireID, ign2, ShortNm) %>% summarise_all(.funs = sum, na.rm = TRUE) %>% 
  # divide by total pixels of fire scar to get percentage for whole fire scar
  mutate(depth_mean = depth_mean/n_depth, 
         tc = tc/n_tc, 
         sp_bs = sp_bs/n_ts, sp_ws = sp_ws/n_ts,
         sp_dc = sp_dc/n_ts, sp_p = sp_p/n_ts, 
         Elev = Elev/n, slope = slope/n,
         soil30 = soil30/n_6) %>% 
  mutate(depth_mean = ifelse(depth_mean == 0, NA, depth_mean)) %>% # where depth is 0, no burning happened
  # select variables of interest
    select(FireID, ShortNm, ign2, depth_mean, sp_bs, sp_ws, sp_dc, sp_p, Elev, slope, soil30, tc) %>% ungroup()

# for burn depth: filter out other dominance than black spruce
burnDepth = spat_vars_mean %>% filter(sp_ws + sp_dc + sp_p < 0.9) %>% 
  filter(!is.na(depth_mean)) %>% mutate(burn_depth = depth_mean*(-1))


### Table 1: t-tests for environmental variables -------------------------------------------------------

## part 1: t-tests for burn depth
# Welch's t-test
for (state in c('AK', 'NT')){
  x = filter(burnDepth, ShortNm == state, ign2 == 'overwinter') %>% select(burn_depth) %>% unlist()
  y = filter(burnDepth, ShortNm == state, ign2 == 'other') %>% select(burn_depth) %>% unlist()
  t = t.test(x, y, alternative = "two.sided")
  print(t)
  print(paste(sd(na.omit(x)), sd(na.omit(y)), length(x)))
}

## part 2: t-test for all other environmental variables
# reshape to long format so we can loop through all variables
test = spat_vars_mean %>% select(-FireID, -depth_mean) %>% gather(key = "variable", value = "value", -one_of('ign2', 'ShortNm'))
# compute t-tests and output means, standard deviations and t-test p-value in a table
for (state in c('NT', 'AK')){
  res = NULL
  dattest = filter(test, ShortNm == state)
  for (var in unique(test$variable)){
    x = filter(dattest, ign2 == 'overwinter', variable == var) %>% select(value) %>% unlist()
    y = filter(dattest, ign2 == 'other', variable == var) %>% select(value) %>% unlist()
    ttest = t.test(x, y, alternative = 'two.sided')
    res = rbind(res, data.frame(var, mean(na.omit(x)), sd(na.omit(x)), 
                                  mean(na.omit(y)), sd(na.omit(y)), round(ttest$p.value, 5)))
  }
  print(state)
  print(res)
}


### Figure 3: boxplots burn depth ---------------------------------------------------

# compute means per group (overwintering fires/other fires)
means = group_by(burnDepth, ShortNm, ign2) %>% summarise(mean = mean(burn_depth))
# create dataset for annotations
anno = means %>% spread(ign2, mean)

# plot
p = ggplot(burnDepth) + theme_bw() + facet_wrap(~ShortNm, ncol = 2) +
  geom_boxplot(aes(x = ign2, y = burn_depth, group = ign2, fill = ign2), 
               alpha = 0.7, lwd = 0.2, outlier.size = 0.1) +
  geom_point(data = means, aes(x = ign2, y = mean, group = ign2, fill = ign2), shape = 3, size = 0.5) +
  geom_text(data = anno, aes(x = Inf, y = -Inf, label = sprintf('mu == "%.1f"', round(overwinter, 1))), parse = TRUE, 
            col = 'Firebrick', size = 2.5, hjust = 1.2, vjust = -1.2) +
  geom_text(data = anno, aes(x = Inf, y = -Inf, label = sprintf('mu == "%.1f"', round(other, 1))), parse = TRUE, 
            col = 'steelblue1', size = 2.5, hjust = 1.2, vjust = -2.5) +
  scale_fill_manual(values = c('steelblue1', 'firebrick'), guide = FALSE) +
  scale_color_manual(values = c('steelblue1','firebrick'), guide = FALSE) +
  labs(y = 'Burn depth (cm)', x = "") + ylim(c(-30,-7)) + 
  theme(text=element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 5/6)
tag_facet(p, open = '', close = '', tag_pool = c('a','b'))



