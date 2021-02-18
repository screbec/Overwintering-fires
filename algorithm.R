### Apply the overwintering algorithm
## includes the following steps:
## find suitable threshold for distance to antecedent year fire scar
## find suitable threshold for day of emergence (relative to regional snowmelt)
## find suitable lag time for lightning ignition
## find suitable distance from road for anthropogenic ignitions
## apply all thresholds and identity overwintering fires
## plot EG Fig. 3
## plot ED Fig. 4

# load required packages
library(tidyverse)
library(sf)


### read and join all datasets -----------------
# ignition point shapefile
shape = rbind(st_read('ignitionsAK.shp'),
              st_read('ignitionsNWT.shp'))

# distance to previous year fire
nndist = read_csv('nndist.csv', col_types = cols('i', 'n', 'i', 'c', 'n')) %>% select(-X1)

# road intersections
road = read_csv('distroads.csv', col_types = cols('i', 'n', 'n', 'c')) %>% 
 select(-X1) %>% filter(!is.na(distRoad)) %>% distinct()

# lightning intersections
lightning = read_csv('lightning.csv') %>% select(-X1) %>%
  mutate(lightning = ifelse(lightning == '[]', NA, as.integer(lightning))) %>% filter(!is.na(lightning))

# regional yearly ldos ## this is the regional snow melt computed from the snow melt raster
snowMean = read_csv('snow_regional_mean_intAK.csv') %>% add_column(ShortNm = 'AK') %>% 
  bind_rows(read_csv('snow_regional_mean_intNWT.csv') %>% add_column(ShortNm = 'NT')) %>% select(-X1)

# join all datasets
combined = shape %>% 
  left_join(lightning, by = 'ID') %>% 
  left_join(road, by = 'ID') %>% 
  left_join(nndist, by = c("ID", 'Year')) %>%
  left_join(snowMean, by = c('Year', 'ShortNm')) %>%
  mutate(doy_reg_snow = doy - ldos)

### find suitable thresholds for distance and date ------------------------

# split off small overwintering fires
overwintering_small = filter(combined, cause == c('Overwinter small'))
data = filter(combined, !cause == c('Overwinter small'))

# check distributions for suitable thresholds
summary(overwintering_small$nndist)
summary(overwintering_small$doy_reg_snow)

# set holdover thresholds based on distributions
HOdist = 1000 # maximum distance from antecedent fire scar
HOdoysnow = 50.0 # maximum temporal lag between regional snowmeltday and ignition

# what % of overwintering fires lies within threshold?
ecdf(overwintering_small$nndist)(HOdist)
ecdf(overwintering_small$doy_reg_snow)(HOdoysnow)

### find thresholds for lightning lag and road distance -----------------------

# split off all non-potential overwinteirng fires
dataOthers = filter(data, !(nndist < HOdist & doy_reg_snow <= HOdoysnow), !is.na(nndist))

# compare the lightning laga dn road distance for all fires 
# that are labelled as lightning/anthropogenic ignitions
# according to the large fire databases (attribute: cause)
summary(filter(dataOthers, cause == 'Lightning')$lightning) # lightning lag time
summary(filter(dataOthers, cause == 'Human')$distRoad) # distance to road

# set threshold accoridng to distribution
lightThresh = 6 # maximum temporal lag between lightning strike and ignition
roadThresh = 1000 # minimum distance to road in m

# what % of lightning/human cause fires lies within threshold?
ecdf(filter(dataOthers, cause == 'Lightning')$lightning)(lightThresh)
ecdf(filter(dataOthers, cause == 'Human')$distRoad)(roadThresh)


### apply algorithm and write to file -------------

# apply filtering to all other fires
data = mutate(data, ign = ifelse(lightning %in% -1:lightThresh, 'lightning', 
                                  ifelse(distRoad < roadThresh, 'human',
                                  ifelse(nndist < HOdist & doy_reg_snow <= HOdoysnow, 'overwinter', 'unknown'))))
write_csv(combined, 'overwintering_out.csv')

# write out a list of all overwintering IDs (used as identifier in spatial drivers preprocessing)
ids = c(overwintering_small$ID, filter(data, ign == 'overwinter')$ID)
write_csv(data.frame(ids), path = 'overwintering_ids.txt')

### ED Figure 3: histograms for distances and doy since snowmelt---------------------------------

# doy since snowmelt 
plotdata1 = filter(data, !cause == 'Overwinter') %>% mutate(nndistlog = ifelse(nndist == 0, -3, log10(nndist/1000)),
                                                          ign2 = 'Other')
plotdata2 = overwinter_small %>% mutate(nndistlog = ifelse(nndist == 0, -3, log10(nndist/1000)),
                                                                 ign2 = 'Overwinter')
p1 = ggplot() + theme_bw() + labs(y = 'Fraction of data', x = 'Days since regional yearly snowmelt') +
  geom_histogram(data = plotdata1, aes(x = doy_reg_snow, y = ..count../sum(..count..)), 
                 col = 'white', fill = 'cornflowerblue', stat = "bin", binwidth = 10) +
  geom_histogram(data = plotdata2, aes(x = doy_reg_snow, y = ..count../sum(..count..)), 
                 col = 'white', fill = 'firebrick', alpha = 0.5, stat = "bin", binwidth = 10) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "a", fontface = 'bold', hjust = -0.5, vjust = 1.1, size = 4) +
  theme(text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 5/6)

# distance to old firescar (logarithmic x axis is used in plotting)
p2 = ggplot() + theme_bw() + labs(y = '', x = 'Distance to old scar (km)', fill = '') +
  geom_histogram(data = plotdata1, aes(x = nndistlog, y = ..count../sum(..count..), fill = ign2), 
                 col = 'white', stat = "bin", bins = 18) + 
  geom_histogram(data = plotdata2, aes(x = nndistlog, y = ..count../sum(..count..), fill = ign2), 
                 col = 'white', stat = "bin", bins = 18) + 
  annotate(geom = "text", x  = -Inf, y = Inf, label = "b", fontface = 'bold', hjust = -0.5, vjust = 1.1) +
  scale_fill_manual(values = c('cornflowerblue', alpha('firebrick', 0.5))) + 
  scale_x_continuous(breaks = seq(-3, 3, 1), labels = 10^seq(-3, 3, 1), expand = c(0.1, 0)) +
  theme(text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 5/6, 
        legend.position = c(0.75, 0.85), 
        legend.background = element_blank())
grid.arrange(p1, p2, ncol = 2


### ED Figure 4: histogram for lightning and human ignitions -------------------------------------------------------
# lightning
p1 = ggplot() + theme_bw() + labs(y = 'Fraction of data', x = 'Lag Days (Ignition - Lightning)') +
  geom_histogram(data = filter(dataOthers, cause == 'Lightning'), 
                 aes(x = lightning, y = ..count../sum(..count..), fill = ShortNm), 
                 bins = 17, position = "stack", col = 'white') +
  scale_fill_manual(values = c('chartreuse4', 'chartreuse3'), guide = F) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "a", fontface = 'bold', hjust = -0.5, vjust = 1.1, size = 4) +
  geom_vline(xintercept = lightThresh, size = 0.5, col = 'gray22') +
  theme(text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 5/6)

# human
p2 = ggplot() + theme_bw() + labs(y = '', x = 'Distance to road (m)', fill = '') +
  geom_histogram(data = filter(dataOthers, cause == 'Human'), 
                 aes(x = distRoad, y = ..count../sum(..count..), fill = ShortNm), 
                 bins = 25, position = "stack", col = 'white') +
  scale_fill_manual(values = c('chartreuse4', 'chartreuse3')) +
  annotate(geom = "text", x = -Inf, y = Inf, label = "b", fontface = 'bold', hjust = -0.5, vjust = 1.1) +
  geom_vline(xintercept = roadThresh, size = 0.5, col = 'gray22') +
  scale_x_continuous(expand = c(0.1,0), labels = comma) +
  theme(text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 5/6, 
        legend.position = c(0.8, 0.8), 
        legend.background = element_blank())
grid.arrange(p1, p2, ncol = 2)


