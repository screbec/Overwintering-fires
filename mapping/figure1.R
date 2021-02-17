# Code for Figure 1: Landsat time series of an overwintering fire in Alaska

### load required packages
library(tidyverse)
library(raster)
library(sf)
library(ggspatial)

# inputs
## for base images
# perimeters of the two fire scars (original scar, new burn scar the year after)
overwinter_perim = 'overwinter_perimeter.shp'
# landsat 8 imagery
imgpaths = c('LC08_072016_20150924_3913.tif',
             'LC08_072016_20160412_3913.tif',
             'LC08_072016_20160530_3913.tif', 
             'LC08_072016_20160615_3913.tif')
## for map inset (Fig. 1b):
# shapefile of country boundaries of Alaska and Canada, i.e. of Alaska, Yukon and Northwest Territories
AKCA = 'AKYTNT.shp'
# shapefile of interior Alaska
akint = 'ak_interior.shp'

### load data --------------------------------------------------------------
### load shape for inset map and transform to UTM
AK = st_read(AKCA) %>% st_transform(crs = 32605)
innerAK = st_read(AKint) %>% st_transform(crs = 32605) %>%
  mutate(NA_L2NAME = 'Interior Alaska')

### create clip features for image
clip_max = as(extent(429433, 443000, 7004200, 7020000), 'SpatialPolygons')
img = stack(imgpaths[1])
crs(clip_max) = crs(img)
clip_max = st_as_sf(clip_max)


### load shape for overlay and clip
perims = st_read(overwinter_perim) %>%
  st_transform(crs = 32605) %>% st_intersection(clip_max)

### Read and crop images
stagecnt = 1
for (stage in imgpaths){
  # read and crop image
  img = stack(stage) %>% crop(clip_max) %>% as.data.frame(xy = TRUE)
  colnames(img) = make.names(c('x', 'y', c(1:(dim(img)[2]-2))))
  
  # assign to new variable name
  varname = paste('stage', stagecnt, sep = '_')
  img = add_column(img, stage = stagecnt)
  assign(varname, img)
  stagecnt = stagecnt + 1
}
rm(stagecnt, img, stage, varname)


### bind to tidy dataframe including all stages
data = dplyr::select(stage_1, x, y, X3, X4, X5, stage) %>%
  rbind(dplyr::select(stage_2, x, y, X3, X4, X5, stage)) %>%
  rbind(dplyr::select(stage_3, x, y, X3, X4, X5, stage)) %>%
  rbind(dplyr::select(stage_4, x, y, X3, X4, X5, stage)) %>%
  rename(b = X3, g = X4, r = X5)


# function for contrast stretch
stretch <- function(x, minval, maxval){
  x = ifelse(x <= minval, 0, ifelse(x > maxval , 1, x/(maxval-minval)-(minval/(maxval-minval))))
}
quantiles = data %>% group_by(stage) %>% summarise(r_quantmin = quantile(r, .02), r_quantmax = quantile(r, .98),
                                                   g_quantmin = quantile(g, .02), g_quantmax = quantile(g, .98),
                                                   b_quantmin = quantile(b, .02), b_quantmax = quantile(b, .98))

# based on the quantiles, select thresholds for min and max:
minmax = cbind(rep(c(1,2,3,4)), 
               c(5000,5000,5000,5000),
               c(3000,7000,3000,3000),
               c(3000,7000,3000,3000),
               c(15000,22000,22000,22000),
               c(20000,22000,40000,40000),
               c(15000,22000,30000,30000)) %>% as.data.frame()
colnames(minmax) = c('stage', 'rmin', 'gmin', 'bmin', 'rmax', 'gmax', 'bmax')

# apply contrast stretch
dat_stretched = left_join(data, minmax, by = 'stage') %>%
  mutate(r = stretch(r, rmin, rmax), g = stretch(g, gmin, gmax), b = stretch(b, bmin, bmax))

### mapping ------------------------------------------------------------------------

plotmap = function(plotdat, perim, perimfill, perimline, perimalpha, datelabel, anno = F, perim2 = F, inset = F){
  text_size_pt = 7 # used for theme
  text_size_mm = text_size_pt/2.8 #used for geom_text and annotate text
  
  p = ggplot() + theme_void() + theme(legend.position = "none") + 
    geom_point(data = plotdat, aes(x = x, y = y, col = rgb(r, g, b))) +
    scale_color_identity() +
    coord_fixed(ratio = 1) + theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "pt")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_sf(data = filter(perims, FIREYEAR == 2015), col = perimline, fill = perimfill, alpha = perimalpha) +
    annotate(geom="label", x = -Inf, y = -Inf, label = datelabel, 
             col = NA, fill = 'white', alpha = 0.5, size = text_size_mm, vjust = -0.2, hjust = -0.08) +
    annotate(geom="text", x = -Inf, y = -Inf, label = datelabel, 
             col = 'black', size = text_size_mm, vjust = -1.2, hjust = -0.2) + 
    theme(text = element_text(size = text_size_pt))
  if (anno){
    p = p + annotation_scale(style = 'ticks', location = 'br', 
                             text_col = 'gray', line_col = 'gray')
  }
  if (perim2){
    p = p + geom_sf(data = filter(perims, FIREYEAR == 2016), col = 'cornflowerblue', fill = NA)
  }
  if (inset){
    # plot an inset map
    pI  = ggplot() + theme_bw() + labs(x = '', y = '') +
      geom_sf(data = AK, fill = "ivory", show.legend = F, color = "gray75", lwd = 0.4) +
      geom_sf(data = innerAK, fill = "gray75", color = NA, show.legend = F, alpha = 0.3) +
      coord_sf(label_axes = "--SWEN", xlim = c(-700000, 1046578), ylim = c(6235000, 7936072)) +
      theme_void() +
      theme(panel.background=element_rect(fill = "white", colour = NA)) +
      geom_point(aes(x = 439716, y = 7011250), colour = "firebrick", size = 2)
    p = p + annotation_custom(gtable::gtable_filter(ggplotGrob(pI), "panel"), 
                                   xmin = 438443, xmax = Inf, ymin = 7015643, ymax = Inf)
  }
  print(p)
}
# maps of separate images
plotmap(plotdat = filter(dat_stretched, stage == 1), 
        perimline = 'gray60', perimfill = 'gray60', perimalpha = 0.05,
        datelabel = '24.09.2015')
plotmap(plotdat = filter(dat_stretched, stage == 2), 
        perimline = 'white', perimfill = 'white', perimalpha = 0.05,
        datelabel = '12.04.2016', inset = T)
plotmap(plotdat = filter(dat_stretched, stage == 3), 
        perimline = 'gray60', perimfill = 'gray60', perimalpha = 0.1,
        datelabel = '30.05.2016')
plotmap(plotdat = filter(dat_stretched, stage == 4), 
        perimline = NA, perimfill = 'white', perimalpha = 0.2,
        datelabel = '15.06.2016', anno = T, perim2 = T)

