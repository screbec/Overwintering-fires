# Figure 4: plot overwintering ignitions on maps of elevation i.e.

### load required packages
library(tidyverse)
library(raster)
library(sf)
library(viridis)
library(ggspatial)


# ignition shapefile of overwintering fires
shape = st_read('overwinter_small.shp') %>% mutate(ign = 'overwinter small') %>%
  rbind(st_read('overwinter.shp') %>% filter(ign == 'overwinter'))

# imagery (DEM)
img = raster('arcticdem.tif')

# AK/CA shapefile for overview map
AKCA = st_read('AKCA.shp')

# path to lakes shapefile for NT
lakepath = 'lakes.shp'


### functions ----

# function creating a sp dataframe from the (buffered) extent of another shapefile
bbox2clip = function(shape, buffer = 50000){
  clip = as(extent(shape %>% st_buffer(buffer)), 'SpatialPolygons')
  crs(clip) = crs(shape)
  clip = st_as_sf(clip) # this is in shapefile projection (for the map)
  return(clip)
}

# function for plotting raster maps with sp dataframe and coordinates
plotimg = function(rasterimg, ignshape, regionshape, region, lakepath, cbminmax){
  ## prepare data fpr plotting ##
  # filter shapefile for overwintering fires and region
  shape_OV = filter(ignshape, ShortNm == region) %>% st_centroid()
  # create clip geometry using bounding box of the shapefile
  clip = bbox2clip(shape_OV)
  # reproject and clip region shapefile
  REG = regionshape %>% st_transform(crs(ignshape)) %>% st_intersection(clip)
  
  ## start plotting ##
  
  # create background image with coordinates
  pcoords = ggplot() + theme_bw() +
    theme(legend.position = "none", text = element_text(size = 8)) +
    geom_sf(data = shape_OV, col = 'mediumblue', fill = 'mediumblue') +
    geom_sf(data = clip, col = 'black', fill = 'NA') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    annotation_scale(style = 'ticks', location = 'bl', height = unit(0.2, "cm"), text_cex = 0.6,
                     text_col = 'black', line_col = 'black')
  print(pcoords) # this can be exported as svg
  
  # map raster
  pmap = ggplot() + theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(0.5,0.5,0.5,0.5), "pt")) +
    geom_raster(data = rasterimg, aes(x = x, y = y, fill = val))
  # for NT map add the lakes as orientation
  if (region == 'NT') {
    lakes = st_read(lakepath) %>% 
      filter(name_en %in% c('Great Bear Lake', 'Great Slave Lake')) %>% 
      st_transform(crs(shape_OV))
    pmap = pmap +
      geom_sf(data = lakes, col = '#91c8cd', fill = '#b4e3f8', lwd = 0.1)
  } 
  pmap = pmap + 
    geom_sf(data = REG, fill = NA, show.legend = F, color = "black", lwd = 0.2) +
    geom_sf(data = clip, col = 'white', fill = 'NA') +
    geom_sf(data = shape_OV, aes(size = ign), fill = 'steelblue1', col = 'black', shape = 21) +
    scale_size_manual(values = c(3, 2)) +
    scale_fill_viridis(option = "B", na.value = NA, direction = -1, limits = cbminmax) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  print(pmap) # this can be exported as tiff
}

### preprocessing for plotting ------------------------------------

# simplify geometries for overview map (pckg rmapshaper)
AKCA_simple = st_simplify(AKCA, dTolerance = 5000)

# create boundary boxes
clip_AK = filter(shape, ShortNm == 'AK') %>% st_centroid() %>% 
  bbox2clip() %>% st_transform(crs(img))
clip_NT = filter(shape, ShortNm == 'NT') %>% st_centroid() %>% 
  bbox2clip() %>% st_transform(crs(img))

# clip for AK and NT and combine for colourbar min and max
demNT = img %>% crop(clip_NT) %>% 
  projectRaster(crs = crs(shape)) %>%
  as.data.frame(xy = TRUE) %>% rename(val = `arcticdem_mosaic_1km_v3.0`) 
demAK = img %>% crop(clip_AK) %>% 
  projectRaster(crs = crs(shape)) %>%
  as.data.frame(xy = TRUE) %>% rename(val = `arcticdem_mosaic_1km_v3.0`)

# set limits for colorbar
colbarmin = min(c(min(na.omit(demNT$val)), min(na.omit(demAK$val))))
colbarmax = max(c(max(na.omit(demNT$val)), max(na.omit(demAK$val))))

### plotting -----------------------------------------

# plot overview image
p_overview = ggplot() + theme_bw() + labs(x = '', y = '') +
  theme(text = element_text(size = 8)) + 
  geom_sf(data = AKCA_simple, fill = "ivory", show.legend = F, color = "gray75", lwd = 0.2) +
  geom_sf(data = clip_AK, col = 'Firebrick', fill = NA) +
  geom_sf(data = clip_NT, col = 'cornflowerblue', fill = NA) +
  coord_sf(crs = st_crs(102001), 
           xlim = c(-3200000, 3000000), ylim = c(1500000, 5000000))
print(p_overview)

# color bar only
p_cb = ggplot() + theme_void() + labs(fill = 'Altitude (m)') +
  geom_raster(data = sample_n(demAK, 1000) , aes(x = x, y = y, fill = val)) +
  scale_fill_viridis(option = "B", na.value = NA, direction = -1, limits = c(colbarmin, colbarmax))
print(p_cb)

# plot maps
plotimg(demAK, shape, filter(AKCA, ShortNm == 'AK'), 'AK', lakepath, 
        c(colbarmin,colbarmax))
plotimg(demNT, shape, filter(AKCA, ShortNm %in% c('NWT', 'NVT')), 'NT', lakepath,
        c(colbarmin,colbarmax))
