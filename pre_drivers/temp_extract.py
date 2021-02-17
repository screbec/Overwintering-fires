# -*- coding: utf-8 -*-
"""
Extract regional temperature statistics (mean, percentile, hot days)
consists of two parts:
    - creating a mask for NARR data based on a rasterised region
    - computing regional averages of statistics
Requirements:
    - daily maximum temperature data (i.e. from NARR)
    - rasterised region (can be different coordinate system and resolution)
Output: csv file with average values for all statistics and regions per year

@author: RCScholten
"""

outfile = 'temp_max_extremes.csv' # output file
climfiles = glob.glob('/clim/*') # folder with yearly files of daily maximum temperature
AKrasterfile = 'akint.tif' # rasterised tif of interior Alaska, side product of burned area normalisation script
NTrasterfile = 'ntint.tif' # rasterised tif of interior NWT, side produict of burned area normalisation script

import gdal
import xarray as xr
import numpy as np
import datetime
import glob
import rtree
from scipy import stats

#%% prepare a mask for the regional mean
## functions for coordinate transformations (NARR data and ecoregion rasters differ in coordinate system
def CreateTransformEPSG(inEPSG, outEPSG):
    '''Creates a coordinate transform between two coord systems 
    specified as EPSG code'''
    # create source csr from EPSG
    src = osr.SpatialReference()                            
    src.ImportFromEPSG(inEPSG)
    # create target csr from epsg
    tgt = osr.SpatialReference()                            
    tgt.ImportFromEPSG(outEPSG)
    
    # create coordinate transform
    transform = osr.CoordinateTransformation(src, tgt)      
    return transform
    
def TransformPoint(transform, x_in, y_in):
    """applies a coordinate transform to an x and y coordinate
    lon = x, lat = y"""
    # create point geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)                      
    point.AddPoint(x_in, y_in)
    point.Transform(transform)
    
    x_out = point.GetX()
    y_out = point.GetY()
    return x_out, y_out
    
def pixel2World(gt, Xpixel, Ypixel):
    '''Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the geospatial coordinate of a pixel location'''
    Xgeo = gt[0] + Xpixel*gt[1] + Ypixel*gt[2]
    Ygeo = gt[3] + Xpixel*gt[4] + Ypixel*gt[5]
    return (Xgeo, Ygeo)


# read x and y files of NARR data and flatten
ds = xr.open_dataset(climfiles)[0], decode_times = False)
orig_shape = ds.lon.shape
lons = np.array(ds.lon).flatten(order = 'C')
lats = np.array(ds.lat).flatten(order = 'C')

# compute all coordinates in CAEA (EPSG 102001)
trans = CreateTransformEPSG(4326, 102001)
coord_tuples = [TransformPoint(trans, float(lon), float(lat)) for lon,lat in zip(lons,lats)]

# read into index with id from flattened df
idx = rtree.index.Index()
for fid, tup in enumerate(coord_tuples):
    idx.insert(fid, (tup[0], tup[1], tup[0], tup[1]), fid)

# read rasterised regions
AKraster = gdal.Open(AKrasterfile)
NTraster = gdal.Open(NTrasterfile)

# now check which x and y are in AKint and NTint
# check for each 255 pixel in mask which NARR center point is closest
for regionRaster in [AKraster, NTraster]:
    # create output list
    nn_narr = []
    # read region raster into array
    regionMask = regionRaster.ReadAsArray()
    # loop over rows and cols
    for Xpx in range(regionRaster.RasterXSize):
        for Ypx in range(regionRaster.RasterYSize):
            # check if pixel is not 0
            if regionMask[Ypx, Xpx] == 255:
                # convert pixel coordinate to CAEA x and y
                x_coord, y_coord = pixel2World(regionRaster.GetGeoTransform(), Xpx, Ypx)
                # find nearest neightbor using index
                nn_narr.append([n.object for n in idx.nearest((x_coord, y_coord, x_coord, y_coord), 1, objects=True)])
    # rempove duplicate indices
    nn_narr = list(set([item for sublist in nn_narr for item in sublist]))
    # create unique masks for AK and NT
    if regionRaster == AKraster:
        narrMask_AK = np.zeros(lons.shape, bool)
        narrMask_AK[nn_narr] = True
    else:
        narrMask_NT = np.zeros(lons.shape, bool)
        narrMask_NT[nn_narr] = True
   
# we also keep the masks in original shape for later
akmask = narrMask_AK.reshape(orig_shape, order = 'C')
ntmask = narrMask_NT.reshape(orig_shape, order = 'C')


#%% this is where the temperature statistics (mean, percentile, hot days) are computed and written to file

# start writing to csv file
f = open(outfile, "w")
f.write('Year,ShortNm,param,temp_perc\n')
# loop through yearly temperature files
for climfile in climfiles:
    ds = xr.open_dataset(climfile, decode_times = False)
    yr = climfile[-7:-3] # extract year form filename

    # clip to may 1 - sept 30
    start = datetime.date(int(yr), 5, 1).timetuple().tm_yday
    end = datetime.date(int(yr), 9, 30).timetuple().tm_yday
    data = ds.air[start:end+1,:,:]

    # calculate percentiles per pixel and reshape to apply masks
    data_perc = np.percentile(data, perc, axis = 0).flatten(order = 'C')   
    
    # compute the statistics
    for stat in ['mean', 'hot_days90', '90']:
        if stat == 'mean':
            # take total mean over area
            AKstat = np.nanmean(np.where(akmask, data, np.nan))
            NTstat = np.nanmean(np.where(ntmask, data, np.nan))  
        elif stat == 'hot_days90':
            # compute days above the long-term 90 percentile (HARDCODED!)
            perc_AK90 = (np.array(data) > 294.2).sum(axis = 0).flatten(order = 'C')
            perc_NT90 = (np.array(data) > 296.6).sum(axis = 0).flatten(order = 'C')
            AKstat = np.mean(perc_AK90[narrMask_AK])
            NTstat = np.mean(perc_NT90[narrMask_NT])
        elif stat == '90':
            # take mean of 90th percentiles over area
            AKstat = np.mean(data_perc[narrMask_AK])
            NTstat = np.mean(data_perc[narrMask_NT])
        # write to file
        outline1 = tuple([yr, 'AK', stat, AKstat])
        outline2 = tuple([yr, 'NT', perc, NT_mean])
        f.write('%i,%s,%s,%f\n' %outline1)
        f.write('%i,%s,%s,%f\n' %outline2)

