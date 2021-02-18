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
from osgeo import gdal, osr, ogr
import xarray as xr
import numpy as np
import datetime
import glob
import rtree

# inputs
outfile = 'temp_max_extremes.csv' # output file
climfiles = glob.glob('/clim/*') # folder with yearly files of daily maximum temperature
ak_rasterfile = 'akint.tif' # rasterised tif of interior Alaska, side product of burned area normalisation script
nt_rasterfile = 'ntint.tif' # rasterised tif of interior NWT, side produict of burned area normalisation script

#%% prepare a mask for the regional mean
## functions for coordinate transformations (NARR data and ecoregion rasters differ in coordinate system
def create_transform_epsg(in_epsg, out_epsg):
    '''Creates a coordinate transform between two coord systems
    specified as EPSG code'''
    # create source csr from EPSG
    src = osr.SpatialReference()
    src.ImportFromEPSG(in_epsg)
    # create target csr from epsg
    tgt = osr.SpatialReference()
    tgt.ImportFromEPSG(out_epsg)

    # create coordinate transform
    transform = osr.CoordinateTransformation(src, tgt)
    return transform

def transform_point(transform, x_in, y_in):
    """applies a coordinate transform to an x and y coordinate
    lon = x, lat = y"""
    # create point geometry from coordinates
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(x_in, y_in)
    point.Transform(transform)

    x_out = point.GetX()
    y_out = point.GetY()
    return x_out, y_out

def pixel_to_world(gt, xpixel, ypixel):
    '''Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate the geospatial coordinate of a pixel location'''
    xgeo = gt[0] + xpixel*gt[1] + ypixel*gt[2]
    ygeo = gt[3] + xpixel*gt[4] + ypixel*gt[5]
    return (xgeo, ygeo)


# read x and y files of NARR data and flatten
ds = xr.open_dataset(climfiles[0], decode_times=False)
orig_shape = ds.lon.shape
lons = np.array(ds.lon).flatten(order='C')
lats = np.array(ds.lat).flatten(order='C')

# compute all coordinates in CAEA (EPSG 102001)
trans = create_transform_epsg(4326, 102001)
coord_tuples = [transform_point(trans, float(lon), float(lat)) for lon, lat in zip(lons, lats)]

# read into index with id from flattened df
idx = rtree.index.Index()
for fid, tup in enumerate(coord_tuples):
    idx.insert(fid, (tup[0], tup[1], tup[0], tup[1]), fid)

# read rasterised regions
ak_raster = gdal.Open(ak_rasterfile)
nt_raster = gdal.Open(nt_rasterfile)

# now check which x and y are in AKint and NTint
# check for each 255 pixel in mask which NARR center point is closest
for region_raster in [ak_raster, nt_raster]:
    # create output list
    nn_narr = []
    # read region raster into array
    region_mask = region_raster.ReadAsArray()
    # loop over rows and cols
    for x_px in range(region_raster.RasterXSize):
        for y_px in range(region_raster.RasterYSize):
            # check if pixel is not 0
            if region_mask[y_px, x_px] == 255:
                # convert pixel coordinate to CAEA x and y
                x_coord, y_coord = pixel_to_world(region_raster.GetGeoTransform(), x_px, y_px)
                # find nearest neightbor using index
                nn_narr.append([n.object for n in idx.nearest((x_coord, y_coord, x_coord, y_coord), 1, objects=True)])
    # rempove duplicate indices
    nn_narr = list(set([item for sublist in nn_narr for item in sublist]))
    # create unique masks for AK and NT
    if region_raster == ak_raster:
        narr_mask_ak = np.zeros(lons.shape, bool)
        narr_mask_ak[nn_narr] = True
    else:
        narr_mask_nt = np.zeros(lons.shape, bool)
        narr_mask_nt[nn_narr] = True

# we also keep the masks in original shape for later
ak_mask = narr_mask_ak.reshape(orig_shape, order='C')
nt_mask = narr_mask_nt.reshape(orig_shape, order='C')


#%% this is where the temperature statistics (mean, percentile, hot days) are computed and written to file

# start writing to csv file
f = open(outfile, "w")
f.write('Year,ShortNm,param,temp_perc\n')
# loop through yearly temperature files
for climfile in climfiles:
    ds = xr.open_dataset(climfile, decode_times=False)
    yr = climfile[-7:-3] # extract year form filename

    # clip to may 1 - sept 30
    start = datetime.date(int(yr), 5, 1).timetuple().tm_yday
    end = datetime.date(int(yr), 9, 30).timetuple().tm_yday
    data = ds.air[start:end+1, :, :]

    # calculate 90th percentile per pixel and reshape to apply masks
    data_perc = np.percentile(data, 90, axis=0).flatten(order='C')

    # compute the statistics
    for stat in ['mean', 'hot_days90', '90']:
        if stat == 'mean':
            # take total mean over area
            ak_stat = np.nanmean(np.where(ak_mask, data, np.nan))
            nt_stat = np.nanmean(np.where(nt_mask, data, np.nan))
        elif stat == 'hot_days90':
            # compute days above the long-term 90 percentile (HARDCODED!)
            perc_AK90 = (np.array(data) > 294.2).sum(axis=0).flatten(order='C')
            perc_NT90 = (np.array(data) > 296.6).sum(axis=0).flatten(order='C')
            ak_stat = np.mean(perc_AK90[narr_mask_ak])
            nt_stat = np.mean(perc_NT90[narr_mask_nt])
        elif stat == '90':
            # take mean of 90th percentiles over area
            ak_stat = np.mean(data_perc[narr_mask_ak])
            nt_stat = np.mean(data_perc[narr_mask_nt])
        # write to file
        outline1 = tuple([yr, 'AK', stat, ak_stat])
        outline2 = tuple([yr, 'NT', 90, nt_stat])
        f.write('%i,%s,%s,%f\n' %outline1)
        f.write('%i,%s,%s,%f\n' %outline2)
