# -*- coding: utf-8 -*-
"""
Normalise burned area from large fire databases using AKFED burned area
while accounting for differences in the minimum mapping unit (MMU)
Steps:
    - create a mask for interior ecoregions of Alaska and Northwest Territories
    - compute ratios of annual burned area from AKFED/annual burned area from fire perimeters for each region
    - for different versions of ratios depending on MMU, fire perimeters below a threshold are removed
    - predict AKFED burned area equivalent from historic fire perimeters using the appropriate ratios
    - write annual burned area to file for both regions
Requirements:
    - AKFED day of burning (or other burned area dataset which considers unburned islands)
    - fire perimeters (i.e. from national fire databases) for 1975-2018
Output: a csv file with annual burned area for both regions

@author: RCScholten
"""
# raster of AKFED burned area
akfed_AK_2001 = 'AKFED_20012018_intAK.tif' # AKFED day of burning for interior Alaska
akfed_NT_2001 = 'AKFED_20012018_intNT.tif' # AKFED day of burning for interior Northwest Territories
perim_2001 = 'burnedarea_20012018' # rasterized fire perimeters 2001-2018 with unique ids (i.e. from official fire databases, here FireIDs > 45000 represent NWT & IDs < 45000 AK)
perim_1975 = 'burnedarea_19752000' # rasterized fire perimeters 1975-2000
outfile = 'burned_area_trend.csv' # output file with yearly burned area for interior AK and interior NWT

# import modules
import sys
from osgeo import gdal, ogr
import numpy as np
import pandas as pd


def readDS(filename, drvstring = 'ESRI Shapefile'):
    try:
        driver = ogr.GetDriverByName(drvstring)
    except:
        print('Driver not found.')
        sys.exit(1)
    ds = driver.Open(filename)
    lyr = ds.GetLayer()
    return(ds, lyr)

def layer2raster(inlyr, res, outname, extent = None):
    '''takes a shapefile layer and transforms it to a raster with a 
    specified resolution and output name 
    optionally the extent of the new raster can be specified
    in the form of (minx, maxx, miny, maxy)'''
    if not extent:
        # get ulx and uly from extent
        minx, maxx, miny, maxy = inlyr.GetExtent() 
    else:
        # take specified extent
        minx, maxx, miny, maxy = extent 
    # calculate number of x and y pixels from extent and resolution
    xpx = int(np.ceil((maxx - minx)/res))
    ypx = int(np.ceil((maxy - miny)/res))
    # create the geotransform
    geotrans = (np.floor(minx), res, 0, np.ceil(maxy), 0, -res)  
    # create output datasource
    driver = gdal.GetDriverByName('GTiff')
    outraster = driver.Create(outname, xpx, ypx, 1, gdal.GDT_Byte)
    outraster.SetGeoTransform(geotrans)
    outraster.SetProjection(inlyr.GetSpatialRef().ExportToWkt())
    gdal.RasterizeLayer(outraster, [1], inlyr, None)
    return outraster

#%% prepare a mask for the regional mean

# get extent and resolution from the raster
in_ds1 = gdal.Open(perim_2001)
geoTransform = in_ds1.GetGeoTransform()
minx = geoTransform[0]
maxy = geoTransform[3]
maxx = minx + geoTransform[1] * in_ds1.RasterXSize
miny = maxy + geoTransform[5] * in_ds1.RasterYSize
ext = (minx, maxx, miny, maxy)
res = geoTransform[1]

# set output path and name for rasterised shapefiles
outfileAK = 'akint.tif'
outfileNT = 'ntint.tif'

# read shapefiles of interior ecoregions
AKint_ds, AKint_lyr = readDS('ak_interior.shp')
NTint_ds, NTint_lyr = readDS('nwt_interior.shp')

# rasterize shapefiles
AKraster = layer2raster(AKint_lyr, res, outfileAK, ext)
NTraster = layer2raster(NTint_lyr, res, outfileNT, ext)

# read rasters as array
akmask = AKraster.ReadAsArray()
ntmask = NTraster.ReadAsArray()    

# close datasets
NTraster = AKraster = in_ds1 = None

  
#%% for overlapping period of burned area and perimeters: compute the ratio of burned area

# open perimeter dataset for 2001-2018
in_ds1 = gdal.Open(perim_2001)
arr1 = in_ds1.ReadAsArray()

# now  we have to include the mmu (minimum mapping unit)
# this is 200 ha for NWT (8 AKFED pixels)
# for AK we have different mmus depending on the year,
# 405 ha (16 pixels) for 75-86, 40.5 ha (2px) for 87-15
# fire scars < 25ha are less then one pixel

# count occurrences of each fire id per year
# has to be done per year since fire ids are not unique...
# ids > 45000 are in NWT
idlist2 = []
idlist8 = []
idlist16 = []
for yr in range(arr1.shape[0]):
    unique, counts = np.unique(arr1[yr,:,:], return_counts=True) 
    idlist2.append(unique[(counts < 3) & (unique < 45000)]) # these ids have to be ignored in AK, 1987-2015
    idlist8.append(unique[(counts < 9) & (unique > 45000)]) # these ids have to be ignored in NWT
    idlist16.append(unique[(counts < 17) & (unique < 45000)]) # these ids have to be ignored in AK, 1975-1986
idlist2 = [item for sublist in idlist2 for item in sublist]
idlist16 = [item for sublist in idlist16 for item in sublist]
idlist8 = [item for sublist in idlist8 for item in sublist]

# broadcast mask to the shape defined by number of years
ak3d_mask = np.broadcast_to(akmask  == 0, arr1.shape)
nt3d_mask = np.broadcast_to(ntmask  == 0, arr1.shape)

# apply mask for AK
ak1 = np.copy(arr1)
ak1[ak3d_mask] = 0
# create datasets with 2px and 16px removed 
ak2 = np.copy(ak1)
ak2[np.in1d(ak2, idlist2).reshape(ak2.shape)] = 0
ak16 = np.copy(ak1)
ak16[np.in1d(ak16, idlist16).reshape(ak2.shape)] = 0
# convert to binary
ak2[ak2 > 0] = 1
ak16[ak16 > 0] = 1

# apply mask for NT
nt1 = np.copy(arr1)
nt1[nt3d_mask] = 0
# for nt8 remove everything below 8px per year
nt8 = np.copy(arr1)
nt8[np.in1d(nt8, idlist2).reshape(nt8.shape)] = 0
# convert to binary
nt8[nt8 > 0] = 1

# count number of ones * area = total 'burned area' from fire perimeters
aksum2 = np.sum(ak2, (1,2)) * 0.25
aksum16 = np.sum(ak16, (1,2)) * 0.25
ntsum8 = np.sum(nt8, (1,2)) * 0.25

## read AKFED burned area data for same time frame
image_ds = gdal.Open(akfed_AK_2001,)
ba = image_ds.ReadAsArray() 
ba[ba > 0.0] = 1 # set all burn day values to 1 to take sume
ak_ba = np.sum(ba, (1,2)) * 0.25 # total burned area from AKFED (km2)

image_ds = gdal.Open(akfed_NT_2001,)
ba = image_ds.ReadAsArray() 
ba[ba > 0.0] = 1
nt_ba = np.sum(ba, (1,2)) * 0.25 # total burned area from AKFED (km2)


# calculate ratios for all versions
akratio2 = np.sum(ak_ba)/np.sum(aksum2) # Alaska, using only perimeters > 40.5 ha (1987-2015)
akratio16 = np.sum(ak_ba)/np.sum(aksum16) # Alaska, using only perimeters > 405ha (1975-1986)
ntratio8 = np.sum(nt_ba)/np.sum(ntsum8) # NTW, using only perimeters > 200 ha


#%% compute yearly sum of archive dataset and normalise using ratios  
# read archive dataset (1975 - 2000)
in_ds = gdal.Open(ba_1975)
aksum = []
ntsum = []
# loop through bands
for band in range(in_ds.RasterCount):
    arr = in_ds.GetRasterBand(band).ReadAsArray()

    # apply mask and convert to binary
    ak = np.copy(arr)
    ak[akmask == 0] = 0
    ak[ak > 0] = 1
    
    nt = np.copy(arr)
    nt[ntmask == 0] = 0
    nt[nt > 0] = 1
    
    # count number of ones
    aksum.append(np.sum(ak))
    ntsum.append(np.sum(nt))

aksum = np.array(aksum)
ntsum = np.array(ntsum)

# apply ratios to historical dataset
# for Alaska the ratios are applied depending on the time frame
ak_ba_pred = np.concatenate((aksum[:12] * akratio16 * 0.25, aksum[12:] * akratio2 * 0.25))
# for northwest territories one ratio is applied to the full dataset
nt_ba_pred = ntsum * ntratio8 * 0.25


# concatenate with new dataset and write to file
xax = np.arange(start = 1975, stop = 2019)
aktotal = np.concatenate((ak_ba_pred, ak_ba), axis=None)
nttotal = np.concatenate((nt_ba_pred, nt_ba), axis=None)
df = pd.DataFrame({'NT' : nttotal, 'AK' : aktotal, 'Year' : xax})
df.to_csv(outfile)

