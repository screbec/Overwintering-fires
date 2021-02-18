# -*- coding: utf-8 -*-
"""
Computes the first snow-free day of the season
of each pixel using MODSCAG snow fraction Tiffs

Requirements: a continuous time series of all days of spring
of all MODSCAG snow fraction tiles needed for the study area
for each year analysed
i.e. here tiles h10v02 - h13v02 between Julian day 80 - 182, 2000-2018 are used
covering northwestern North America
MODSCAG data should be in original folder structure and with original file names: modscag-historic/year/doy/
Output part 1: Tiles of the first snow-free day of the season (beginning of spring)
Output part 3: csv of regional snowmelt day for each year and region

processing consists of three parts:
- Part 1: computing the first snow-free day (in this script)
- Part 2: merging, reprojecting and clipping tiles to the area of interest (can be done with any software)
- Part 3: computing regional means for each year while excluding NA pixels and burned pixels (in this script)


@author: RCScholten
"""

# import modules
from osgeo import gdal
import numpy as np
import glob
import pandas as pd

# inputs
path_to_modscag = 'path/to/modsag'
path_to_output = 'path/to/output'
outfile = 'snow_regional_mean.csv' # this is the output file which is later read in 04_algorithm.R for the algorithm
doy_file = "doy_Akfed.tif" # this is a tif combining all AKFED day of burning tifs for the study area, one band for each year
snowfile = 'snowfree_merge_reproj_clip.tif' # this is the first snow-free day image, one band for each year (output of part 3, merged, reprojected and clipped)

#%% functions
def first_day(imgfile):
    '''this returns the first snow-free day each year
    does not account for snowfall after this date'''
    if len(imgfile) == 105:
        doy = imgfile[59:62]
    else:
        doy = imgfile[63:66]                        # extract doy from filename
    image_ds = gdal.Open(imgfile)
    arr = image_ds.GetRasterBand(1).ReadAsArray().astype(np.int)
    arr[arr > 100] = 999                        # na pixels get value 999
    arr[(arr > 15) & (arr < 101)] = 888         # snow covered pixels get value 888
    arr[arr < 16] = doy                         # pixels with 15 % or less snow cover get doy
    image_ds = None

    return arr

def common_pixels(img1, img2, res, pxloc=True):
    '''takes two overlapping gdal rasters and checks for common pixels
    returns first common pixel location and width and height'''
    gt_img1 = img1.GetGeoTransform()
    gt_img2 = img2.GetGeoTransform()

    xmin = gt_img1[0] if gt_img1[0] > gt_img2[0] else gt_img2[0]
    ymin = gt_img1[3] if gt_img1[3] < gt_img2[3] else gt_img2[3]
    xmin_px1 = round((gt_img1[0]-xmin)/res)
    xmin_px2 = round((gt_img2[0]-xmin)/res)
    ymin_px1 = round((gt_img1[3]-ymin)/res)
    ymin_px2 = round((gt_img2[3]-ymin)/res)

    xpx_img1 = img1.RasterXSize - xmin_px1
    xpx_img2 = img2.RasterXSize - xmin_px1
    ypx_img1 = img1.RasterYSize - ymin_px1
    ypx_img2 = img2.RasterYSize - ymin_px2

    Xsize = xpx_img1 if xpx_img1 < xpx_img2 else xpx_img2
    Ysize = ypx_img1 if ypx_img1 < ypx_img2 else ypx_img2

    if pxloc:
        return xmin_px1, xmin_px2, ymin_px1, ymin_px2, Xsize, Ysize
    else:
        return xmin, ymin, Xsize, Ysize

#%% Optional: check which tiles are missing between 2000 and 2017
for tile in ['h10v02', 'h11v02', 'h12v02', 'h13v02']:
    for y in range(2000, 2018):
        imgs = glob.glob(path_to_modscag + '/modscag-historic/' + str(y) + '/*/*%s.006*' %tile)
        print(len(imgs))
        if not (len(imgs) > 101) and not (y in [2001, 2002]):
            print(tile, y)
            doys = [int(imgfile[69:72]) for imgfile in imgs]
            missing = [doy for doy in list(range(80, 183)) if doy not in doys]
            print(missing)



#%% Part 1 - this only creates the first snow-free day tiles
for tile in ['h10v02', 'h11v02', 'h12v02', 'h13v02']:       # for each tile

    # read example dataset
    in_ds1 = gdal.Open(glob.glob(path_to_modscag + '/modscag-historic/2009/*/*%s.006*' %tile)[0])

    # create output datasource for each tile
    drv = gdal.GetDriverByName('GTiff')
    outds = drv.Create(path_to_output + '/%s.tif' %(tile),
                       in_ds1.RasterXSize, in_ds1.RasterYSize, 19, gdal.GDT_Byte)
    outds.SetGeoTransform(in_ds1.GetGeoTransform())     # sets geotransform from input
    outds.SetProjection(in_ds1.GetProjection())         # sets projection from input
    in_ds1 = None

    for y in range(2000, 2019):
        # set folder structure to data
        files = path_to_modscag + '/modscag-historic/' + str(y) + '/*/*%s.006*'

        # gather filenames, check if time series is complete (less than 3 days missing)
        imgs = glob.glob(files %tile)
        if not (len(imgs) > 101) and not (y in [2001, 2002]):
            print(y, 'could not be processed. Incomplete time series.')
            continue

        # compute first snow-free day
        doy_arr = [first_day(file) for file in imgs]
        doy_arr = np.asarray(doy_arr)
        flat_arr = np.amin(np.ma.masked_greater(doy_arr, 800), axis=0)
        na_val = 999

        # for each pixel check the preceding 4 days: if NA, then set doy to 0
        for row in range(flat_arr.shape[0]):
            for col in range(flat_arr.shape[1]):
                if flat_arr[row, col]:
                    doy_idx = flat_arr[row, col]-80 # doy idx = doy - 80 (since 80 is starting point)
                    # print(flat_arr[row,col], doy_arr[doy_idx, row, col]) # check if index is correct
                    days_before = doy_arr[doy_idx-4:doy_idx, row, col]
                    # if any of those days is not na (has data), the date can be used
                    if np.all(days_before == na_val):
                        flat_arr[row, col] = 0

        # replace NAs (masked) with 0
        flat_arr = flat_arr.filled(fill_value=0)

        # write out new array to rasterband (band for each year)
        outband = outds.GetRasterBand(y-1999)
        outband.WriteArray(flat_arr)
        outband.SetNoDataValue(0)
        outds.FlushCache()                          # saves to disk
        outband = None

    outds = None

#%% Part 2


# next steps:
# - Merge and reproject all tiles to required coordinate system (i.e. in QGIS)
# - Clip to study area (i.e. interior ecoregions)



#%% Part 3
#!!! after merging, reprojecting and clipping:
# take regional average while excluding:
# - pixels with missing data for any of the years and
# - pixels which have burned within 2000-2018
# this has to be done for all regions that are studied separately

# read image files and define output

ba_img = gdal.Open(doy_file)
snow_image = gdal.Open(snowfile)

# first compare the geotransforms of both images to pick overlapping pixels
xmin_ba, xmin_snow, ymin_ba, ymin_snow, xn, yn = common_pixels(ba_img, snow_image, 500)

# build a mask to exclude all previous burned area for each year (only include value 0 = unburned)
for bnd in range(ba_img.RasterCount)[:-1]:
    band = ba_img.GetRasterBand(bnd + 1)
    arr1 = band.ReadAsArray(xmin_ba, ymin_ba, xn, yn)
    try:
        mask = mask & (arr1 == 0)
    except:
        mask = (arr1 == 0)
ba_img = None

# extend mask by excluding all NA (NA = 0) pixels of snow images
for bnd in range(snow_image.RasterCount):
    band = snow_image.GetRasterBand(bnd + 1)
    arr1 = band.ReadAsArray(xmin_snow, ymin_snow, xn, yn)
    mask = mask & (arr1 != 0)

# then loop through images, apply the mask and return the np.nanmean
res1 = []
for bnd in range(snow_image.RasterCount):
    out = [bnd + 2000]
    band = snow_image.GetRasterBand(bnd + 1)
    arr1 = band.ReadAsArray(xmin_snow, ymin_snow, xn, yn)
    out.append(np.nanmean(arr1[mask]))
    res1.append(out)
snow_image = None

df = pd.DataFrame(res1, columns=['Year', 'ldos'])
df.to_csv(outfile)
