# -*- coding: utf-8 -*-
"""
Clip burn depth and environmental variables
(tree cover, tree species, elevation, slope, topsoil carbon content)
to fire perimeters and assign class overwinter/other to each perimeter
based on the id

Requirements:
    - list of overwinteirng fire ids (output of algorithm.R)
    - and parent fires (extracted from nndist)
    - fire perimeters (i.e. rasterised AKFED burned area with unique IDs for each perimeter)
    - AKFED burn depth data
    - tree cover data (i.e. from MODIS continuous fields)
    - other geographic data (i.e. here tree species, elevation, slope, soil carbon)

Output: csv file of mean values of each geographic variable for each fire perimeter

@author: RCScholten
"""

# import required modules
import rasterio
import fiona
import numpy as np
import pandas as pd

# inputs (partly set below...)
outfile = 'spatial.csv'
overwinter_ids = 'overwintering_ids.txt' # list of overwintering ids (output from algorithm)
nndist = 'nndist.csv'
depth_path = 'depth_20012018.tif' # burn depth raster for clipping
tc_path = 'tc_500_20002017.tif' # take the tree cover of the year before the fire

#%% main
# load overwintering ids
with open(overwinter_ids, 'r') as f:
    ho_ids = f.read().splitlines()
ho_ids = [int(hoid) for hoid in ho_ids[1:]]

# extract list of their original fire perimeter ID
id2fireid = pd.read_csv(nndist)
fire_ids = list(set(id2fireid[id2fireid['ID'].isin(ho_ids)].FireID_tgt))

# start writing csv file
f = open(outfile, "w")
f.write('FireID,Year,ShortNm,ign2,depth_mean,n_depth,tc,n_tc,sp_bs,sp_ws,sp_dc,sp_gs,sp_nv,sp_p,n_ts,Elev,n,slope,n,soil30,soil100,soil200,soil300,n\n')

# loop through regions
for region in ['NT', 'AK']:

    # inputs and outputs according to region (these are specific to each region)
    if region == 'NT':
        shape_path = 'perimeters_NT.shp'
        img_paths = ['treespecies_nwt.tif',
                     'arcticdem_NWT.tif',
                     'arcticdem_slope_NWT.tif',
                     'ncscd_ytnt.tif']

    elif region == 'AK':
        shape_path = 'perimeters_AK.shp'
        img_paths = ['treespecies_ak.tif',
                     'arcticdem_AK.tif',
                     'arcticdem_slope_AK.tif',
                     'ncscd_ak.tif']

    # start clipping
    # open images with rasterio
    depth_img = rasterio.open(depth_path)
    tc_img = rasterio.open(tc_path)

    # iteration through features
    with fiona.open(shape_path, "r") as shapefile1:
        for feat1 in shapefile1:
            # extract relevant information from feature
            geom = feat1["geometry"]
            obj_id = feat1["properties"]["FireID"]
            year = int(feat1["properties"]["Year"])
            # determine whether the feature produced an overwintering fire or not
            ign2 = 'overwinter' if obj_id in fire_ids else 'other'

            # year 2018 has to be skipped since we don't have burned area for 2019
            # and consequently don't know if fire scars in 2018 caused holdovers in 2019
            if year > 2017:
                continue

            # extract tiff values masking with feature
            arr_depth, trans = rasterio.mask.mask(depth_img, [geom], crop=True) #, all_touched=True)
            arr_tc, trans = rasterio.mask.mask(tc_img, [geom], crop=True, filled=False) #, all_touched=True)

            # convert NoData to NAN
            arr_depth = np.where(arr_depth == 0.0, np.nan, arr_depth)

            # calculate stats on target year
            mean_depth = np.nanmean(arr_depth[year-2001, :, :])     # 2001 = year0
            n_depth = np.count_nonzero(~np.isnan(arr_depth[year-2001, :, :]))
            # rasterio assigns 0 as nodataval if the nodatavalue is missing in a raster
            # since 0 is a valid value for our rasters we need a workaround
            # we return the masked array (l 108, filled = False)
            # and count the number of 'False' = not masked values in the mask for n
            # if the whole array is masked (geom does not hit the center point of a raster pixel)
            # it cannot calculate the nanmean, therefore the if-clause
            if False in arr_tc.mask:
                mean_tc = np.nanmean(arr_tc[year-2001, :, :])           # 2000 = year0
                n_tc = np.count_nonzero(~arr_tc.mask[year-2001, :, :])
            else:
                mean_tc = np.nan
                n_tc = 0

            # add all these to a list
            out = [obj_id, year, region, ign2, mean_depth, n_depth, mean_tc, n_tc]

            # loop throug other images and add to list
            for img_path in img_paths:
                img1 = rasterio.open(img_path)
                if not img1.nodata:
                    arr1, trans = rasterio.mask.mask(img1, [geom], crop=True, filled=False)
                    if False in arr_tc.mask:
                        mean = np.nanmean(arr1, axis=(1, 2))
                        n = np.count_nonzero(~arr1.mask[0, :, :])
                    else:
                        mean = [np.nan] * arr1.shape[0]
                        n = 0
                else:
                    arr1, trans = rasterio.mask.mask(img1, [geom], crop=True)
                    arr1 = np.where(arr1 == img1.nodata, np.nan, arr1)
                    mean = np.nanmean(arr1, axis=(1, 2))
                    n = np.count_nonzero(~np.isnan(arr1[0, :, :]))
                out.extend(mean)
                out.append(n)

            # write data into csv
            out = tuple(out)
            line = '%s,%i,%s,%s,%f,%i,%f,%i,%f,%f,%f,%f,%f,%f,%i,%f,%i,%f,%i,%f,%f,%f,%f,%i\n' %out
            f.write(line)
# close file
f = None
