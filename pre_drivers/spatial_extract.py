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
from rasterio.mask import mask
import fiona
import numpy as np
import pandas as pd

#%% main
# load list of overwintering ids (output from algorithm)
with open('overwintering_ids.txt', 'r') as f:
   hoids = f.read().splitlines()
hoids = [int(hoid) for hoid in hoids[1:]]

# extract list of their original fire perimeter ID
id2fireid = pd.read_csv('nndist.csv')
fireids  = list(set(id2fireid[id2fireid['ID'].isin(hoids)].FireID_tgt))

# raster files for clipping (these are the same for AK and NWT)
depthpath = 'depth_20012018.tif'
tcpath = 'tc_500_20002017.tif' # take the tree cover of the year before the fire

# start writing csv file
csv = 'spatial.csv'
f = open(csv, "w")
f.write('FireID,Year,ShortNm,ign2,depth_mean,n_depth,tc,n_tc,sp_bs,sp_ws,sp_dc,sp_gs,sp_nv,sp_p,n_ts,Elev,n,slope,n,soil30,soil100,soil200,soil300,n\n')

# loop through regions   
for region in ['NT', 'AK']:

    # inputs and outputs according to region (these are specific to each region)
    if region == 'NT':
        shapepath = 'perimeters_NT.shp'
        imgpaths = ['treespecies_nwt.tif', 
                    'arcticdem_NWT.tif',
                    'arcticdem_slope_NWT.tif',
                    'ncscd_ytnt.tif']

    elif region == 'AK':
        shapepath = 'perimeters_AK.shp'
        imgpaths = ['treespecies_ak.tif',
                    'arcticdem_AK.tif',
                    'arcticdem_slope_AK.tif',
                    'ncscd_ak.tif']
           
    # start clipping
    # open images with rasterio
    depthimg = rasterio.open(depthpath)
    tcimg = rasterio.open(tcpath)

    # iteration through features
    with fiona.open(shapepath, "r") as shapefile:
        for feat in shapefile:
            # extract relevant information from feature
            geom = feat["geometry"]
            objID = feat["properties"]["FireID"]
            year = int(feat["properties"]["Year"])
            # determine whether the feature produced an overwintering fire or not
            ign2 = 'overwinter' if objID in fireids else 'other'
            
            # year 2018 has to be skipped since we don't have burned area for 2019
            # and consequently don't know if fire scars in 2018 caused holdovers in 2019
            if year > 2017:
                continue                

            # extract tiff values masking with feature
            arrdepth, trans1 = rasterio.mask.mask(depthimg, [geom], crop=True) #, all_touched=True)
            arrtc, trans3 = rasterio.mask.mask(tcimg, [geom], crop=True, filled=False) #, all_touched=True)

            # convert NoData to NAN
            arrdepth = np.where(arrdepth==0.0, np.nan, arrdepth)

            # calculate stats on target year
            mean_depth = np.nanmean(arrdepth[year-2001,:,:])     # 2001 = year0
            n_depth = np.count_nonzero(~np.isnan(arrdepth[year-2001,:,:]))
            # rasterio assigns 0 as nodataval if the nodatavalue is missing in a raster
            # since 0 is a valid value for our rasters we need a workaround
            # we return the masked array (l 108, filled = False) 
            # and count the number of 'False' = not masked values in the mask for n
            # if the whole array is masked (geom does not hit the center point of a raster pixel)
            # it cannot calculate the nanmean, therefore the if-clause
            if False in arrtc.mask:
                mean_tc = np.nanmean(arrtc[year-2001,:,:])           # 2000 = year0
                n_tc = np.count_nonzero(~arrtc.mask[year-2001,:,:])
            else:
                mean_tc = np.nan
                n_tc = 0
            
            # add all these to a list
            out = [objID, year, region, ign2, mean_depth, n_depth, mean_tc, n_tc]
            
            # loop throug other images and add to list
            for imgpath in imgpaths:
                img = rasterio.open(imgpath)
                if not img.nodata:
                    arr, trans = rasterio.mask.mask(img, [geom], crop=True, filled=False)
                    if False in arrtc.mask:
                        mean = np.nanmean(arr, axis = (1,2))
                        n = np.count_nonzero(~arr.mask[0,:,:])
                    else:
                        mean = [np.nan] * arr.shape[0]
                        n = 0
                else:
                    arr, trans = rasterio.mask.mask(img, [geom], crop=True)
                    arr = np.where(arr==img.nodata, np.nan, arr)
                    mean = np.nanmean(arr, axis = (1,2))
                    n = np.count_nonzero(~np.isnan(arr[0,:,:]))
                out.extend(mean)
                out.append(n)

            # write data into csv
            out = tuple(out)
            line = '%s,%i,%s,%s,%f,%i,%f,%i,%f,%f,%f,%f,%f,%f,%i,%f,%i,%f,%i,%f,%f,%f,%f,%i\n' %out
            f.write(line)
# close file
f = None











