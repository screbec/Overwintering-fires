# -*- coding: utf-8 -*-
"""
Intersection of ignition points with lightning strikes
All lightning strikes with a specified buffer distance (i.e. 2500 m) 
which occur within one day after and 16 days before the ignition
are recorded including their temporal lag time (lightning strike - ignition)
This serves for setting a temporal threshold for lightning ignitions.
Requirements:
    - shapefile with lightning strikes
    - shapefile with ignition points (i.e. from AKFED), including positional uncertainty
    
Output: csv file of all ignition points that intersect with a lightning strike
and the temporal lag between lightning strike and ignition

@author: RCScholten
"""

outfile = 'lightning.csv' # this is the output file which is later read in 04_algorithm.R for the algorithm
fileIP = 'ignitions.shp'  # ignitions shapefile
fileLight = 'lightning.shp' # lightning shapefile

# import required modules
import pandas as pd
from osgeo import ogr
import rtree

#%% functions
def BuildTreeGdal(lyr, IDfld = 'FID', *args):
    '''builds Rtree from an ogr layer (can be any geometry type)
    buffer distance and buffer field can be specified as args[0,1]'''
    idx = rtree.index.Index()
    for fid1, feat in enumerate(lyr):
        if args: # buffer specifications are given in args
            if len(args) > 1:
                # extract buffer distance target feature
                bufferDist = feat.GetFieldAsDouble(args[1])
            else:
                bufferDist = 0
            # apply buffer
            geom = feat.GetGeometryRef().Buffer(args[0] + bufferDist) 
        else:
            geom = feat.GetGeometryRef()
        xmin, xmax, ymin, ymax = geom.GetEnvelope()
        if IDfld =='FID':
            fldval = feat.GetFID()
        else: 
            fldval = feat.GetField(IDfld)
        idx.insert(fid1, (xmin, ymin, xmax, ymax), fldval)
    # reset reading of layer in case it is used at a later stage
    lyr.ResetReading()
    return idx

def filterShp(value, filter_field, in_shapefile):
    '''reads and filters a shapefile
    value and filter_field can be lists if multiple filters should be applied'''
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(in_shapefile, 0)
    layer = dataSource.GetLayer()
    # Filter by our query
    if type(filter_field) == list:
        i = 0
        while(i < len(filter_field)):
            if type(value[i]) == str:
                query_temp = "{} = '{}'".format(filter_field[i], value[i])
            else:
                query_temp = "{} = {}".format(filter_field[i], value[i])
            try:
                query_str = query_str + ' and ' + query_temp
            except:
                query_str = query_temp
            i += 1      
    else:
        if type(value) == str:
            query_temp = "{} = '{}'".format(filter_field, value)
        else:
            query_str = '{} = {}'.format(filter_field, value)
    layer.SetAttributeFilter(query_str)
    return dataSource, layer


def intersectTree(ipArgs, LightArgs):
    '''Takes a target layer and a rtree index and its corresponding layer, 
    a field in the feature that is used for calculating the buffer distance (dE), and
    the attribute names of the target attributes of both layers.
    Returns 1 if there is an intersection for the date and 0 if there isn't
    Caution: does not check for real intersection, just intersection of Envelopes!'''
    # read all needed parameters into variables
    feat, doyAtt = ipArgs # ipArgs contains the feature and the name of the day of burning attribute field
    idxLight, Lightlyr, LightAtt = LightArgs # LightArgs contains the rtree index, the lightning layer and the name of the day of lightning strike attribute field
    
    # extract doy and geometry envelope from ignition points
    ip_doy = feat.GetFieldAsInteger(doyAtt)
    geom = feat.GetGeometryRef()
    xmin, xmax, ymin, ymax = geom.GetEnvelope()

    # check for spatiotemporal overlap with lightning data
    items = idxLight.intersection((xmin, ymin, xmax, ymax), objects=True)
    fidsFilter = [n.object for n in items] 
    
    # for each lightning fid check if doys overlap
    result = None
    lightdoys = [Lightlyr.GetFeature(fid).GetFieldAsInteger(LightAtt) for fid in fidsFilter]
    if len(lightdoys) > 0:
        result = [ip_doy - lightdoy for lightdoy in lightdoys 
                 if (lightdoy < (ip_doy + 2) and lightdoy > (ip_doy - 16))]
        if len(result) > 0:
            if (-1) in result and len(result) > 1:
                result.remove(-1)
            result = min(result)

    return result

#%% main
if __name__ == '__main__':
    
    # create output list
    res = []
    # loop through the years and assign lightning ignitions
    for y in range(2001,2019):
        
        # load filtered datasets
        ip_ds, ip_lyr = filterShp(y, 'Year', fileIP)
        l_ds, l_lyr = filterShp(y, 'year', fileLight)
        
        # build rtrees, for lightning points apply buffer 
        ltree = BuildTreeGdal(l_lyr, 'FID', 2000) 
        
        # loop though fireID layer, do intersections
        for feat in ip_lyr:      
            resLight = intersectTree([feat, 'ID', 'doy'], [ltree, l_lyr, 'doy'])
            res.append([feat.GetFieldAsInteger('ID'), resLight])

    # convert to pandas dataframa and write to csv
    df = pd.DataFrame(res, columns=['ID','lightning'])
    df.to_csv(outfile)





