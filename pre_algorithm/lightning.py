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

# import required modules
import pandas as pd
from osgeo import ogr
import rtree

# inputs
outfile = 'lightning.csv' # this is the output file which is later read in 04_algorithm.R for the algorithm
file_ip = 'ignitions.shp'  # ignitions shapefile
file_light = 'lightning.shp' # lightning shapefile

#%% functions
def build_tree_gdal(lyr, id_fld='FID', *args):
    '''builds Rtree from an ogr layer (can be any geometry type)
    buffer distance and buffer field can be specified as args[0,1]'''
    idx = rtree.index.Index()
    for fid1, feat in enumerate(lyr):
        if args: # buffer specifications are given in args
            if len(args) > 1:
                # extract buffer distance target feature
                buffer_dist = feat.GetFieldAsDouble(args[1])
            else:
                buffer_dist = 0
            # apply buffer
            geom = feat.GetGeometryRef().Buffer(args[0] + buffer_dist)
        else:
            geom = feat.GetGeometryRef()
        xmin, xmax, ymin, ymax = geom.GetEnvelope()
        if id_fld == 'FID':
            fldval = feat.GetFID()
        else:
            fldval = feat.GetField(id_fld)
        idx.insert(fid1, (xmin, ymin, xmax, ymax), fldval)
    # reset reading of layer in case it is used at a later stage
    lyr.ResetReading()
    return idx

def filter_shp(value, filter_field, in_shapefile):
    '''reads and filters a shapefile
    value and filter_field can be lists if multiple filters should be applied'''
    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.Open(in_shapefile, 0)
    layer = data_source.GetLayer()
    # Filter by our query
    if isinstance(filter_field, list):
        i = 0
        while i < len(filter_field):
            if isinstance(value[i], str):
                query_temp = "{} = '{}'".format(filter_field[i], value[i])
            else:
                query_temp = "{} = {}".format(filter_field[i], value[i])
            try:
                query_str = query_str + ' and ' + query_temp
            except:
                query_str = query_temp
            i += 1
    else:
        if isinstance(value, str):
            query_temp = "{} = '{}'".format(filter_field, value)
        else:
            query_str = '{} = {}'.format(filter_field, value)
    layer.SetAttributeFilter(query_str)
    return data_source, layer

def intersect_tree(ip_args, light_args):
    '''Takes a target layer and a rtree index and its corresponding layer,
    a field in the feature that is used for calculating the buffer distance (dE), and
    the attribute names of the target attributes of both layers.
    Returns 1 if there is an intersection for the date and 0 if there isn't
    Caution: does not check for real intersection, just intersection of Envelopes!'''
    # read all needed parameters into variables
    feat, doy_att = ip_args # ipArgs contains the feature and the name of the day of burning attribute field
    idx_light, light_lyr, light_att = light_args # LightArgs contains the rtree index, the lightning layer and the name of the day of lightning strike attribute field

    # extract doy and geometry envelope from ignition points
    ip_doy = feat.GetFieldAsInteger(doy_att)
    geom = feat.GetGeometryRef()
    xmin, xmax, ymin, ymax = geom.GetEnvelope()

    # check for spatiotemporal overlap with lightning data
    items = idx_light.intersection((xmin, ymin, xmax, ymax), objects=True)
    fids_filter = [n.object for n in items]

    # for each lightning fid check if doys overlap
    result = None
    light_doys = [light_lyr.GetFeature(fid).GetFieldAsInteger(light_att) for fid in fids_filter]
    if light_doys:
        result = [ip_doy - light_doy for light_doy in light_doys
                  if (light_doy < (ip_doy + 2) and light_doy > (ip_doy - 16))]
        if result:
            if -1 in result and len(result) > 1:
                result.remove(-1)
            result = min(result)

    return result

#%% main
if __name__ == '__main__':

    # create output list
    res1 = []
    # loop through the years and assign lightning ignitions
    for y in range(2001, 2019):

        # load filtered datasets
        ip_ds, ip_lyr = filter_shp(y, 'Year', file_ip)
        l_ds, l_lyr = filter_shp(y, 'year', file_light)

        # build rtrees, for lightning points apply buffer
        ltree = build_tree_gdal(l_lyr, 'FID', 2000)

        # loop though fireID layer, do intersections
        for feat1 in ip_lyr:
            res_light = intersect_tree([feat1, 'ID', 'doy'], [ltree, l_lyr, 'doy'])
            res1.append([feat1.GetFieldAsInteger('ID'), res_light])

    # convert to pandas dataframa and write to csv
    df = pd.DataFrame(res1, columns=['ID', 'lightning'])
    df.to_csv(outfile)
