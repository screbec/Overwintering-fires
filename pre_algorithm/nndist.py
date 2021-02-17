# -*- coding: utf-8 -*-
"""
Finds the nearest neighboring fire perimeter of the antecedent fire seasone 
for each ignition point using an RTree index
and calculates the nearest distance between both polygons
returns results as csv file

Requirements:
- shapefile of fire perimeters
- shapefile of ignition points

Output: csv file listing the closest fire perimeter of teh antecedent year
for each ignition point, as well as the distance between both

@author: RCScholten
"""

perimeter_file = "perimeter.shp" # this is the filename and path of the fire perimeters
ip_file = "ignitions.shp" # filename and path to the ignitions shapefile
outpath = 'nndist.csv' # this is the output file which is later read in 04_algorithm.R for the algorithm

# import modules
import pandas as pd
import fiona
from shapely.geometry import shape
import rtree

#%% functions
def readFilter(file, attribute, value, IDattribute):
    '''reads in a shapefile using fiona, filters it by the attribute and value 
    given and returns a list of ids and geometries'''
    with fiona.open(file) as src:
        filtered = filter(lambda f: f['properties'][attribute]==value, src)
        out_geoms = [(poly['properties'][IDattribute], shape(poly['geometry'])) for poly in filtered]
    return(out_geoms)

def BuildTree(geoms, FID=False):
    '''Builds Rtree from a shapely multipolygon shape
    if FID is true, the FID is used as object identifier
    else the idattribute in geoms is used'''
    idx = rtree.index.Index()
    for ind, poly_tuple in enumerate(geoms):
        oid, poly = poly_tuple
        if FID:
            idx.insert(ind, poly.bounds, ind)
        else:
            idx.insert(ind, poly.bounds, oid)
    return idx

def nndistpoly(old_geoms, new_geoms, spatial_index, year):
    ''' 
    Find Polygons/MultiPolygons that are the three nearest neighbors of a target 
    polygon using a RTree spatial index (based on the bounding boxes) and 
    calculate the distance between the nearest points of the nearest polygons 
    '''
    # create a list for the output
    result = list()
    # create dictionary from old_geoms to retrieve geometry by Fireid
    for idx, poly_tuple in enumerate(new_geoms):      
        srcid, poly = poly_tuple
        # return three nearest objects by fireid
        nn = list(spatial_index.nearest(poly.bounds, 3, objects=True))
        fids, bbox = zip(*[(item.object, item.bbox) for item in nn])
        dist = list()
        # calculate nearest points and distance between nearest points for the polygons
        for fid in fids:
            tgtid, multipoly = old_geoms[fid]
            dist.append((multipoly.distance(poly), tgtid))
            nndist = min(dist, key = lambda t: t[0])
        # return the source fireid, fireid of nn and nearest distance of nn
        result.append([srcid, year, nndist[1], nndist[0]])   
    return(result)

#%% nearest distance to fire of previous year
if __name__ == '__main__':
    # create output list
    res = list()
    resIP = list()
    # loop through the years
    for y in range(2002,2019):
        newgeoms = readFilter(ip_file, 'Year', y, 'ID')
        oldgeoms = readFilter(perimeter_file, 'Year', y - 1, 'FireID')
        # Create spatial index
        spatidx = BuildTree(oldgeoms, FID = True)
        # find nearest neighbour and calculate distance
        res.append(nndistpoly(oldgeoms, newgeoms, spatidx, y))
    # remove 'None' entries
    res = [x for x in res if x is not None]
    res = [item for sublist in res for item in sublist]
    # convert list to pd.dataframe and join with shapefile
    df = pd.DataFrame(res, columns=['ID','Year','FireID_tgt','nndist'])
    # save to csv
    df.to_csv(outpath)





