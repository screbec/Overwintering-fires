# -*- coding: utf-8 -*-
"""
Takes a target shapefile, a line shapefile and its buffered version,
intersects the buffer with the polygons using an rtree index
if there is an intersection: calculates minimal distance and returns it
together with an ID attribute value (Infrastruture type)
Requirements:
    - ignition points (i.e. from AKFED), including positional uncertainty
    - infrastructure shapefile (unbuffered and buffered to i.e. 5km)
Output: csv of all ignitions that intersect with an infrastructure element,
the type of infrastructure it intersects with, and the distance between both

@author: RCScholten
"""


outname = 'distroads.csv' # this is the output file which is later read in 04_algorithm.R for the algorithm
roads = "infrastructure.shp" # this is the infrastructure shapefile
roadsbuf = 'infra_buffer5km.shp' # this is a polygon shapefile of the infrastructures with a 5km buffer
roadvar = 'INFRA_TYPE' # the script also outputs the type of infrastructure of every intersection. This is the attribute of the shapefile pointing to the infrastructure type
ignpoints = "ignitions.shp" # this is the ignition point shapefile


# import modules
import fiona
from shapely.geometry import shape
import pandas as pd
import rtree

#%% functions
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

def nndistline(line_geoms, poly_geoms, buffer_index):
    ''' 
    Find Polys that intersect the line buffer,
    return the minimal distance between line and poly
    and the infratype (ID variable of line_geoms)
    '''
    # create a list for the output
    result = list()
    for poly_tuple in poly_geoms: 
        srcid, poly = poly_tuple
        # return fids of buffer intersecting with polygon
        nn = list(buffer_index.intersection(poly.bounds, objects=True))
        if len(nn) > 0:
            fids, bbox = zip(*[(item.object, item.bbox) for item in nn])
            dist = list()
            for fid in fids:
                infratype, line = line_geoms[fid]
                dist.append((line.distance(poly), infratype))
            nndist = min(dist, key = lambda t: t[0])
        else:
            nndist = (None, None)
        # return the source id, nearest distance, and infratype
        result.append([srcid, nndist[0], nndist[1]]) 
    return(result)
    
    
#%% main
if __name__ == '__main__':
    
    # read inputs
    road_geoms = [(poly['properties'][roadvar], shape(poly['geometry'])) for poly in fiona.open(roads)]
    buffer_geoms = [(poly['properties'][roadvar], shape(poly['geometry'])) for poly in fiona.open(roadsbuf)]
    ip_geoms = [(poly['properties']['ID'], shape(poly['geometry'])) for poly in fiona.open(ignpoints)]
    
    # find intersections and calculate nndist
    spatidx = BuildTree(buffer_geoms, FID=True)
    distIP = nndistline(road_geoms, ip_geoms, spatidx)
        
    # convert to pandas dataframa and write to csv
    dfIP = pd.DataFrame(distIP, columns=['ID','distRoad','Infratype'])
    dfIP.to_csv(outname)
    