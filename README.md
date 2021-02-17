# Overwintering fires algorithm & analyses (Scholten et al. 2021)
Python and R code for the identification of overwintering fires from AKFED burned area and emissions

Data requirements:
- ignitions (shapefile), i.e. from AKFED
- fire perimeters (shapefile), i.e. from AWFM or CNFD (Alaska and Canada Fire Databases)
- MODSCAG snow fraction data (tiffs)
- infrastructure data & 5 km buffer (shapefiles)
- lightning strikes (shapefiles)
- temperature and fire weather data, i.e. from NARR
- burn depth, i.e. from AKFED
- environmental variables, i.e. tree species, tree cover, elevation/slope, soil carbon, etc.
- shapefiles of Alaska, Northwest Territories and interior AK/NWT

The scripts include the preprocessing steps (P) and the analysis (algorithm/statistics/figures) (A) as shown in the workflow below.
For easier readibility functions used in each script are directly pasted in the script.

Workflow: 
(A marks the actual analysis script, P marks preprocessing steps that generate the data needed for A)
1. Algorithm for identifying overwintering fires:
P Compute first snowfree day from modscag data (pre_algorithm/modscag_ldos.py) 
P Find nearest antecedent fire scar for each ignition point and calculate distance (pre_algorithm/nndist.py)
P Spatio-temporal intersection of ignition points with lightning strikes (pre_algorithm/lightning.py)
P Spatial intersection and distance calculation with roads within 5 km buffer (pre_algorithm/roads.py)
A Identify and apply the thresholds for the identification of overwintering fires (algorithm.R)
2. Temporal drivers
P generate 1975-2018 burned area time series from fire perimeters (pre_drivers/ba_long.py)
P extract regional MJJAS temperature indices (mean, hot days, 90th percentile) (pre_drivers/temp_extract.py)
P extract burn depth (pre_drivers/spatial_extract.py)
A temporal drivers analysis (including Fig. 2, ED Fig. 6, ED Fig.7) (temporal_drivers.R)
3. Spatial drivers
P extract burn depth & environmental variables for fire scars (pe_drivers/spatial_extract.py)
A spatial drivers analysis (incl. Fig. 3, Table 1) (spatial_drivers.R)
4. Maps
A code for Figure 1 (mapping/figure1.R)
A code for Figure 4 (mapping/figure4.R)

Software requirements:
- python code tested with Python 3.7.8, required packages:
    gdal 2.3.3
    fiona 1.8.4
    rtree 0.9.4
    numpy 1.19.1
    pandas 1.1.0
    shapely 1.6.4
- R code tested with R 4.0.2, required packages:
    tidyverse 1.3.0
    raster 3.0-12
    sf 0.8-1
    egg 0.4.5 (for labelled multi-facet figures)
    ggspatial 1.1.4 (for scale bar)
    viridis 0.5.1 (for colour palette of Figure 4)
