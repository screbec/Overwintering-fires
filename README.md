# Overwintering fires algorithm & analyses
Python and R code for the identification and analysis of overwintering fires

This code uses ground truth and and satellite data to identify fires that have overwintered and subsequently analyse spatial and temporal drivers of overwintering. The results of this work are curently under review (Scholten et al., 2021). The datasets required for the code are listed in the data requirements. They are all freely available and linked in the data availability statement of Scholten et al. (2021). 
  
The scripts include the preprocessing steps (P) and the analyses (A) as shown in the workflow below. The source code for all display items in the article (Figures 1-4 and Table 1) is also included. All preprocessing was done in python, whereas analyses and plots were coded in R. For easier readibility functions used in each script are pasted directly in the corresponding script. Lastly, please note that this code was not written by a professional software developer, so it may not be written in the most beautiful and effective way possible. If you have any comments or suggestions regarding the code, please share them with us. Also feel free to contact us if you have any questions about the code, data or the analysis in general.

### Workflow: 
(A marks the actual analysis script, P marks preprocessing steps that generate the data needed for A)
#### 1. Algorithm for identifying overwintering fires:
- P Compute first snowfree day from modscag data (pre_algorithm/modscag_ldos.py) 
- P Find nearest antecedent fire scar for each ignition point and calculate distance (pre_algorithm/nndist.py)
- P Spatio-temporal intersection of ignition points with lightning strikes (pre_algorithm/lightning.py)
- P Spatial intersection and distance calculation with roads within 5 km buffer (pre_algorithm/roads.py)
- A Identify and apply the thresholds for the identification of overwintering fires (algorithm.R)
#### 2. Analysis of temporal drivers
- P generate 1975-2018 burned area time series from fire perimeters (pre_drivers/ba_long.py)
- P extract regional MJJAS temperature indices (mean, hot days, 90th percentile) (pre_drivers/temp_extract.py)
- P extract burn depth (pre_drivers/spatial_extract.py)
- A temporal drivers analysis (including Fig. 2, ED Fig. 6, ED Fig.7) (temporal_drivers.R)
#### 3. Analysis of spatial drivers
- P extract burn depth & environmental variables for fire scars (pe_drivers/spatial_extract.py)
- A spatial drivers analysis (incl. Fig. 3, Table 1) (spatial_drivers.R)
#### 4. Maps
- code for Figure 1 (mapping/figure1.R)
- code for Figure 4 (mapping/figure4.R)

### Data requirements:
- ignitions (shapefile), i.e. from AKFED v3
- fire perimeters (shapefile), i.e. from Alaska Wildfire Maps and Canada National Fire Databases
- MODSCAG snow fraction data (tiffs)
- infrastructure data & 5 km buffer (shapefiles)
- lightning strikes (shapefiles, i.e. from Alaska lightning detection network)
- temperature and fire weather data, i.e. from NARR
- burn depth, i.e. from AKFED
- environmental variables, i.e. tree species, tree cover, elevation/slope, soil carbon, etc.
- shapefiles of Alaska, Northwest Territories and interior AK/NWT

### Software requirements:
- Python code tested with Anaconda Python 3.7.8, required packages:
  - gdal 2.3.3
  - fiona 1.8.4
  - rtree 0.9.4
  - numpy 1.19.1
  - pandas 1.1.0
  - shapely 1.6.4
- R code tested with R 4.0.2, required packages:
  - tidyverse 1.3.0
  - raster 3.0-12
  - sf 0.8-1
  - egg 0.4.5 (for labelled multi-facet figures)
  - ggspatial 1.1.4 (for scale bar)
  - viridis 0.5.1 (for colour palette of Figure 4)
