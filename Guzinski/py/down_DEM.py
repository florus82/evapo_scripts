import geopandas as gpd
import pandas as pd
import shapely
import os
import openeo
import matplotlib.pyplot as plt
import xarray as xr
import time
import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *

connection = openeo.connect("openeo.dataspace.copernicus.eu").authenticate_oidc()

while len(getFilelist('/data/Aldhani/eoagritwin/et/Auxiliary/DEM/raw_tiles/', '.tif')) < 99:
    for i in range(0,11,1):
        for j in [0,1]:# in range(0,10,1):
            aoi = {'west': 5+i, 'south': 46+j, 'east': 6+i, 'north': 47+j}
            storPath = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/raw_tiles/DEM_GER_Last_{j}_{i}.tif'
            print(storPath)
            
            if os.path.exists(storPath):
                    t = time.localtime()
                    ti = time.strftime("%H:%M:%S", t)
                    print(f"already exists - next one at {ti}")

            else:
                try:
                    dem_cube = connection.load_collection(
                                    "COPERNICUS_30",
                                    spatial_extent = aoi,
                                    bands=["DEM"]
                                    )                               
                    dem_cube.download(storPath)

                except Exception as e:
                    print(e)
                    t = time.localtime()
                    ti = time.strftime("%H:%M:%S", t)
                    print(f"thrown at {ti}")
                    continue

