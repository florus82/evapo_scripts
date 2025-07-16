import geopandas as gpd
import pandas as pd
import shapely
import os
import openeo
import matplotlib.pyplot as plt
import xarray as xr
import time



connection = openeo.connect("openeo.dataspace.copernicus.eu").authenticate_oidc()

connection.describe_collection("SENTINEL3_SLSTR_L2_LST")

years = [year for year in range(2017, 2025, 1)]
months = [f'{month:02d}' for month in range(1, 13, 1)]
dates = []
for year in years:
    for month in months:
        next_month = int(month) + 1
        if next_month > 12:
            next_month = 1
            next_year = year + 1
        else:
             next_year = year
        dates.append([f'{year}-{month}-01', f'{year}-{month}-15'])
        dates.append([f'{year}-{month}-15', f'{next_year}-{next_month:02d}-01'])

germany = {"west": 5.592041, "south": 47.129951, "east": 15.26001, "north": 55.09723}

while dates:

    for date in dates:  # Iterates over a copy of the list to allow removal at the end without messing up the order

            storPath = f"/data/Aldhani/eoagritwin/et/Auxiliary/VZA/raw/Germany_{date[0]}.nc"
            print(storPath)
            
            if os.path.exists(storPath):
                t = time.localtime()
                ti = time.strftime("%H:%M:%S", t)
                print(f"already exists - next one at {ti}")

            else:
                try:
                    sentinel3_cube = connection.load_collection(
                    "SENTINEL3_SLSTR_L2_LST",
                    spatial_extent = germany,
                    temporal_extent = date,
                    bands=["viewZenithAngles"]#, "exception", "confidence_in"]
                    )
                    
                    # band1 = sentinel3_cube.band('confidence_in')
                    # band2 = sentinel3_cube.band('viewZenithAngles')
                    # band3 = sentinel3_cube.band('exception')

                    # cloud_mask = (band1 >= 16384) | (band2 <= -32768) | (band3 != 0)
                    # sentinel3_cube_masked = sentinel3_cube.mask(cloud_mask)
                    # sentinel3_cube_masked = sentinel3_cube_masked.filter_bands(['viewZenithAngles'])
                    # sentinel3_cube_masked.download(storPath)
                    sentinel3_cube.download(storPath)
                    dates.remove(date)

                except Exception as e:
                    print(e)
                    t = time.localtime()
                    ti = time.strftime("%H:%M:%S", t)
                    print(f"thrown at {ti}")
                    continue