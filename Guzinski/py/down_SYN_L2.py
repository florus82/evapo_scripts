import geopandas as gpd
import pandas as pd
import shapely
import os
import openeo
import matplotlib.pyplot as plt
import xarray as xr
import time
connection = openeo.connect("openeo.dataspace.copernicus.eu").authenticate_oidc()

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

    for date in dates[:]:  # Iterates over a copy of the list to allow removal at the end without messing up the order

            storPath = f"/data/Aldhani/eoagritwin/et/Auxiliary/S3_SYN_L2_NDVI/raw/Germany_{date[0]}.nc"
            print(storPath)
            
            if os.path.exists(storPath):
                t = time.localtime()
                ti = time.strftime("%H:%M:%S", t)
                print(f"already exists - next one at {ti}")

            else:
                try:
                    cube = connection.load_collection(
                    "SENTINEL3_SYN_L2_SYN",
                    spatial_extent = germany,
                    temporal_extent = date,
                    bands=["Syn_Oa08_reflectance","Syn_Oa17_reflectance", "CLOUD_flags"]
                    )
                    
                    # clouds = cube.band("CLOUD_flags")
                    # clear_mask = (clouds.bitwise_and(1) == 0)
                    # masked = cube.mask(clear_mask)
    
                    # # calculate ndvi according to https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel-3/ndvi/
                    # ndvi = (masked.band("Syn_Oa17_reflectance") - masked.band("Syn_Oa08_reflectance")) / (masked.band("Syn_Oa17_reflectance") + masked.band("Syn_Oa08_reflectance"))
                    # ndvi = ndvi.save_result(format="NETCDF")
                    cube.download(storPath)
                    dates.remove(date)

                except Exception as e:
                    print(e)
                    t = time.localtime()
                    ti = time.strftime("%H:%M:%S", t)
                    print(f"thrown at {ti}")
                    continue