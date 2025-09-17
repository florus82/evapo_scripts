import cdsapi
import os 
import time 

base_path = '/data/Aldhani/eoagritwin/et/Auxiliary/AOD_CAMS/'
client = cdsapi.Client()

dataset = "cams-global-reanalysis-eac4"
# Define the dataset and request parameters

years = [str(y) for y in range(2017,2025,1)]
months = [m for m in range(1,13,1)]
for year in years:
    for month in months:
        try:
            request = {
                "variable": ["total_aerosol_optical_depth_550nm"],
                "date": [f'{year}-{month:02d}-01/{year}-{month:02d}-29'],

                "time": [
                    "00:00", "03:00", "06:00", "09:00",
                    "12:00", "15:00", "18:00", "21:00", 
                ],
                'data_format': 'grib',
                "area": [56, 5, 47, 16]
            }


            client.retrieve(dataset, request, f'{base_path}AOD550_{year}_{month}.grib')
        
        except Exception as e:
            print(e)
            t = time.localtime()
            ti = time.strftime("%H:%M:%S", t)
            print(f"thrown at {ti}")
            continue


