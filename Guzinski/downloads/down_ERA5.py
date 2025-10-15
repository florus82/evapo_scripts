import cdsapi
import os 
import time 

base_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/'
client = cdsapi.Client()

dataset = "reanalysis-era5-single-levels"

variables = [             
        "2m_dewpoint_temperature", # vapor pressure
        "2m_temperature",
        "surface_pressure",
        "100m_u_component_of_wind",
        "100m_v_component_of_wind",
        "total_column_water_vapour",
        "geopotential",
        "surface_solar_radiation_downward_clear_sky"]

years = [str(y) for y in range(2017,2025,1)]

months = [str(y) for y in range(1,13,1)]


for variable in variables:
    for year in years:
        for month in months:
            varPath = f'{base_path}{variable}'
            if not os.path.exists(varPath):
                os.makedirs(varPath)
            storPath = f'{varPath}/{variable}_{year}_{month}.grib'

            if os.path.exists(storPath):
                t = time.localtime()
                ti = time.strftime("%H:%M:%S", t)
                print(f"already exists - next one at {ti}")

            else:
                print()
                try:
                    request = {
                        "product_type": ["reanalysis"],
                        "variable": [
                            variable
                        ],
                        "year": [
                            year    
                        ],
                        "month": [
                            month
                        ],
                        "day": [
                            "01", "02", "03",
                            "04", "05", "06",
                            "07", "08", "09",
                            "10", "11", "12",
                            "13", "14", "15",
                            "16", "17", "18",
                            "19", "20", "21",
                            "22", "23", "24",
                            "25", "26", "27",
                            "28", "29", "30",
                            "31"
                        ],
                        "time": [
                            "00:00", "01:00", "02:00",
                            "03:00", "04:00", "05:00",
                            "06:00", "07:00", "08:00",
                            "09:00", "10:00", "11:00",
                            "12:00", "13:00", "14:00",
                            "15:00", "16:00", "17:00",
                            "18:00", "19:00", "20:00",
                            "21:00", "22:00", "23:00"
                        ],
                        "data_format": "grib",
                        "download_format": "unarchived",
                        "area": [56, 5, 47, 16]
                    }


                    target = storPath

                    client.retrieve(dataset, request, target)
                
                except Exception as e:
                    print(e)
                    t = time.localtime()
                    ti = time.strftime("%H:%M:%S", t)
                    print(f"thrown at {ti}")
                    continue