import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.guzinski import *
from helperToolz.dicts_and_lists import INT_TO_MONTH
import cfgrib
from joblib import Parallel, delayed
import time

os.environ["PROJ_LIB"] = "/data/Aldhani/users/potzschf/conda/envs/cds_era5/share/proj"

# get variables downloaded
storPath = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/'
base_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/grib/'
dummy_path = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/daily_observations_all/2019/Daily_LST_max_2019_August.tif'
dem_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/DEM_GER_LST_WARP.tif'
geopot_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/low_res/geopotential/geopotential_low_res.tif'
variables = [name for name in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, name))]
# tiles = getFilelist(dem_path, '.tif')

#for variable in variables:

for variable in ['2m_temperature']:
    if variable != 'geopotential': # geopotential is constant over time --> only a single layer needed
        for year in [2019]:#range(2017, 2025, 1):
            files = [file for file in getFilelist(os.path.join(base_path, variable), '.grib') if str(year) in file]
            for file in files:
                # make folder to store data and subset files only to growing season
                outDir = f'{storPath}{'low_res'}/{variable}/{year}/'
                os.makedirs(outDir, exist_ok=True)
                
                m = int(file.split('_')[-1].split('.')[0])
                if growingSeasonChecker(m):
                    month = INT_TO_MONTH[f'{m:02d}']
                    outPath = f'{outDir}{variable}_{'low_res'}_{year}_{month}.tif'
                    if os.path.exists(outPath):
                        print(f'{variable} already processes for {month}/{year}')
                        continue
                    else:
                        warp_ERA5_to_reference(file, dummy_path, outPath, sharp_DEM=dem_path,
                                               sharp_geopot=geopot_path,
                                               sharp_blendheight=100, sharp_rate=STANDARD_ADIABAT)
    
    else:
        files = [file for file in getFilelist(os.path.join(base_path, variable), '.grib')]
        outDir = f'{storPath}{'low_res'}/{variable}/'
        os.makedirs(outDir, exist_ok=True)
        outPath = f'{outDir}{variable}_{'low_res'}.tif'
        if os.path.exists(outPath):
            print(f'{variable} already processes')
            continue
        else:
            warp_ERA5_to_reference(files[0], dummy_path, outPath, bandL=[1])
