import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.evapo import *
from helperToolz.dicts_and_lists import INT_TO_MONTH
import cfgrib
os.environ["PROJ_LIB"] = "/data/Aldhani/users/potzschf/conda/envs/cds_era5/share/proj"

# get variables downloaded
storPath = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/'
base_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/grib/'
dummyTIFF = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/daily_observations_all/2020/Daily_LST_medians_2020_September.tif'
variables = [name for name in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, name))]


# check out 2m_temperature
#for variable in variables:
variable = '2m_temperature'
for year in range(2017, 2025, 1):
    print(f'processing year {year}')
    files = [file for file in getFilelist(os.path.join(base_path, variable), '.grib') if str(year) in file]
    outDir = f'{storPath}{variable}/{year}/'
    os.makedirs(outDir, exist_ok=True)

    for file in files:
        month = INT_TO_MONTH[f'{int(file.split('_')[-1].split('.')[0]):02d}']
        warp_ERA5_to_reference(file,
                               dummyTIFF,
                               f'{outDir}{variable}_{month}.tif')