import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.guzinski import * 
from helperToolz.dicts_and_lists import INT_TO_MONTH
import cfgrib
import pvlib
from pvlib.irradiance import get_total_irradiance
from pvlib.location import Location

os.environ["PROJ_LIB"] = "/data/Aldhani/users/potzschf/conda/envs/cds_era5/share/proj"


# define era5 paths for loading and exporting
storPath = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/'
base_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/grib/'
# get variables downloaded
dem_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/DEM/'
slope_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/' # probably needed
aspect_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/'
lat_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LAT/'
lon_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LON/'
time_path = # maybe get that rather from the file itself!!!
variables = [name for name in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, name))]
tiles = getFilelist(dem_path, '.tif')


variable = 'surface_solar_radiation_downward_clear_sky'#'2m_temperature'



for year in [2019]:#range(2017, 2025, 1):
    files = [file for file in getFilelist(os.path.join(base_path, variable), '.grib') if str(year) in file]
    outDir = f'{storPath}{variable}/{year}/'
    
    os.makedirs(outDir, exist_ok=True)

    for file in files:
        m = int(file.split('_')[-1].split('.')[0])
        if growingSeasonChecker(m):
            month = INT_TO_MONTH[f'{m:02d}']
            warp_ERA5_to_reference(file,
                                dem_path,
                                f'{outDir}{variable}_{year}_{month}.tif')
    else:
        continue


    from pvlib.irradiance import get_total_irradiance
from pvlib.location import Location

# --- Inputs ---
# ERA5 variable: ssrdcs in J/m^2 over 1 hour (convert to W/m^2)
ssrdcs_j = 3600000  # e.g. 3.6 MJ/m^2 for 1 hour
ssrdcs = ssrdcs_j / 3600  # convert to W/m²

# Time and location
times = pd.date_range('2022-06-21 12:00', periods=1, freq='1H', tz='UTC')
latitude = 37.0
longitude = -3.5
altitude = 1500  # meters

# Terrain
slope_deg = 30     # degrees
aspect_deg = 180   # south-facing

# Create a pvlib location
site = Location(latitude, longitude, tz='UTC', altitude=altitude)
solpos = site.get_solarposition(times)

# Extraterrestrial radiation (horizontal)
dni_extra = pvlib.irradiance.get_extra_radiation(times)

# Compute clearness index
ghi = ssrdcs  # ERA5 clear-sky horizontal global radiation
cos_zenith = np.cos(np.radians(solpos['zenith'].values[0]))
ghi_clear = dni_extra.values[0] * cos_zenith
kt = ghi / ghi_clear

# Decompose GHI to DNI and DHI using Erbs model
dni, dhi = pvlib.irradiance.erbs(ghi, solpos['zenith'], times)['dni'], \
           pvlib.irradiance.erbs(ghi, solpos['zenith'], times)['dhi']

# Now compute radiation on the tilted terrain
tilt = slope_deg
azimuth = aspect_deg  # 180 = south
irradiance_tilted = get_total_irradiance(
    surface_tilt=tilt,
    surface_azimuth=azimuth,
    dni=dni,
    ghi=ghi,
    dhi=dhi,
    solar_zenith=solpos['zenith'],
    solar_azimuth=solpos['azimuth']
)

print("GHI (horizontal):", ghi, "W/m²")
print("Irradiance on slope:")
print("  Beam:", irradiance_tilted['beam'].values[0])
print("  Diffuse:", irradiance_tilted['diffuse'].values[0])
print("  Ground-reflected:", irradiance_tilted['grounddiffuse'].values[0])
print("  Total:", irradiance_tilted['poa_global'].values[0])