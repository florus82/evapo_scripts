import sys
import rasterio
from datetime import datetime
from pvlib import solarposition
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *

tiles = [file.split('SLOPE_')[-1].split('.')[0] for file in getFilelist('/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/', '.tif')]

# here we don't need to care about if we actually use the related LST value to the time of observation. if we dont use it later, the solar value wont be used

for year in [2019]:# range(2017, 2025, 1):
    for tile in tiles:

        # File paths for slope and aspect rasters (in degrees)
        slope_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/SLOPE_{tile}.tif'
        aspect_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/ASPECT_{tile}.tif'
        dem_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/DEM/DEM_{tile}.tif'
        lat_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LAT/Latitude_{tile}.tif'
        lon_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LON/Longitude_{tile}.tif'
        acq_time_list = getFilelist(f'/data/Aldhani/eoagritwin/et/Sentinel3/tiffs/Acq_time/{year}', '.tif')
        stor_dir = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE/{tile}/{year}/'

        os.makedirs(stor_dir, exist_ok=True)

        for file in acq_time_list:
            month = file.split('.tif')[0].split('_')[-1]
            if month in ['April', 'May', 'June', 'July', 'August', 'September', 'October']:
                # load all data and convert if needed
                with rasterio.open(slope_path) as slope_src:
                    slope = slope_src.read(1)  # Read first band
                    
                with rasterio.open(aspect_path) as aspect_src:
                    aspect = aspect_src.read(1)

                with rasterio.open(lon_path) as lon_src:
                    lon = lon_src.read(1)

                with rasterio.open(lat_path) as lat_src:
                    lat = lat_src.read(1)

                with rasterio.open(dem_path) as dem_src:
                    dem = dem_src.read(1)

                # Replace no data or negative values with nan if needed
                slope = np.where(slope < 0, np.nan, slope)
                aspect = np.where(aspect < 0, np.nan, aspect)

                # Convert degrees to radians for trigonometric calculations
                slope_rad = np.deg2rad(slope)
                aspect_rad = np.deg2rad(aspect)

                # warp S3 dates into tile and read-in
                warped_ds = warp_raster_to_reference(file, reference_path=slope_path, output_path='MEM', resampling='near')
                days = warped_ds.RasterCount

                for day in range(days):
                    if os.path.exists(f'{stor_dir}INCIDENCE_{tile}_{year}_{month}_{(day+1):02d}.tif'):
                        t = time.localtime()
                        ti = time.strftime("%H:%M:%S", t)
                        print(f"already exists - next one at {ti}")
                    else:
                        try:
                            
                            print(f'Working on tile {tile} on {(day+1):02d}/{month}/{year}')

                            time_warp = warped_ds.GetRasterBand(day+1).ReadAsArray()

                            # Flatten Unix time and convert to datetime
                            timestamps_flat = pd.to_datetime(time_warp.ravel(), unit='s', utc=True)

                            # Flatten lat/lon and DEM
                            lat_flat = lat.ravel()
                            lon_flat = lon.ravel()
                            dem_flat = dem.ravel()
                            # Compute solar position
                            solpos = solarposition.get_solarposition(time=timestamps_flat, latitude=lat_flat, longitude=lon_flat, altitude=dem_flat)

                            # Convert degrees to radians and reshape to original 2D
                            zenith_rad = np.deg2rad(solpos['zenith'].values).reshape(time_warp.shape)
                            azimuth_rad = np.deg2rad(solpos['azimuth'].values).reshape(time_warp.shape)

                            # cos(theta_i) = cos(theta_z)*cos(beta) + sin(theta_z)*sin(beta)*cos(gamma_s - gamma)

                            cos_theta_i = (np.cos(zenith_rad) * np.cos(slope_rad) +
                                        np.sin(zenith_rad) * np.sin(slope_rad) * 
                                        np.cos(azimuth_rad - aspect_rad))

                            # convert incidence angle in degrees
                            incidence_angle = np.rad2deg(np.arccos(cos_theta_i))

                            ds = gdal.Open(slope_path)
                            gt = ds.GetGeoTransform()
                            prj = ds.GetProjection()

                            export_intermediate_products('0_0', incidence_angle, gt, prj, stor_dir, f'INCIDENCE_{tile}_{year}_{month}_{(day+1):02d}.tif',typ='float')
                        
                        except Exception as e:
                                print(e)
                                t = time.localtime()
                                ti = time.strftime("%H:%M:%S", t)
                                print(f"thrown at {ti}")
                                continue
            else:
                print(f'{month} is not needed at the moment')