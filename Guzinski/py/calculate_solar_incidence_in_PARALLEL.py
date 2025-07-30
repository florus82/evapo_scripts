import sys
import rasterio
from datetime import datetime
import time
from pvlib import solarposition
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from joblib import Parallel, delayed
from shapely.geometry import Point
import geopandas as gpd

def calc_Incidence_per_tile(tile, year):

    # File paths for slope and aspect rasters (in degrees)
    slope_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/SLOPE_{tile}.tif'
    aspect_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/ASPECT_{tile}.tif'
    dem_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/DEM/DEM_{tile}.tif'
    lat_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LAT/Latitude_{tile}.tif'
    lon_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LON/Longitude_{tile}.tif'
    acq_time_list = getFilelist(f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/int_format/{year}', '.tif')
    stor_dir = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE/minVZA/{year}/{tile}/'

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
                    # print(f"already exists - next one at {ti}")
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

                        # export npz dump
                        # np.savez_compressed(f'{stor_dir}INCIDENCE_{tile}_{year}_{month}_{(day+1):02d}.npz', incidence_angle=incidence_angle) # np.load("incidence_angle.npz")['incidence_angle']
                        # df = pd.DataFrame({'incidence_angle': incidence_angle.ravel()})
                        # df.to_parquet(f'{stor_dir}INCIDENCE_{tile}_{year}_{month}_{(day+1):02d}.parquet', index=False)
                        # np.save(f'{stor_dir}SHAPE_{tile}_{year}_{month}_{(day+1):02d}.npy', original_shape)

                        export_intermediate_products('0_0',
                                                     incidence_angle, gt, prj, stor_dir,
                                                     f'INCIDENCE_{tile}_{year}_{month}_{(day+1):02d}.tif',typ='float',
                                                     comp=True)
                        
                        print('export done!')

                    except Exception as e:
                            print(e)
                            t = time.localtime()
                            ti = time.strftime("%H:%M:%S", t)
                            print(f"thrown at {ti}")
                            continue
        else:
            continue
            # print(f'{month} is not needed at the moment')


ncores = 75
year = 2019

tiles = [file.split('SLOPE_')[-1].split('.')[0] for file in getFilelist('/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/', '.tif')]

# filter for tiles that have already been processed completely
# tiles2 = [tile for tile in tiles if len(getFilelist(f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE/minVZA/{year}/{tile}/', '.tif')) != 214]
# jobs = [[tiles2[i], 2019]  for i in range(len(tiles2))]
# print(f'\n{len(tiles2)} tiles will be processed\n')

jobs = [[tiles[i], 2019]  for i in range(len(tiles))]
print(f'\n{len(tiles)} tiles will be processed\n')


###### needed to initiate while to run forever until every files were calculated (needed due to weird multiply overflow error)
bad_tiles = [1]


if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time:" + starttime)
    print("")
    # while len(bad_tiles) > 0 :
    Parallel(n_jobs=ncores)(delayed(calc_Incidence_per_tile)(i[0], i[1]) for i in jobs)
        # base_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE/'
        # tiles = get_forceTSI_output_Tiles(getFilelist(base_path, '.tif', deep=True))
        # bad_tiles = []
        # for tile in tiles:
        #     t_len = len(getFilelist(f'{base_path}{tile}/{year}', '.tif'))
        #     if t_len != 214: # arbitrary value, number of files for April - October
        #         bad_tiles.append(tile)

    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start : " + starttime)
    print("end: " + endtime)
    print("")
