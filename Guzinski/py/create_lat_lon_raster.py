import sys
import rasterio
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from rasterio.transform import xy
from pyproj import Transformer



tiles = [file.split('SLOPE_')[-1].split('.')[0] for file in getFilelist('/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/', '.tif')]

for tile in tiles:
    slope_raster_path = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/SLOPE_{tile}.tif'

    # get metadata
    with rasterio.open(slope_raster_path) as src:
        width = src.width
        height = src.height
        transform = src.transform
        crs_src = src.crs 

    # create grid
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))

    # get center coordinates of pixel
    xs, ys = rasterio.transform.xy(transform, rows, cols, offset='center')
    xs = np.array(xs).reshape((height, width))
    ys = np.array(ys).reshape((height, width))

    # transform to wgs84
    transformer = Transformer.from_crs("EPSG:3035", "EPSG:4326", always_xy=True)
    lons, lats = transformer.transform(xs, ys)

    # export
    out_meta = src.meta.copy()
    out_meta.update({
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:3035",
        "transform": src.transform
    })

    with rasterio.open(f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LON/Longitude_{tile}.tif', 'w', **out_meta) as dst:
        dst.write(lons.astype('float32'), 1)

    # Example to save lat raster
    with rasterio.open(f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LAT/Latitude_{tile}.tif', 'w', **out_meta) as dst:
        dst.write(lats.astype('float32'), 1)