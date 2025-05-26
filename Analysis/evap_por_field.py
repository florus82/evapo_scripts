import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.evapo import *
import geopandas as gpd
import tarfile
from pyproj import CRS
from affine import Affine
from rasterio.enums import Resampling

workhorse = True

if workhorse:
    origin = 'Aldhani/eoagritwin/'
else:
    origin = ''


# load evapo composites
evap_path = f'/data/{origin}et/Landsat/composites/Brandenburg/'
evaps = getFilelist(evap_path, '.tif', deep=True)

# get fields
fields = f'/data/{origin}fields/Auxiliary/grid_search/Brandenburg/quick_n_dirty/Fields_polygons.shp'
shp_ds = ogr.Open(fields)
layer = shp_ds.GetLayer()

# rasterize fields
dummy = gdal.Open(evaps[0])
dummy_x_size = dummy.RasterXSize
dummy_y_size = dummy.RasterYSize
dummy_gt = dummy.GetGeoTransform()
dummy_proj = dummy.GetProjection()

sub = gdal.GetDriverByName('MEM').Create('', dummy_x_size, dummy_y_size, 1, gdal.GDT_UInt32)
sub.SetGeoTransform(dummy_gt)
sub.SetProjection(dummy_proj)
band = sub.GetRasterBand(1)
band.SetNoDataValue(0)
gdal.RasterizeLayer(sub, [1], layer, options=["ATTRIBUTE=FieldID"])
sub_arr = sub.ReadAsArray()


# load in composites and stack em
ds_list = [gdal.Open(evap) for evap in evaps]
raster_list = [ds.GetRasterBand(1).ReadAsArray() for ds in ds_list]
block = np.dstack(raster_list)


years = [evap.split('Brandenburg')[1].split('/')[1] for evap in evaps]
months = [evap.split('_')[-2] for evap in evaps]

zone_ids = np.unique(sub_arr)
zone_ids = zone_ids[zone_ids != 0]  # remove background if 0
results = []

print(f'total zones: {zone_ids}')
for zone_id in zone_ids:

    if zone_id % 1000 == 0:
        print(zone_id)

    mask = sub_arr == zone_id

    for band_idx in range(1, block.shape[2]):

        band_stats = {
            "zone_id": zone_id,
            "Year": years[band_idx],
            "month": months[band_idx]
        }

        band = block[:, :, band_idx]
        values = band[mask]
        values = values[~np.isnan(values)]  # optional: remove NaNs

        if values.size > 0:
            band_stats["band_mean"] = values.mean()
            band_stats["band_std"] = values.std()
            band_stats["band_min"] = values.min()
            band_stats["band_max"] = values.max()
        else:
            band_stats["band_mean"] = np.nan
            band_stats["band_std"] = np.nan
            band_stats["band_min"] = np.nan
            band_stats["band_max"] = np.nan
        
        results.append(band_stats)



df = pd.DataFrame(results)
df.to_csv(f'/data/{origin}et/Landsat/composites/Brandenburg/field_stats.csv', index=False)