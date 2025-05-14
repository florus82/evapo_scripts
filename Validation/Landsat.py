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

utm_to_epsg = {
    '28': 32628,  # Western Portugal, Azores
    '29': 32629,  # Western Spain, Portugal
    '30': 32630,  # Spain, France, UK
    '31': 32631,  # France, Benelux, Germany, Western Norway
    '32': 32632,  # Germany, Denmark, Switzerland, Italy, Austria
    '33': 32633,  # Central Europe: Poland, Czechia, Hungary, Croatia, Sweden, Norway
    '34': 32634,  # Eastern Europe: Finland, Baltic States, Romania
    '35': 32635,  # Western Russia, Ukraine
    '36': 32636,  # Russia, Black Sea region
}

def get_union_bounds(bounds_list):
    """Compute the union of multiple bounding boxes."""
    minx = min([b[0] for b in bounds_list])
    miny = min([b[1] for b in bounds_list])
    maxx = max([b[2] for b in bounds_list])
    maxy = max([b[3] for b in bounds_list])
    return (minx, miny, maxx, maxy)


#### get path and rows of scenes that have data for the chosen AOI (e.g. Brandenburg)
aoi_set_man = 'Brandenburg'

# load shapefiles and check projections
ger = gpd.read_file(f'/data/{origin}misc/gadm41_DEU_shp/gadm41_DEU_1.shp')
aoi = ger[ger['NAME_1'] == aoi_set_man]

orbits = gpd.read_file(f'/data/{origin}misc/WRS2_descending_0/WRS2_descending.shp')

if aoi.crs != orbits.crs:
    aoi = aoi.to_crs(orbits.crs)

# find overlapping paths/rows
#intersecting = orbits[orbits.intersects(aoi.unary_union)]
intersecting = gpd.sjoin(orbits, aoi, how="inner", predicate="intersects")
path_rows = [[p, r] for p, r in zip(intersecting['PATH'], intersecting['ROW'])]
# make sure, paths and rows have the correct format
path_rows = [f'{str(p).zfill(3)}{str(r).zfill(3)}' for p, r in path_rows]

# get all paths from downloaded products --> subsetted to paths and rows
landsat_files = getFilelist(f'/data/{origin}et/Landsat/raw/', '.tar.gz', deep=True)

# create a look-up dictionary for time subsets
lookUp = LandsatETFileManager(landsat_files)

#### do the compositing monthly
for year in range(2020,2025,1):
    for month in range(1, 13, 1):
        print(f'{year}_{month}')
        # check if temp_folder is empty and delete everything if not
        tempF = f'/data/{origin}et/Landsat/extracts/'
        if len(getFilelist(tempF, '.nc')) > 0:
            for file in getFilelist(tempF, '.nc'):
                os.remove(file)
            print('kill complete')

        # subset data and extract
        year_month = lookUp.get_by_year_and_month(year, month)
        year_month_path_row = [scene for scene in year_month for pr in path_rows if pr in scene]

        for landsat_file in year_month_path_row:
            with tarfile.open(landsat_file, 'r:gz') as tar:
                tar.extractall(tempF)

                # List of file paths
        files_nc = getFilelist(f'/data/{origin}et/Landsat/extracts/', '.nc')
        files_xml = getFilelist(f'/data/{origin}et/Landsat/extracts/', '.xml')
        datasets = []

        for f_nc in files_nc:
            
            utm_zone, w, e, n, s = get_UTM_zone_and_corners_from_xml(f_nc, files_xml)
            ds = xr.open_dataset(f_nc)
            da = ds['ETA']
            da.rio.set_spatial_dims(x_dim="XDim_ETA", y_dim="YDim_ETA", inplace=True)
            da.rio.write_crs(f'EPSG:{utm_to_epsg[utm_zone]}', inplace=True)
            # da.rio.write_nodata(-9999, inplace=True)
            datasets.append(da)

        # Reproject to common CRS 
        target_crs = 'EPSG:32633'
        datasets_reprojected = [ds.rio.reproject(target_crs) for ds in datasets]

        # get common bounds
        all_bounds = [ds.rio.bounds() for ds in datasets_reprojected]
        union_bounds = get_union_bounds(all_bounds)

        # define resolution
        res = datasets_reprojected[0].rio.resolution()
        res_x, res_y = abs(res[0]), abs(res[1])

        # Create output grid shape
        minx, miny, maxx, maxy = union_bounds
        width = int((maxx - minx) / res_x)
        height = int((maxy - miny) / res_y)

        # Create a template with correct transform
        transform = Affine.translation(minx, maxy) * Affine.scale(res_x, -res_y)

        # Reproject all datasets to this common grid
        aligned = []
        for da in datasets:
            aligned_da = da.rio.reproject(
                dst_crs=target_crs,
                shape=(height, width),
                transform=transform,
                resampling=Resampling.cubic  # or bilinear/cubic
            )
            aligned.append(aligned_da)

        # Stack and composite (e.g., using nanmean for mosaic)
        stacked = xr.concat(aligned, dim="stack")
        composite = stacked.median(dim="stack", skipna=True)

        # mask composite
        aoi_m = aoi.to_crs(composite.rio.crs)
        composite_clipped = composite.rio.clip(aoi.geometry.values, aoi.crs, drop=True, invert=False)
        composite_clipped.rio.to_raster(f'/data/{origin}et/Landsat/composites/{year}_{month}_median.tif')
