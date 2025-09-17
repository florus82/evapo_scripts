import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.dicts_and_lists import INT_TO_MONTH, REAL_INT_TO_MONTH
from helperToolz.guzinski import growingSeasonChecker
from helperToolz.helpsters import *
from helperToolz.evapo import *
import geopandas as gpd
import tarfile
from rasterio.transform import from_bounds
import shutil
import gzip

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

dir_for_corrupt_files = f'/data/{origin}et/Landsat/raw/bad_files/'
os.makedirs(dir_for_corrupt_files, exist_ok=True)

#### get path and rows of scenes that have data for the chosen AOI (e.g. Brandenburg)

# load shapefiles and 
ger = gpd.read_file(f'/data/{origin}misc/gadm41_DEU_shp/gadm41_DEU_1.shp')
states = [state for state in ger['NAME_1']]
for i in [3,0,1,6,7,8,9,10,11,12,13,14,15,2,4,5]: # determine order in which states are processed
    
    state = states[i]
    print(state)
    aoi = ger[ger['NAME_1'] == state]

    orbits = gpd.read_file(f'/data/{origin}misc/WRS2_descending_0/WRS2_descending.shp')
    # check projections
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
    for year in range(2019,2020):#(2024, 2017, -1):
        print(year)
        # outdir = f'/data/{origin}et/Landsat/composites/max/{state}/{year}'
        # if os.path.exists(outdir):
        #     fList = getFilelist(outdir, '.tif')
        #     if len(fList) == 12:
        #         continue
        # else:
        #     os.makedirs(outdir, exist_ok=True)
        
        for month in range(1, 13, 1):
            if growingSeasonChecker(month):
                outdir = f'/data/{origin}et/Landsat/daily_extracts/{state}/{year}/{REAL_INT_TO_MONTH[month]}/'
                os.makedirs(outdir, exist_ok=True) 
                # check if temp_folder is empty and delete everything if not
                tempF = f'/data/{origin}et/Landsat/extracts/'
                if len(getFilelist(tempF, '.nc')) > 0:
                    for end in ['.nc', '.xml', '.txt']:
                        for file in getFilelist(tempF, end):
                            os.remove(file)
                    print('kill complete')

                # subset data for year and month and extract
                year_month = lookUp.get_by_year_and_month(year, month)
                year_month_path_row = [scene for scene in year_month for pr in path_rows if pr in scene]
                 
                for landsat_file in year_month_path_row:
                    try:
                        with tarfile.open(landsat_file, 'r:gz') as tar:
                            tar.extractall(tempF)
                    except (tarfile.ReadError, gzip.BadGzipFile, OSError):
                        basename = os.path.basename(landsat_file)
                        corrupt_path = os.path.join(dir_for_corrupt_files, basename)
                        shutil.move(landsat_file, corrupt_path)


                # List of file paths
                files_nc = getFilelist(f'/data/{origin}et/Landsat/extracts/', '.nc')
                if len(files_nc) == 0:
                    continue
                files_xml = getFilelist(f'/data/{origin}et/Landsat/extracts/', '.xml')
                datasets = []

                bound_coords = {
                    'UTM': [],
                    'UL_X': [],
                    'UL_Y': [],
                    'LR_X': [],
                    'LR_Y': []
                }

                # loop over all extracted files and give spatial ref 
                for f_nc in files_nc:

                    ds = xr.open_dataset(f_nc)
                    da = ds['ETA']
                    if len(np.unique(da.values)) == 1:
                        continue
                    # find UL corner coordinates to find the most outer one for common grid
                    utm_zone, ul_x, ul_y, lr_x, lr_y = get_UTM_zone_and_corners_from_xml(f_nc, files_xml)
                
                    da.rio.set_spatial_dims(x_dim="XDim_ETA", y_dim="YDim_ETA", inplace=True)
                    da.rio.write_crs(f'EPSG:{utm_to_epsg[utm_zone]}', inplace=True)
                    datasets.append(da)

                if len(datasets) == 0:
                    continue

                warped_arrays = []

                for counti, xr_raster in enumerate(datasets):
                    src_ds = xarray_to_gdal_mem(xr_raster)
                    warped_ds = warp_to_template(src_ds, f'/data/{origin}et/Auxiliary/Landsat_GER_mask/states/{state}.tif', outType=gdal.GDT_Float32)
                    # arr = warped_ds.GetRasterBand(1).ReadAsArray()
                    warped_arrays.append(warped_ds)#arr)

                warped_stack = np.dstack(warped_arrays)
                # makeTif_np_to_matching_tif(warped_stack, 
                #                             f'/data/{origin}et/Auxiliary/Landsat_GER_mask/states/{state}.tif',
                #                         f'/data/{origin}et/Landsat/stacks/ETA_{state}_{year}_{month:02d}_median.tif', 0, gdalType=gdal.GDT_Float32, bands=warped_stack.shape[2])
                warped_stack = warped_stack.astype(float)
                warped_stack[warped_stack == 0.0] = np.nan
                
                
                # mask output
                mask_ds = gdal.Open(f'/data/{origin}et/Auxiliary/Landsat_GER_mask/states/{state}.tif')
                mask_arr = mask_ds.GetRasterBand(1).ReadAsArray()
                
                # single day export
                warped_masked = np.where(mask_arr[:, :, None], warped_stack, np.nan)

                for day in range(warped_masked.shape[2]):
                    makeTif_np_to_matching_tif(warped_masked[:, :, day], 
                                                f'/data/{origin}et/Auxiliary/Landsat_GER_mask/states/{state}.tif',
                                                f'{outdir}/Landsat_ETA_{state}_{year}_{month:02d}_{day+1}.tif',
                                                0,
                                                gdalType=gdal.GDT_Float32)

                # monthly composite exports
                # median  = np.nanmax(warped_stack, axis=2)

                # # mask output
                # mask_ds = gdal.Open(f'/data/{origin}et/Auxiliary/Landsat_GER_mask/states/{state}.tif')
                # mask_arr = mask_ds.GetRasterBand(1).ReadAsArray()
                # median_masked = median * mask_arr

                # makeTif_np_to_matching_tif(median_masked, 
                #                             f'/data/{origin}et/Auxiliary/Landsat_GER_mask/states/{state}.tif',
                #                         f'{outdir}/Landsat_ETA_{state}_{year}_{month:02d}_max.tif', 0, gdalType=gdal.GDT_Float32)