import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
import geopandas as gpd
from collections import defaultdict
from joblib import Parallel, delayed
from datetime import datetime, timedelta
import time


# state the compositing scheme, the folder, where evap files were stored and number of cores to be used during parellel running parts
comp = 'minVZA'
pathi = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/SharpEvap/Brandenburg/SecondShot/evap/'
ncores = 6

# create sub-folder
evap_outFolder = f"{pathi}{comp}/"
vrt_folder= path_safe(f"{evap_outFolder}vrt/")

# separate for 4 different evap products (soil vs canopy, ssrd calc vs ssrd func)
all_files = getFilelist(evap_outFolder, '.tif', deep=True)

soil_func_files = [file for file in all_files if 'Soil_func' in file]
soil_calc_files = [file for file in all_files if 'Soil_calc' in file]

canopy_func_files = [file for file in all_files if 'Canopy_func' in file]
canopy_calc_files = [file for file in all_files if 'Canopy_calc' in file]

# get the max extent for vrts. this is based on the mask (which actually was dervied from a vrt created and then deleted as there are some with smaller extent)
mask_ds = gdal.Open('/data/Aldhani/eoagritwin/fields/Auxiliary/grid_search/Brandenburg/quick_n_dirty/Fields_as_mask_pixel_20.tif')
gt = mask_ds.GetGeoTransform()
mask_proj = mask_ds.GetProjection()
xmin = gt[0]
ymax = gt[3]
px_size_x = gt[1]
px_size_y = abs(gt[5])
xres = mask_ds.RasterXSize
yres = mask_ds.RasterYSize
xmax = xmin + xres * px_size_x
ymin = ymax - yres * px_size_y

vrt_options = gdal.BuildVRTOptions(
    outputBounds=[xmin, ymin, xmax, ymax],
    xRes=px_size_x,       
    yRes=px_size_y, 
    outputSRS=mask_proj,
    resampleAlg='nearest',
    separate=False
)

months = [REAL_INT_TO_MONTH[m] for m in range(1,13) if growingSeasonChecker(m)]

date_reg = re.compile(r'20\d{2}_(?:' + '|'.join(months) + r')_\d{1,2}_')
date_reg2 = re.compile(r'20\d{2}_(?:' + '|'.join(months) + r')_\d{1,2}.')

filesList = [soil_func_files, soil_calc_files, canopy_func_files, canopy_calc_files]
filesNames = ['Soil_func', 'Soil_calc', 'Canopy_func', 'Canopy_calc']

# create daily vrts mosaics for each type of evap output (soil, canopy)
for fname, fileL in zip(filesNames, filesList):
    dicti = defaultdict(list)
    for file in fileL:
        match = date_reg.search(file)
        if match:
            date = match.group()[:-1]  
            year, month, day = date.split('_')
            day = day.zfill(2)  # "1" → "01", "13" → "13"
            date_key = f"{year}_{month}_{day}"  
            dicti[date_key].append(file)
            continue
    
    for key, files in dicti.items():
        
        outfolder = f'{vrt_folder}{fname}/'

        vrt_path = path_safe(f'{outfolder}{fname}_{key}.vrt')
        vrt = gdal.BuildVRT(vrt_path, files, options=vrt_options)
        vrt = None

        convertVRTpathsTOrelative(vrt_path=vrt_path)
        pass

# load the comp dates
compDates = pd.read_csv(getFilelist(pathi, '.csv', deep=True)[0])['compdates'].tolist()

starts = [datetime.strptime(str(compDate), '%Y%m%d') - timedelta(days=4) for compDate in compDates]
ends = [datetime.strptime(str(compDate), '%Y%m%d') + timedelta(days=4) for compDate in compDates]

# check lowest and highest start date in vrts and adapt start and end if needed
vrt_dates = list(dicti.keys())
vrt_datesL = []
for vrt_date in vrt_dates:
    year, month, day = vrt_date.split('_') 
    vrt_datesL.append(f'{year}{MONTH_TO_02D[month]}{day}')
vrt_start = datetime.strptime(str(min(vrt_datesL)), '%Y%m%d')
vrt_end = datetime.strptime(str(max(vrt_datesL)), '%Y%m%d')

if vrt_start < starts[0]:
    starts[0] = vrt_start

if vrt_end > ends[-1]:
    ends[-1] = vrt_end

# load mask for vrts (based on the polygonized delineated fields (@20m!!!!))
mask_arr = mask_ds.GetRasterBand(1).ReadAsArray()
mask_arr[mask_arr>0] = 1

outpath_single = path_safe(f"{evap_outFolder}median_9day/ET_single/")
outpath_sum = path_safe(f"{evap_outFolder}median_9day/ET_sum/")

joblist = []

for et_var in ['Soil_func', 'Soil_calc', 'Canopy_func', 'Canopy_calc']:
    for date_start, date_end in zip(starts, ends):
        joblist.append([vrt_folder, et_var, date_start, date_end, date_reg2, outpath_single, mask_arr])
print(f'\n{len(joblist)} jobs will be processed\n')


if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time:" + starttime)
    print("")

    Parallel(n_jobs=ncores)(delayed(make_ET_median)(job[0], job[1], job[2], job[3], job[4], job[5], job[6]) for job in joblist)

    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start : " + starttime)
    print("end: " + endtime)
    print("")


print('all median composites calculated - now canopy and soil will be sumed up')

medList = getFilelist(outpath_single, '.tif')
for date_start in starts:

    dt = datetime.strftime(date_start, '%Y_%m_%d')
    date_sub_calc = [med_arr for med_arr in medList if dt in med_arr and 'calc' in med_arr]
    date_sub_func = [med_arr for med_arr in medList if dt in med_arr and 'func' in med_arr]

    arrL = []
    for pathi in date_sub_calc:
        ds = gdal.Open(pathi)
        arrL.append(ds.GetRasterBand(1).ReadAsArray())
    arr_sum = np.nansum(np.dstack([arrL[0], arrL[1]]),axis=2)
    arr_sum[arr_sum >= 10] = np.nan
    # we cut off the highest values at the Polish border for now
    cutOff = np.nanpercentile(arr_sum, [99.9])[0]
    arr_sum[arr_sum > cutOff] = np.nan
    makeTif_np_to_matching_tif(arr_sum, pathi, f"{outpath_sum}{dt}_median_ET.tif", 0)