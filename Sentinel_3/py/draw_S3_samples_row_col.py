import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.evapo import *
from datetime import datetime
workhorse = True

if workhorse:
    origin = 'Aldhani/eoagritwin/'
else:
    origin = ''

def get_folders_in_dir(dir):
    return [f for f in os.listdir(dir) if os.path.isdir(os.path.join(dir, f))]

trash = '/data/Aldhani/eoagritwin/et/Auxiliary/trash/band_intermediate/'
temp = '/data/Aldhani/eoagritwin/et/Auxiliary/trash/vrt/'
os.makedirs(trash, exist_ok=True)
os.makedirs(temp, exist_ok=True)

# read csv for valid row_cols for samples to draw. They are based on the share of agriculture (HR Landcover maps) within a S3 pixel
thresh_csv = pd.read_csv(f'/data/{origin}et/Auxiliary/landcover/csv/row_cols.csv')



########################################## get the row_col combinations for thresholds of share of agriculture in pixel
for col in thresh_csv.columns:
    #print(f'finding indices for {col}')
    nested = [entry.split('_') for entry in thresh_csv[col] if type(entry) == str]
    rows, cols = zip(*nested)
    #print(len(rows))
# only use Thresh50, as all valid row_col will be here. For other THresh sets, a subset from this can be taken
# this are the last instances of rows and cols
rows = [int(row) for row in rows]
cols = [int(col) for col in cols]

number_of_indices = 6

############################################# check if all data is there and get the the different stacks (indices)
for year in [2019, 2017]:
  
    force_path = f'/data/Aldhani/eoagritwin/force/output/S3/{year}/'
    folders = get_folders_in_dir(force_path) # folders of all tiles within year
    folder_to_clean = [os.path.join(force_path, folder) for folder in folders if len(getFilelist(os.path.join(force_path, folder), '.tif', deep=True)) != number_of_indices]
    if len(folder_to_clean) > 0: # a folder to clean would be one with more than 6 tif files (6 VIs)
        raise ValueError('the folders do not contain the same number of images...')
    # for folder in folder_to_clean:
    #     [os.remove(file) for file in getFilelist(folder, '.tif', deep=True) if '20181001' in file]
    indices = list(set([file_annual.split('_')[-2] for file_annual in getFilelist(os.path.join(force_path, folders[0]), '.tif', deep=True)]))
    indices.sort() # unique list of indices


############################################### search/define the bands in raster stacks that belong to the year of interest
    ds = gdal.Open(getFilelist(os.path.join(force_path, folders[0]), '.tif', deep=True)[0], 0) # open the first tile-stack
    if ds.RasterCount > 366:
        bandname = str(ds.GetRasterBand(93).GetDescription()) # hard-coded numbers, but works 
        if bandname == f'{year}0101':
            start_index = 93
            end_index = 457
            if year in [2016,2020,2024]:
                end_index +=1
        else:
            raise ValueError('band for Jan 01 not found!')

############################################### get all files for the year
    files = [getFilelist(os.path.join(force_path, folder), '.tif', deep=True) for folder in folders]
    files = [file for list in files for file in list]

    # loop over indices
    
    for index in indices:
        print("Current time:", datetime.now().strftime("%H:%M:%S"),'   ', index)
        
        index_files = [file for file in files if index in file]

        # check if output already created
        outPath = f'/data/{origin}et/Auxiliary/dumps_for_training_collecting/{year}_{index}.npy'
        if os.path.exists(outPath):
            print('file exists already!')
        else:
            conti = []
            for band_number in range(start_index,end_index+1,1):
                # create daily vrts and extract

                for i, tif in enumerate(index_files):
                    output = trash + tif.split('/')[-1].split('.')[0] + '_' + str(band_number) + '_' + str(i) + '.tif'
                    gdal.Translate(output, tif, bandList=[band_number])

                # Now build VRT
                vrt_options = gdal.BuildVRTOptions(separate=False)
                vrt_ds = gdal.BuildVRT(f'{temp}{year}_{index}_{band_number}.vrt', getFilelist(trash, '.tif'), options=vrt_options)
                # vrt_ds.FlushCache()
                # vrt = gdal.Open(f'{temp}{year}_{index}_{band_number}.vrt')
                warpi = warp_to_template(vrt_ds, 
                        '/data/Aldhani/eoagritwin/et/Auxiliary/S3_tif_GER_maks/powermask.tif',
                        mask_path='/data/Aldhani/eoagritwin/et/Auxiliary/S3_tif_GER_maks/powermask.tif',
                        #outPath=f'/data/Aldhani/eoagritwin/et/FORCE/vrt_dumps/{year}_{index}_{band_number}.tif',
                        outType=gdal.GDT_Int16)
                conti.append(warpi[rows,cols])
                [os.remove(file) for file in getFilelist(trash, '.tif')]

            arr = np.stack(conti)
            np.save(outPath, arr)
            
            print("Current time:", datetime.now().strftime("%H:%M:%S"))
