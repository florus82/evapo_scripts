import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.evapo import *
import re

workhorse = True
clean_up = True
clean_up_deep = False
do_S3 = True

if workhorse:
    origin = 'Aldhani/eoagritwin/'
else:
    origin = ''

indices = ['EVI', 'NDM', 'NDV', 'TCB', 'TCG', 'TCW']
metrics = ['mean', 'median']

Month_to_DOY_start_reg = {'January': 1, 'February': 32, 'March': 60, 'April': 91, 'May': 121, 'June': 152, 'July': 182,
    'August': 213, 'September': 244, 'October': 274, 'November': 305, 'December': 335}

Month_to_DOY_start_gap = {'January': 1, 'February': 32, 'March': 61, 'April': 92, 'May': 122, 'June': 153, 'July': 183,
    'August': 214, 'September': 245, 'October': 275, 'November': 306, 'December': 333}

Month_to_Int = {'January': 1, 'February': 2,  'March': 3,  'April': 4, 'May': 5,  'June': 6,  'July': 7, 
    'August': 8, 'September': 9,  'October': 10,  'November': 11, 'December': 12}

def extract_month(file_path):
    match = re.search(r'(S3|S2)_\d{4}_(\w+)_(mean|median)', file_path)
    if match:
        return match.group(2)
    else:
        return None



######################################################################### get the indices
# read csv for valid row_cols for samples to draw. They are based on the share of agriculture (HR Landcover maps) within a S3 pixel
thresh_csv = pd.read_csv(f'/data/{origin}et/Auxiliary/landcover/csv/row_cols.csv')
rows_cols = thresh_csv['Thresh50'].str.split('_', expand=True).astype(int)
ind_rows, ind_cols = rows_cols[0].tolist(), rows_cols[1].tolist()


if do_S3:
    ########################################################################## Sentinel 3
    S3_path = f'/data/{origin}et/Sentinel3/tiffs/LST/daily_observations_all/'
    years = [year for year in os.listdir(S3_path) if os.path.isdir(os.path.join(S3_path, year))]
    years.sort()

    for year in years:
        print(f'Extracting S3 data for year {year}')
        # get the right dict --> needed to get the proper DOYs from the monthly files
        if int(year) in range(2000,2040,4):
            conv = Month_to_DOY_start_gap
        else:
            conv = Month_to_DOY_start_reg

        year_files = getFilelist(os.path.join(S3_path, year), '.tif')

        for metric in metrics:
            print(f'....for metric {metric}')
            metric_paths = [file for file in year_files if metric in file]

            for metric_path in metric_paths:

                # get the right doy starting point --> needed because files are stored monthly and doy needed to merge with S2
                for month, doy_start_ind in conv.items():
                    if month in metric_path:
                        doy_start = doy_start_ind
                        monthi = month
                        break
                print(f' work on {monthi}')
                outPath1 = f'/data/{origin}et/Training_ML/training_data/S3_{year}_{month}_{metric}.parquet'

                if os.path.isfile(outPath1):
                    print(f'{outPath1} already exists - next file')
                    continue
                else:
                    # read in stack
                    ds = gdal.Open(metric_path)
                    bands = ds.RasterCount
                    conti = []

                    for band in range(bands):
                        conti.append(ds.GetRasterBand(band+1).ReadAsArray())
                    arr = np.dstack(conti)
                    # write in dicts
                    s3_dicts = []
                    sub_arr = arr[ind_rows, ind_cols,:].T

                    for row in range(sub_arr.shape[0]):
                        for col in range(sub_arr.shape[1]): 
                                s3_dicts.append({
                                    'row': ind_rows[col],
                                    'col': ind_cols[col],
                                    'doy': doy_start + row,
                                    f'S3_{metric}': sub_arr[row,col]
                                })

                    s3_dict = pd.DataFrame(s3_dicts)
                    s3_dict.to_parquet(outPath1, index=False)


########################################################################## Sentinel 2
# get FORCE extracts
npy_files = getFilelist(f'/data/{origin}et/Auxiliary/dumps_for_training_collecting/', '.npy')

for year in range(2017, 2025, 1):
    year_list = [npy_file for npy_file in npy_files if str(year) in npy_file]

    if not year_list:
        print(f'nothing for S2 in year {year} yet')
        continue
    else:
        for index in indices:
            dat = [np.load(file) for file in year_list if index in file]

            if not dat:
                print(f'nothing for S2 in {year} for index {index}')
                continue
            else:
                outPath2 = f'/data/{origin}et/Training_ML/training_data/S2_{year}_{index}.parquet'
                if os.path.isfile(outPath2):
                    print(f'file {outPath2.split('/')[-1]} already exists :) .... not doing it again')
                    continue
                else:
                    print(f'work on file S2_{year}_{index}')
                    npy_dicts = []
                    dat = dat[0]
                    for row in range(dat.shape[0]):
                        for col in range(dat.shape[1]): 
                            npy_dicts.append({'row': ind_rows[col],
                                            'col': ind_cols[col],
                                            'index': index,
                                            'doy': row+1,
                                            'S2': dat[row,col]})
                            
                    npy_dict = pd.DataFrame(npy_dicts)
                    npy_dict.to_parquet(outPath2, index=False)


if clean_up:
    ########################################################################## clean-up Sentinel 3 (make months to single year files)
    S3_files = [file for file in getFilelist(f'/data/{origin}et/Training_ML/training_data/', '.parquet') if 'S3' in file]
    years = sorted(list(set([re.search(r'_(\d{4})_', S3_file).group(1) for S3_file in S3_files])))
    for year in years:
        for metric in metrics:
            # get files for one year and one metric
            annual_files = [S3_file for S3_file in S3_files if year in S3_file and metric in S3_file]
            # sort them chronological
            annual_files = sorted(annual_files, key=lambda f: Month_to_Int.get(extract_month(f), 0))
            # load and merge them
            merged_files = pd.concat([pd.read_parquet(annual_file) for annual_file in annual_files], ignore_index=False)
            # export
            merged_files.to_parquet(f'/data/{origin}et/Training_ML/training_data/S3_{year}_{metric}.parquet', index=False)
            _ = [os.remove(file) for file in annual_files]


if clean_up_deep:
    ########################################################################## clean-up even further (drop metrices)
    S3_files = [file for file in getFilelist(f'/data/{origin}et/Training_ML/training_data/', '.parquet') if 'S3' in file]
    years = sorted(list(set([re.search(r'_(\d{4})_', S3_file).group(1) for S3_file in S3_files])))
    for year in years:
        # get files for one year and one metric
        annual_files = [S3_file for S3_file in S3_files if year in S3_file]
        # load and join them
        df1 = pd.read_parquet(annual_files[0])
        df2 = pd.read_parquet(annual_files[1])
        
        merged_df = pd.merge(df1, df2, on=['row', 'col', 'doy'], how='left')
        # export
        merged_df.to_parquet(f'/data/{origin}et/Training_ML/training_data/S3_{year}.parquet', index=False)
        _ = [os.remove(file) for file in annual_files]