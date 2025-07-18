import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import INT_TO_MONTH

files = sorted(getFilelist('/data/Aldhani/eoagritwin/et/Sentinel3/raw', '.nc'))

for year in [i for i in range(2017,2025,1)]:

    print(f'Start processing .nc files for the year {year}')

    # get a subset of files for that year
    yearFiles = [file for file in files if int(file.split('/')[-1].split('_')[-1][0:4]) == year]

    # create a maks for germany
    mask = makeGermanyMaskforNC('/data/Aldhani/eoagritwin/misc/gadm41_DEU_shp/gadm41_DEU_0.shp', yearFiles[0])

    # set storPath for exported tiffs
    storPath = '/data/Aldhani/eoagritwin/et/Sentinel3/tiffs/'
    LST_path = f'{storPath}LST/daily_observations_all/{year}/'
    Time_path = f'{storPath}LST/Acq_time/int_format/{year}/'
    monthly_composites_path = f'{storPath}LST/monthly_composites/{year}/'
    # ensure that directories exist
    [os.makedirs(dir_path, exist_ok=True) for dir_path in [LST_path, Time_path, monthly_composites_path]]

    yearCont = []# for collecting number of observations per year

    # loop over files and export to .tif at Path locations
    for i, file in enumerate(yearFiles):
        
        print(f'start on file {file}')
        
        accDateTimes = getAllDatesS3(file) # possible to take annual subset if entire files list would be passed here
    #     convertNCtoTIF(file, LST_path, file.split('/')[-1].split('.')[0] + '.tif', accDateTimes, False, True)


        dat = getDataFromNC_LST(file)
        monthCont = [] # for collecting number of observations per month
        dailyCont = [] # for collecting number of observations per day
        dailyVals_median = [] # for collection the actual daily LST values (daily median)
        dailyVals_mean = [] # for collection the actual daily LST values (daily mean)
        dailyVals_max = [] # for collection the actual daily LST values (daily mean)
        bnames = []
        df = pd.Series(accDateTimes)
        counts_per_day = df.dt.floor("D").value_counts().sort_index()
        # vectors for indexing over days
        cumulative_day_counts_end = np.asarray(np.cumsum(counts_per_day))
        cumulative_day_counts_start = np.insert(cumulative_day_counts_end, 0 ,0)

        # cumulative_day_counts_start = np.array(cumulative_day_counts_start)
        # cumulative_day_counts_end = np.array(cumulative_day_counts_end)

        # get accquisition times ready for aggregation and export
        time_cube = np.full(dat.shape, np.nan)#, dtype='float64')
        for i in range(len(df)):
            mask = ~np.isnan(dat[:,:,i])
            time_cube[:,:,i][mask] = convertTimestamp_to_INT(df[i])# df[i].timestamp() # needed conversion as tif export won't work with datetimeobject
        # time_cube = np.where(time_cube == None, 'None', time_cube) # only needed to have dtype=str after running convertTimestamp_to_STR
        dailyTimeDates_max = []
        ################################################ gives monthly min, max, median composites
        # aggreagate by median
        # stack_list = [
        #     np.nanmedian(dat[:, :, start:end], axis=2)
        #     for start, end in zip(cumulative_day_counts_start[:-1], cumulative_day_counts_end)
        # ] 
        # fin_block = np.dstack(stack_list)

        MM = INT_TO_MONTH[file.rsplit('-', maxsplit=1)[-1].split('.')[0]]
        # bands = [f'{MM}_Day_{b+1}' for b in range(fin_block.shape[2])]
        # fin_block = fin_block * mask[:, :, np.newaxis]
        # fin_block[fin_block == 0] = np.nan
        # exportNCarrayDerivatesComp(file, monthly_composites_path, f'Germany_{year}_{MM}_mean.tif', bands, fin_block)

        # # aggreagate by min
        # stack_list = [
        #     np.nanmin(dat[:, :, start:end], axis=2)
        #     for start, end in zip(cumulative_day_counts_start[:-1], cumulative_day_counts_end)
        # ] 
        # fin_block = np.dstack(stack_list)
        # fin_block = fin_block * mask[:, :, np.newaxis]
        # fin_block[fin_block == 0] = np.nan
        # exportNCarrayDerivatesComp(file, monthly_composites_path, f'Germany_{year}_{MM}_min.tif', bands, fin_block)

        # # aggreagate by max
        # stack_list = [
        #     np.nanmax(dat[:, :, start:end], axis=2)
        #     for start, end in zip(cumulative_day_counts_start[:-1], cumulative_day_counts_end)
        # ] 
        # fin_block = np.dstack(stack_list)
        # fin_block = fin_block * mask[:, :, np.newaxis]
        # fin_block[fin_block == 0] = np.nan
        # exportNCarrayDerivatesComp(file, monthly_composites_path, f'Germany_{year}_{MM}_max.tif', bands, fin_block)

        ########################################################### creates metadata raster
        for l in range(len(counts_per_day)):
            # number of observations per month
            monthCont.append(np.any(~np.isnan(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]),axis=2)) # minimum dail obs
            
            # collect the dates to use as bandnames for exported tif stacks
            bnames.append(str(counts_per_day.index[l].date()))

            # collect number of observations per day ( count only one per day!)
            dailyCont.append(np.sum(~np.isnan(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]),axis=2))
            
            # # collect actual LST values
            # dailyVals_median.append(np.nanmedian(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]], axis = 2))
            # dailyVals_mean.append(np.nanmean(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]], axis = 2))
            dailyVals_max.append(np.nanmax(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]], axis = 2))

            # collect the acquisition time and dates in timestamp format
            dailyTimeDates_max.append(getValsatMaxIndex(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]],
                                                        time_cube[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]))
   
        # # export daily values
        # exportNCarrayDerivatesInt(file, LST_path, f'Daily_LST_means_{year}_{MM}.tif', bnames, np.dstack(dailyVals_mean), make_uint16=False, numberOfBands=len(dailyVals_mean))
        # exportNCarrayDerivatesInt(file, LST_path, f'Daily_LST_medians_{year}_{MM}.tif', bnames, np.dstack(dailyVals_median), make_uint16=False, numberOfBands=len(dailyVals_median))
        # exportNCarrayDerivatesInt(file, LST_path, f'Daily_LST_max_{year}_{MM}.tif', bnames, np.dstack(dailyVals_max), make_uint16=False, numberOfBands=len(dailyVals_max))
        exportNCarrayDerivatesInt(file, Time_path, f'Daily_Time_max_{year}_{MM}.tif', bnames, np.dstack(dailyTimeDates_max), make_uint16=False, numberOfBands=len(dailyTimeDates_max))
        # # export day counts
        # exportNCarrayDerivatesInt(file, storPath + 'Analytics/Count_obs_per_day/', f'Daily_obs_for_{year}_{MM}.tif', bnames, np.dstack(dailyCont), True, numberOfBands=len(dailyCont))
        # # export month counts
        # exportNCarrayDerivatesInt(file, storPath + 'Analytics/Count_obs_per_month/', f'Monthly_Min_DailyObs_{('_').join(file.split('_')[-1].split('-')[:2])}.tif', 'monthly_sum_of_daily_obs', np.nansum(np.dstack((monthCont)), axis = 2), True)
        
        # # collect number of observations per year ( count only one per day!)
        # yearCont.append(np.nansum(np.dstack((monthCont)), axis = 2))

    # # export year counts
    # exportNCarrayDerivatesInt(file, storPath + 'Analytics/Count_obs_per_year/', f'Annual_Min_DailyObs_{file.split('_')[-1].split('-')[0]}.tif', 'annual_sum_of_daily_obs', np.nansum(np.dstack((yearCont)), axis = 2), True)