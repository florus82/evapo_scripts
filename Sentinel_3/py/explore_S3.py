import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
# import packages 

int_to_Month = {
    '01': 'January',
    '02': 'February',
    '03': 'March',
    '04': 'April',
    '05': 'May',
    '06': 'June',
    '07': 'July',
    '08': 'August',
    '09': 'September',
    '10': 'October',
    '11': 'November',
    '12': 'December'
    }

files = sorted(getFilelist('/data/Aldhani/eoagritwin/et/Sentinel3/raw', '.nc'))

for year in [i for i in range(2017,2025,1)]:

    # get a subset of files for that year
    yearFiles = [file for file in files if int(file.split('/')[-1].split('_')[-1][0:4]) == year]

    # create a maks for germany
    mask = makeGermanyMaskforNC('/data/Aldhani/eoagritwin/misc/gadm41_DEU_shp/gadm41_DEU_0.shp', yearFiles[0])

    # set storPath for exported tiffs
    storPath = '/data/Aldhani/eoagritwin/et/Sentinel3/tiffs/'
    LST_path = f'{storPath}LST/daily_observations_all/{year}/'
    monthly_composites_path = f'{storPath}LST/monthly_composites/{year}/'
    os.makedirs(LST_path, exist_ok=True)
    os.makedirs(monthly_composites_path, exist_ok=True)
    yearCont = []
    # loop over files and export to .tif at location storPath
    for i, file in enumerate(yearFiles):
        print(file)
        accDateTimes = getAllDatesS3(file) # possible to take annual subset if entire files list would be passed here
    #     convertNCtoTIF(file, LST_path, file.split('/')[-1].split('.')[0] + '.tif', accDateTimes, False, True)

        dat = getDataFromNC(file)
        monthCont = []
        dailyCont = []
        bnames = []
        df = pd.Series(accDateTimes)
        counts_per_day = df.dt.floor("D").value_counts().sort_index()
        # vectors for indexing over days
        cumulative_day_counts_end = np.asarray(np.cumsum(counts_per_day))
        cumulative_day_counts_start = np.insert(cumulative_day_counts_end, 0 ,0)

        # cumulative_day_counts_start = np.array(cumulative_day_counts_start)
        # cumulative_day_counts_end = np.array(cumulative_day_counts_end)

        # aggreagate by median
        # stack_list = [
        #     np.nanmedian(dat[:, :, start:end], axis=2)
        #     for start, end in zip(cumulative_day_counts_start[:-1], cumulative_day_counts_end)
        # ] 
        # fin_block = np.dstack(stack_list)

        MM = int_to_Month[file.rsplit('-', maxsplit=1)[-1].split('.')[0]]
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

        for l in range(len(counts_per_day)):
            # monthCont.append(np.any(~np.isnan(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]),axis=2)) # minimum dail obs
            dailyCont.append(np.sum(~np.isnan(dat[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]),axis=2))
            bnames.append(str(counts_per_day.index[l].date()))
        exportNCarrayDerivatesInt(file, storPath + 'Analytics/Count_obs_per_day/', f'Daily_obs_for_{year}_{MM}.tif', bnames, np.dstack(dailyCont), True, numberOfBands=len(dailyCont))
        # exportNCarrayDerivatesInt(file, storPath + 'Analytics/', f'Minimum_DailyObservations_{('_').join(file.split('_')[-1].split('-')[:2])}.tif', 'monthly_sum_of_daily_obs', np.nansum(np.dstack((monthCont)), axis = 2), True)
        # yearCont.append(np.nansum(np.dstack((monthCont)), axis = 2))
    # exportNCarrayDerivatesInt(file, storPath + 'Analytics/', f'Minimum_DailyObservations_{file.split('_')[-1].split('-')[0]}.tif', 'annual_sum_of_daily_obs', np.nansum(np.dstack((yearCont)), axis = 2), True)