import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.guzinski import *
from helperToolz.dicts_and_lists import INT_TO_MONTH

# set storPath for exported tiffs
LST_stor_Path = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/minVZAcomp/'
LST_path = '/data/Aldhani/eoagritwin/et/Sentinel3/raw_LST/'
VZA_path = '/data/Aldhani/eoagritwin/et/Sentinel3/VZA/monthly_tiff_values/'
AcqTime_stor_path = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/int_format/'
AirTemp_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/low_res/2m_temperature/'

for year in [2019]:
    LST_stor_Path_year = f'{LST_stor_Path}{year}/' 
    AcqTime_stor_path_year = f'{AcqTime_stor_path}{year}/'

    [os.makedirs(path, exist_ok=True) for path in [LST_stor_Path_year, AcqTime_stor_path_year]]

    # get a subset of LST and VZA files for that year
    files = sorted(getFilelist(LST_path, '.nc'))
    yearFiles_LST = [file for file in sorted(getFilelist(LST_path, '.nc')) if int(file.split('/')[-1].split('_')[-1][0:4]) == year]
    yearFiles_VZA = getFilelist(f'{VZA_path}{year}/', 'tif')
    yearFiles_2mT = getFilelist(f'{AirTemp_path}{year}', '.tif')

    # loop over files and export to .tif at Path locations
    for month in [f'{i:02d}' for i in range(1,13)]:
        if growingSeasonChecker(int(month)):
            
            # subset LST to month and get acquisition time and calculate observations per day
            file_LST = [yearfile_LST for yearfile_LST in yearFiles_LST if f'{month}.nc' == yearfile_LST.split('-')[-1]][0]
            accDateTimes = getAllDatesS3(file_LST) 
            df = pd.Series(accDateTimes)
            counts_per_day = df.dt.floor("D").value_counts().sort_index()
            # make iterables from counts per day that catch starting and ending indices to subset all obs per day
            cumulative_day_counts_end = np.asarray(np.cumsum(counts_per_day))
            cumulative_day_counts_start = np.insert(cumulative_day_counts_end, 0 ,0)

            # load data
            dat_LST = getDataFromNC_LST(file_LST)
            
            # apply the temperature threshold
            dat_LST[dat_LST<273.15] = np.nan # LST_MASKING check!

            # get VZA stack and 2m airtemperature mask
            file_VZA = [yearfile_VZA for yearfile_VZA in yearFiles_VZA if f'{month}.tif' == yearfile_VZA.split('_')[-1]][0]
            dat_VZA = stackReader(file_VZA)

            # sanity check
            if (dat_LST.shape == dat_VZA.shape):
                
                # check air temperature (2m ERA5)
                file_2mT = [yearFile_2mT for yearFile_2mT in yearFiles_2mT if f'{INT_TO_MONTH[month]}.tif' == yearFile_2mT.split('_')[-1]][0]
                dat_2mT, time_2mT = stackReader(file_2mT, bands=True)
                
                #### get ERA5 AirTemp (interpolated from both modelled values that are closest to LST)
                bands_low = []
                minutes = [] # get the minutes to interpolate ERA5 temp values to the exact minute of LST acquisition
                for accDT in accDateTimes: # search for each LST observation 
                    for count, air_time in enumerate(time_2mT): # the two neighbouting ERA5 air temp values
                        if accDT.astype('datetime64[h]')== pd.Timestamp(air_time): # this will get the hourly value before the acquisition
                            bands_low.append(count)
                            minutes.append(pd.Timestamp(accDT).minute)
                bands_up = [band + 1 for band in bands_low]# this get the hourly value after the acquisition

                # interpolate
                air_temp_intpol = dat_2mT[:,:,bands_low] - (dat_2mT[:,:,bands_low] - dat_2mT[:,:,bands_up]) * (np.array(minutes, dtype=np.float32) / 60).reshape(1,1,-1)

                # apply air threshold
                dat_LST = np.where((dat_LST - air_temp_intpol) < -2, np.nan, dat_LST)
    
                # now get VZA minimum composite
                minVZAL = []
                minACQL = [] # collect acquisition times of minVZA pixel
                minACQL_read = [] # collect readable acquisition times of minVZA pixel 
                doyL = [] # for band names when exporting
                
                for l in range(len(counts_per_day)):

                    ################## LST values
                    # Select the slices for the day:
                    LST_slice = dat_LST[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]  # shape (X,Y,Z)
                    VZA_slice = dat_VZA[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]  # shape (X,Y,Z)

                    # Create mask where LST is valid and VZA < 45
                    valid_mask = (~np.isnan(LST_slice)) & (VZA_slice < 45)

                    # For each (x,y), set VZA invalid points to a large number so they don't become min
                    vza_for_min = np.where(valid_mask, VZA_slice, np.inf)  # shape (X,Y,Z)

                    # Find index of minimal VZA along axis=2 (time/bands) for each pixel
                    min_vza_idx = np.argmin(vza_for_min, axis=2)  # shape (X,Y)

                    # Now use advanced indexing to get the corresponding LST values:
                    x_indices = np.arange(LST_slice.shape[0])[:, None]  # shape (X,1)
                    y_indices = np.arange(LST_slice.shape[1])[None, :]  # shape (1,Y)

                    best_LST = LST_slice[x_indices, y_indices, min_vza_idx]  # shape (X,Y)

                    # Take care of all invalid pixel that might have sneaked in through np.argmin
                    no_valid_points = ~np.any(valid_mask, axis=2)  # shape (X,Y)
                    best_LST[no_valid_points] = np.nan
                    minVZAL.append(best_LST)
                    doyL.append(f'DOY_{l+1}')

                    ################# Time of observation of selected pixel --> needed for ERA5 stuff
                    time_slice = df[cumulative_day_counts_start[l]:cumulative_day_counts_end[l]].values
                    timestamp_array = np.tile(time_slice, dat_LST.shape[:2] + (1,))
                    acq_time = timestamp_array[x_indices, y_indices, min_vza_idx]  
                    acq_time_unix = acq_time.astype('datetime64[s]').astype(int) # convert back with pd.to_datetime(best_time_unix, unit='s')
                    acq_time_unix[no_valid_points] = 0 # use 0 as na for export
                    minACQL.append(acq_time_unix)
                    # and also as readable tiffs
                    datetimes = time_slice.astype('datetime64[m]').astype('O')
                    time_arr = np.array([int(dt.strftime("%H%M")) for dt in datetimes])
                    timestamp_array = np.tile(time_arr, dat_LST.shape[:2] + (1,))
                    acq_time = timestamp_array[x_indices, y_indices, min_vza_idx]  
                    acq_time[no_valid_points] = 0 # use 0 as na for export
                    minACQL_read.append(acq_time)

                # export minVZA LST composite
                # exportNCarrayDerivatesInt(file_LST, LST_stor_Path_year, f'Daily_LST_VZAmincomp_{year}_{INT_TO_MONTH[month]}.tif',
                #                           doyL, np.dstack(minVZAL), numberOfBands=len(minVZAL))
                exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_VZAmincomp_{year}_{INT_TO_MONTH[month]}.tif',
                                        doyL, np.dstack(minACQL), numberOfBands=len(minACQL), noData=0)
                exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_VZAmincomp_{year}_{INT_TO_MONTH[month]}_readable.tif',
                                        doyL, np.dstack(minACQL_read), numberOfBands=len(minACQL), noData=0)
            else:
                raise ValueError(f'S3 LST stack differs from VZA stack. Something is seriously wrong\nLST:{dat_LST.shape} vs VZA{dat_VZA.shape}')