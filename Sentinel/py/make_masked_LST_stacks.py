import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.guzinski import *
from helperToolz.dicts_and_lists import INT_TO_MONTH

# set storPath for exported tiffs
LST_stor_Path_minVZA = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/LST_composites/minVZA/'
LST_stor_Path_maxLST = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/LST_composites/maxLST/'
LST_stor_Path_top3 = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/LST_composites/top3/'

LST_path = '/data/Aldhani/eoagritwin/et/Sentinel3/raw_LST/'
VZA_path = '/data/Aldhani/eoagritwin/et/Sentinel3/VZA/monthly_tiff_values/'

VZA_stor_Path_minVZA = '/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/minVZA/'
VZA_stor_Path_maxLST = '/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/maxLST/'
VZA_stor_Path_top3 = '/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/top3/'

AcqTime_stor_path = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/int_format/'
AirTemp_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/low_res/2m_temperature/'

for year in [2020]:

    # set storage paths to the current year
    LST_stor_Path_minVZA_year = f'{LST_stor_Path_minVZA}{year}/'
    LST_stor_Path_maxLST_year = f'{LST_stor_Path_maxLST}{year}/' 
    LST_stor_Path_top3_year =  f'{LST_stor_Path_top3}{year}/' 

    VZA_stor_Path_minVZA_year = f'{VZA_stor_Path_minVZA}{year}/'
    VZA_stor_Path_maxLST_year = f'{VZA_stor_Path_maxLST}{year}/'
    VZA_stor_Path_top3_year = f'{VZA_stor_Path_top3}{year}/' 

    AcqTime_stor_path_year = f'{AcqTime_stor_path}{year}/'

    # make paths if they don't exist
    [os.makedirs(path, exist_ok=True) for path in [LST_stor_Path_minVZA_year, LST_stor_Path_maxLST_year, LST_stor_Path_top3_year,
                                                    AcqTime_stor_path_year, VZA_stor_Path_maxLST_year, VZA_stor_Path_minVZA_year, VZA_stor_Path_top3_year]]

    # get a temporal subset of LST, VZA and air temp files for that year
    files = sorted(getFilelist(LST_path, '.nc'))
    yearFiles_LST = [file for file in sorted(getFilelist(LST_path, '.nc')) if int(file.split('/')[-1].split('_')[-1][0:4]) == year]
    yearFiles_VZA = getFilelist(f'{VZA_path}{year}/', 'tif')
    yearFiles_2mT = getFilelist(f'{AirTemp_path}{year}', '.tif')
    mask = makeGermanyMaskforNC('/data/Aldhani/eoagritwin/misc/gadm41_DEU_shp/gadm41_DEU_0.shp', yearFiles_LST[0])

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

            # load data (/all observations for that month)
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

                # interpolate to the minute of observation
                air_temp_intpol = dat_2mT[:,:,bands_low] - (dat_2mT[:,:,bands_low] - dat_2mT[:,:,bands_up]) * (np.array(minutes, dtype=np.float32) / 60).reshape(1,1,-1)

                # apply air threshold
                dat_LST = np.where((dat_LST - air_temp_intpol) < -2, np.nan, dat_LST)
    
                # now get composites (minVZA, maxLST)
                minVZA_LST = [] # collects 2D numpy arrays with masked LST values from minVZA compositiing
                maxLST_LST = [] # collects 2D numpy arrays with masked LST values from maxLST compositiing
                minVZA_VZA = [] # collects 2D numpy arrays with masked VZA values from minVZA compositiing
                maxLST_VZA = [] # collects 2D numpy arrays with masked VZA values from maxLST compositiing
                
                minACQL = [] # collect acquisition times of minVZA pixel
                minACQL_read = [] # collect readable acquisition times of minVZA pixel
                maxACQL = [] # collect acquisition times of maxLST pixel
                maxACQL_read = [] # collect readable acquisition times of maxLST pixel 

                lst_order1, lst_order2, lst_order3 = [], [], [] # order1 holds most pixel within scene
                vza_order1, vza_order2, vza_order3 = [], [], []
                order1_time_int, order2_time_int, order3_time_int = [], [], []  
                order1_time_read, order2_time_read, order3_time_read = [], [], []
                doyL = [] # for band names when exporting
                
                for l in range(len(counts_per_day)):

                    ################## LST values
                    # Select the slices for the day:
                    LST_slice = dat_LST[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]  # shape (X,Y,Z)
                    VZA_slice = dat_VZA[:, :, cumulative_day_counts_start[l]:cumulative_day_counts_end[l]]  # shape (X,Y,Z)

                    # Create mask where LST is valid and VZA < 45
                    valid_mask = (~np.isnan(LST_slice)) & (VZA_slice < 45)
                    
                    # For each (x,y), set VZA/LST invalid points to a large number so they don't become min
                    vza_for_minVZA = np.where(valid_mask, VZA_slice, np.inf)  # shape (X,Y,Z)
                    lst_for_maxLST = np.where(valid_mask, LST_slice, -np.inf)  # shape (X,Y,Z)

                    # Find index of minimal VZA/max LST along axis=2 (time/bands) for each pixel
                    minVZA_idx = np.argmin(vza_for_minVZA, axis=2)  # shape (X,Y)
                    maxLST_idx = np.argmax(lst_for_maxLST, axis=2)  # shape (X,Y)

                    # Now use advanced indexing to get the corresponding LST values:
                    x_indices = np.arange(LST_slice.shape[0])[:, None]  # shape (X,1)
                    y_indices = np.arange(LST_slice.shape[1])[None, :]  # shape (1,Y)

                    # fill 
                    best_LST_minVZA = LST_slice[x_indices, y_indices, minVZA_idx]  # shape (X,Y)
                    best_LST_maxLST = LST_slice[x_indices, y_indices, maxLST_idx]  # shape (X,Y)
                    best_VZA_minVZA = VZA_slice[x_indices, y_indices, minVZA_idx]  # shape (X,Y)
                    best_VZA_maxLST = VZA_slice[x_indices, y_indices, maxLST_idx]  # shape (X,Y)

                    # Take care of all invalid pixel that might have sneaked in through np.argmin
                    no_valid_points_minVZA = ~np.any(valid_mask, axis=2)  # shape (X,Y)
                    no_valid_points_maxLST = ~np.any(valid_mask, axis=2)  # shape (X,Y)
                    
                    best_LST_minVZA[no_valid_points_minVZA] = np.nan
                    best_LST_maxLST[no_valid_points_maxLST] = np.nan

                    best_VZA_minVZA[no_valid_points_minVZA] = np.nan
                    best_VZA_maxLST[no_valid_points_maxLST] = np.nan

                    minVZA_LST.append(best_LST_minVZA * mask)
                    maxLST_LST.append(best_LST_maxLST * mask)
                    
                    minVZA_VZA.append(best_VZA_minVZA * mask)
                    maxLST_VZA.append(best_VZA_maxLST * mask)

                    doyL.append(f'DOY_{l+1}')



                    # ################# Time of observation of selected pixel --> needed for ERA5 stuff
                    time_slice = df[cumulative_day_counts_start[l]:cumulative_day_counts_end[l]].values
                    timestamp_array = np.tile(time_slice, dat_LST.shape[:2] + (1,)) # don’t repeat along the time axis, just preserve it (1,) 
                    
                    acq_time = timestamp_array[x_indices, y_indices, minVZA_idx]  
                    acq_time_unix = acq_time.astype('datetime64[s]').astype(int) # convert back with pd.to_datetime(best_time_unix, unit='s')
                    acq_time_unix[no_valid_points_minVZA] = 0 # use 0 as na for export
                    minACQL.append(acq_time_unix * mask)
                    
                    acq_time = timestamp_array[x_indices, y_indices, maxLST_idx]  
                    acq_time_unix = acq_time.astype('datetime64[s]').astype(int) # convert back with pd.to_datetime(best_time_unix, unit='s')
                    acq_time_unix[no_valid_points_maxLST] = 0 # use 0 as na for export
                    maxACQL.append(acq_time_unix * mask)

                    # and also as readable tiffs
                    datetimes = time_slice.astype('datetime64[m]').astype('O')
                    time_arr = np.array([int(dt.strftime("%H%M")) for dt in datetimes])
                    timestamp_array_read = np.tile(time_arr, dat_LST.shape[:2] + (1,))
                    
                    acq_time = timestamp_array_read[x_indices, y_indices, minVZA_idx]  
                    acq_time[no_valid_points_minVZA] = 0 # use 0 as na for export
                    minACQL_read.append(acq_time * mask)

                    acq_time = timestamp_array_read[x_indices, y_indices, maxLST_idx]  
                    acq_time[no_valid_points_maxLST] = 0 # use 0 as na for export
                    maxACQL_read.append(acq_time * mask)


                    # single scene insert
                    lst_order = [lst_order1, lst_order2, lst_order3]
                    vza_order = [vza_order1, vza_order2, vza_order3]
                    order_time_int = [order1_time_int, order2_time_int, order3_time_int]
                    order_time_read = [order1_time_read, order2_time_read, order3_time_read]

                    # search scenes with largest coverage
                    val_pixels, val_counter = [], []
                    for sdx, scence in enumerate(range(valid_mask.shape[2])):
                        val_pixels.append(np.sum(valid_mask[:,:,sdx]))
                        val_counter.append(sdx)
                    snums, val_ind = sortListwithOtherlist(val_pixels, val_counter, rev=True)

                    for ldx, vali in enumerate(val_ind[:3]):
                        lst_order[ldx].append(LST_slice[:,:,vali] * mask * valid_mask[:,:,vali].astype(int))
                        vza_order[ldx].append(VZA_slice[:,:,vali] * mask * valid_mask[:,:,vali].astype(int))

                        acq_time_unix = acq_time = np.where(valid_mask[:,:,vali], timestamp_array[:,:,vali], np.datetime64('1970-01-01T00:00:00')) 
                        order_time_int[ldx].append(acq_time.astype('datetime64[s]').astype(int) * mask)
                        order_time_read[ldx].append(np.where(valid_mask[:,:,vali], timestamp_array_read[:,:,vali], 0) * mask)


                ################## export minVZA LST composite
                # LST masked
                exportNCarrayDerivatesInt(file_LST, LST_stor_Path_minVZA_year, f'Daily_LST_minVZA_{year}_{INT_TO_MONTH[month]}.tif',
                                          doyL, np.dstack(minVZA_LST), numberOfBands=len(minVZA_LST), noData=0)
                # # VZA masked
                exportNCarrayDerivatesInt(file_LST, VZA_stor_Path_minVZA_year, f'Daily_VZA_minVZA_{year}_{INT_TO_MONTH[month]}.tif',
                                          doyL, np.dstack(minVZA_VZA), numberOfBands=len(minVZA_VZA), noData=0)
                # time
                exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_minVZA_{year}_{INT_TO_MONTH[month]}.tif',
                                        doyL, np.dstack(minACQL), datType=gdal.GDT_Int64, numberOfBands=len(minACQL), noData=0)
                exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_minVZA_{year}_{INT_TO_MONTH[month]}_readable.tif',
                                        doyL, np.dstack(minACQL_read), datType=gdal.GDT_Int64, numberOfBands=len(minACQL), noData=0)
                


                ################# export max LST composite
                # LST masked
                exportNCarrayDerivatesInt(file_LST, LST_stor_Path_maxLST_year, f'Daily_LST_maxLST_{year}_{INT_TO_MONTH[month]}.tif',
                                          doyL, np.dstack(maxLST_LST), numberOfBands=len(maxLST_LST), noData=0)
                # VZA masked
                exportNCarrayDerivatesInt(file_LST, VZA_stor_Path_maxLST_year, f'Daily_VZA_maxLST_{year}_{INT_TO_MONTH[month]}.tif',
                                          doyL, np.dstack(maxLST_VZA), numberOfBands=len(maxLST_VZA), noData=0)
                # time
                exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_maxLST_{year}_{INT_TO_MONTH[month]}.tif',
                                        doyL, np.dstack(maxACQL), datType=gdal.GDT_Int64 ,numberOfBands=len(maxACQL), noData=0)
                exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_maxLST_{year}_{INT_TO_MONTH[month]}_readable.tif',
                                        doyL, np.dstack(maxACQL_read), datType=gdal.GDT_Int64 ,numberOfBands=len(maxACQL), noData=0)
                


                ################# export rank "composites"
                for ord in range(3):
                    print(ord)
                    # LST
                    exportNCarrayDerivatesInt(file_LST, LST_stor_Path_top3_year, f'Daily_LST_order{ord + 1}_{year}_{INT_TO_MONTH[month]}.tif',
                                                doyL, np.dstack(lst_order[ord]), numberOfBands=len(lst_order[ord]), noData=0)
                    # VZA 
                    exportNCarrayDerivatesInt(file_LST, VZA_stor_Path_top3_year, f'Daily_VZA_order{ord + 1}_{year}_{INT_TO_MONTH[month]}.tif',
                                                doyL, np.dstack(vza_order[ord]), numberOfBands=len(vza_order[ord]), noData=0)
                    # time
                    exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_order{ord + 1}_{year}_{INT_TO_MONTH[month]}.tif',
                                            doyL, np.dstack(order_time_int[ord]), datType=gdal.GDT_Int64 ,numberOfBands=len(order_time_int[ord]), noData=0)
                    exportNCarrayDerivatesInt(file_LST, AcqTime_stor_path_year, f'Daily_AcqTime_order{ord + 1}_{year}_{INT_TO_MONTH[month]}_readable.tif',
                                            doyL, np.dstack(order_time_read[ord]), datType=gdal.GDT_Int64 ,numberOfBands=len(order_time_read[ord]), noData=0)
            else:
                raise ValueError(f'S3 LST stack differs from VZA stack. Something is seriously wrong\nLST:{dat_LST.shape} vs VZA{dat_VZA.shape}')