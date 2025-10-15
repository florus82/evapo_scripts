import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
from helperToolz.mirmazloumi import *
from other_repos.pyDMS.pyDMS.pyDMS import *
from other_repos.pyTSEB.pyTSEB import meteo_utils
from other_repos.pyTSEB.pyTSEB import resistances
from other_repos.pyTSEB.pyTSEB import net_radiation
from other_repos.pyTSEB.pyTSEB import clumping_index 
from other_repos.pyTSEB.pyTSEB import TSEB
import time
from joblib import Parallel, delayed


def runSharpi(highResFilename, lowResFilename, lowResMaskFilename, cv, movWin, regrat, outputFilename, useDecisionTree = True):
    commonOpts = {"highResFiles":               [highResFilename],
                    "lowResFiles":              [lowResFilename],
                    "lowResQualityFiles":         [lowResMaskFilename],
                    "lowResGoodQualityFlags":     [1],
                    "cvHomogeneityThreshold":     cv,
                    "movingWindowSize":           movWin,
                    "disaggregatingTemperature":  True}
    dtOpts =     {"perLeafLinearRegression":    True,
                    "linearRegressionExtrapolationRatio": round(regrat, 2)}
    sknnOpts =   {'hidden_layer_sizes':         (10,),
                    'activation':                 'tanh'}
    nnOpts =     {"regressionType":             REG_sklearn_ann,
                    "regressorOpt":               sknnOpts}

    start_time = time.time()
    opts = commonOpts.copy()
    if useDecisionTree:
        opts.update(dtOpts)
        disaggregator = DecisionTreeSharpener(**opts)
    else:
        opts.update(nnOpts)
        disaggregator = NeuralNetworkSharpener(**opts)


    disaggregator.trainSharpener()

    downscaledFile = disaggregator.applySharpener(highResFilename, lowResFilename)

    residualImage, correctedImage = disaggregator.residualAnalysis(downscaledFile, lowResFilename,
                                                                    lowResMaskFilename,
                                                                    doCorrection=True)


    if correctedImage is not None:
        outImage = correctedImage
    else:
        outImage = downscaledFile
    # outData = utils.binomialSmoother(outData)
    outFile = utils.saveImg(outImage.GetRasterBand(1).ReadAsArray(),
                            outImage.GetGeoTransform(),
                            outImage.GetProjection(),
                            f'{os.path.split(outputFilename)[0]}/Values/{os.path.split(outputFilename)[1]}')
    residualFile = utils.saveImg(residualImage.GetRasterBand(1).ReadAsArray(),
                                residualImage.GetGeoTransform(),
                                residualImage.GetProjection(),
                                f'{os.path.split(outputFilename)[0]}/Residuals/{os.path.split(outputFilename)[1]}_resid{os.path.splitext(outputFilename)[1]}')

    outFile = None
    residualFile = None
    downsaceldFile = None

    # print(time.time() - start_time, "seconds")


def runEvapi(year, month, day, comp, sharp, s2Mask, lstMask, tile, tempDir, path_to_temp, path_to_sharp, mvwin, cv, regrat, evap_outFolder):

    storPath_c = f'{evap_outFolder}{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET_Canopy_calc.tif'
    storPath_s = f'{evap_outFolder}{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET_Soil_calc.tif'

    storPath_c_f = f'{evap_outFolder}{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET_Canopy_func.tif'
    storPath_s_f = f'{evap_outFolder}{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET_Soil_func.tif'
    
    # path to era5 raw data
    era5_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/grib/'
    ssrd_mean_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/ssrd_mean_calc/'

    # the DEM, SLOPE, ASPECT, LAT, LON will be used to sharpen some of the era5 variables (the the resolution of the DEM)
    dem_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/DEM_GER_FORCE_WARP.tif' # epsg 4326
    slope_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/SLOPE_GER_FORCE_WARP.tif' # epsg 4326
    aspect_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/ASPECT_GER_FORCE_WARP.tif' # epsg 4326
    lat_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/LAT_GER_FORCE_WARP.tif' # epsg 4326
    lon_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/LON_GER_FORCE_WARP.tif' # epsg 4326

    # the geopotential is needed for the sharpening as well
    geopot_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/low_res/geopotential/geopotential_low_res.tif' # epsg 4326
    

    # path_base to sharpenend folder and S2_comp
    sharp_pathbase = f'{path_to_sharp}Values/'
    s2_pathbase = path_to_temp

    # the LST acquisition time should determine which sharpened LST files are associatedto be processed (as they are associated with it)
    LST_acq_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/int_format/{year}/Daily_AcqTime_{comp}_{year}_{month}.tif' # epsg 4326

    # the VZA at the time of LST acquisition is need
    VZA_at_acq_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/{comp}/{year}/Daily_VZA_{comp}_{year}_{month}.tif' # epsg 4326

    # sharpened LST
    LST_file = f'{sharp_pathbase}{comp}_{year}_{month}_{day:02d}_{mvwin}_{cv}_{regrat}_{s2Mask}_{sharp}_{lstMask}_{tile}.tif' 
    # for NDVI calculation (estimating LAI and others) and warping to S2 resolution, we use the S2 composite used for sharpening
    # S2_file = [file for file in getFilelist(s2_pathbase, 'vrt', deep=False) if f'HIGHRES_{comp}_{year}_{month}_{day:02d}' in file][0]
    S2_file = [file for file in getFilelist(s2_pathbase, 'vrt', deep=True) if 'S2' in file][0]

    # find era5 file that matches the month of LST observation
    valid_variables = sorted(list(dict.fromkeys(file.split('/')[-2] for file in getFilelist(era5_path, '.grib', deep=True) \
                                    if not any(var in file for var in ['geopotential', 'total_column_water_vapour']))))

    # get a list for those era5 files that match the year and month of the provided LST acquisition file
    era5_path_list = find_grib_file(getFilelist(era5_path, '.grib', deep=True), LST_acq_file)
    era5_path_list = [path for path in era5_path_list if any(variable in path for variable in valid_variables)] # era5 are epsg 4326 and still will be after warping to doy
    temp_pressure_checker(era5_path_list)

    # warp datasets needed for calculations to the spatial extent of the sharpened LST
    LST_acq_spatial_sub = warp_raster_to_reference(source_path=LST_acq_file, reference_path=S2_file, output_path='MEM', resampling='near')
    VZA_at_acq_file_sub = warp_raster_to_reference(source_path=VZA_at_acq_file, reference_path=S2_file, output_path='MEM', resampling='near')
    dem_sub = warp_raster_to_reference(source_path=dem_path, reference_path=S2_file, output_path='MEM', resampling='bilinear')
    slope_sub  = warp_raster_to_reference(source_path=slope_path, reference_path=S2_file, output_path='MEM', resampling='bilinear')
    aspect_sub = warp_raster_to_reference(source_path=aspect_path, reference_path=S2_file, output_path='MEM', resampling='bilinear')
    lat_sub = warp_raster_to_reference(source_path=lat_path, reference_path=S2_file, output_path='MEM', resampling='bilinear')
    lon_sub = warp_raster_to_reference(source_path=lon_path, reference_path=S2_file, output_path='MEM', resampling='bilinear')
    geopot_sub = warp_raster_to_reference(source_path=geopot_path, reference_path=S2_file, output_path='MEM', resampling='bilinear')


    # load the era5 variable into cache at LST resolution and read-in the modelled times (one time step per band)
    for path in era5_path_list:
        # print(f'processing {path}')
        # check if DEM sharpener needs to be applied
        if '100m_u_component_of_wind' in path:
            # do the warping without sharpening
            try:
                wind100_u = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, lst_acq_file=LST_acq_spatial_sub, doy=day)
            except Exception as e:
                with open(f'{tempDir}ERROR_{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return

        elif '100m_v_component_of_wind' in path:
                # do the warping without sharpening
            try:
                wind100_v = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, lst_acq_file=LST_acq_spatial_sub, doy=day)
            except Exception as e:
                with open(f'{tempDir}ERROR_{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return
        # elif 'geopotential' in path:
        #     # do the warping without sharpening
        #     geopot = get_warped_ERA5_at_doy(path_to_era_grib=path, lst_acq_file=LST_acq_file, doy=day)

        elif 'downward' in path: # terrain correction included
            try:
                ssrd, szenith, sazimuth, ssrd_nc, ssrd_mean_func = get_ssrdsc_warped_and_corrected_at_doy(path_to_ssrdsc_grib=path, reference_path=LST_acq_spatial_sub, 
                                                                                lst_acq_file=LST_acq_spatial_sub, doy=day, 
                                                                                slope_path=slope_sub,
                                                                                aspect_path=aspect_sub,
                                                                                dem_path=dem_sub,
                                                                                lat_path=lat_sub,
                                                                                lon_path=lon_sub)
            except Exception as e:
                with open(f'{tempDir}ERROR_{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return
            
        elif '2m_temperature' in path: # DEM and adiabatic sharpening, following Guzinski 2021
            try:
                air_temp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                                lst_acq_file=LST_acq_spatial_sub, doy=day,
                                                sharp_blendheight=100,
                                                sharp_DEM=dem_sub,
                                                sharp_geopot=geopot_sub,
                                                sharp_rate=STANDARD_ADIABAT,
                                                sharpener='adiabatic')
            except Exception as e:
                with open(f'{tempDir}ERROR_{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return
            
        elif '2m_dewpoint_temperature' in path: # DEM and adiabatic sharpening, following Guzinski 2021
            try:
                dew_temp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                                lst_acq_file=LST_acq_spatial_sub, doy=day,
                                                sharp_blendheight=100,
                                                sharp_DEM=dem_sub,
                                                sharp_geopot=geopot_sub,
                                                sharp_rate=MOIST_ADIABAT,
                                                sharpener='adiabatic')
            except Exception as e:
                with open(f'{tempDir}ERROR_{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return
            
        else: 
            # do warping with DEM sharpening only
            # sanity check
            if not 'surface_pressure' in path:
                raise ValueError('There is and unattended ERA5 variable in the loop - CHECK!!!!')
            else:
                try:
                    sp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                                lst_acq_file=LST_acq_spatial_sub, doy=day,
                                                sharp_DEM=dem_sub,
                                                sharp_blendheight=100,
                                                sharp_geopot=geopot_sub,
                                                sharp_temp=air_temp,
                                                sharpener='barometric') / 100
                except Exception as e:
                    with open(f'{tempDir}ERROR_{comp}_{year}_{month}_{day}_{mvwin}_{cv}_{regrat}_{lstMask}_{s2Mask}_{sharp}_{tile}_ET.log', 'a') as f:
                        f.write(f'{e}')
                    return
                

    wind_speed_20 = calc_wind_speed(wind100_u, wind100_v) # check wind_u
    

    ds = gdal.Open(f'{ssrd_mean_path}surface_solar_radiation_downward_clear_sky_{year}_{int(MONTH_TO_02D[month])}')
    ssrd_mean = ds.GetRasterBand(day).ReadAsArray() / 3600
    
    ssrd_mean_calc_20 = warp_np_to_reference(ssrd_mean, f'{ssrd_mean_path}surface_solar_radiation_downward_clear_sky_{year}_{int(MONTH_TO_02D[month])}', LST_file) # check this too!!!!!
    ssrd_mean_func_20 = ssrd_mean_func
    ssrd_20 = ssrd
    air_temp_20 = air_temp
    dew_temp_20 = dew_temp
    sp_20 = sp
    szenith_20 = szenith
    sazimuth_20 = sazimuth

    # calculate windspeed

    # load vza
    vza_ds = VZA_at_acq_file_sub
    vza_20 = vza_ds.GetRasterBand(day).ReadAsArray()

    # load sharpened LST
    lst_ds = gdal.Open(LST_file)
    lst_20 =lst_ds.GetRasterBand(1).ReadAsArray()

    del wind100_u, wind100_v, ssrd, air_temp, dew_temp, sp, szenith, sazimuth, ssrd_mean, ssrd_nc, ssrd_mean_func


    condition = (air_temp_20 > 0) & (dew_temp_20 > 0)  & (sp_20 > 0) & (szenith_20 > 0) & (sazimuth_20 > 0) & (wind_speed_20 > 0) & (lst_20 > 0) & (vza_20 > 0)
    ssrd_20[~condition] = np.nan
    ssrd_mean_calc_20[~condition] = np.nan
    ssrd_mean_func_20[~condition] = np.nan
    air_temp_20[~condition] = np.nan
    dew_temp_20[~condition] = np.nan
    sp_20[~condition] = np.nan
    szenith_20[~condition] = np.nan
    sazimuth_20[~condition] = np.nan
    wind_speed_20[~condition] = np.nan
    # lst_20 = np.ma.masked_where(~condition, lst_20)
    # vza_20 = np.ma.masked_where(~condition, vza_20)
    lst_20[~condition] = np.nan
    vza_20[~condition] = np.nan

    # calculate the NDVI from the S2 composite (following formula from force --> bandnames: (NIR - RED) / (NIR + RED))
    S2_ds = gdal.Open(S2_file)
    for idx, bname in enumerate(getBandNames(S2_file)):
        if bname == 'RED':
            red = S2_ds.GetRasterBand(1 + idx).ReadAsArray()
        elif bname == 'NIR':
            nir = S2_ds.GetRasterBand(1 + idx).ReadAsArray()
        else:
            continue
    ndvi_20 = (nir - red) / (nir + red)
    ndvi_20_ma = np.where(ndvi_20 < 0, np.nan, ndvi_20)
    # ndvi_20_ma = np.ma.masked_invalid(ndvi_20)
    # ndvi_20_ma = np.ma.masked_where(ndvi_20_ma < 0, ndvi_20_ma)
    LAI_np = 0.57*np.exp(2.33*ndvi_20)
    LAI_pos = np.where(LAI_np < 0, np.nan, LAI_np)

    # estimate canopy height from estimated LAI
    hc = hc_from_lai(LAI_pos, hc_max = 1.2, lai_max = np.nanmax(LAI_np), hc_min=0)

    # estimate long wave irradiance
    ea = meteo_utils.calc_vapor_pressure(T_K=dew_temp_20)
    L_dn = calc_longwave_irradiance(ea = ea, t_a_k = air_temp_20, p = sp_20, z_T = 100, h_C = hc) # ## does that make sense with the 100m!!!!!!!!!!!!!!!!!!!
    d_0_0 = resistances.calc_d_0(h_C=hc)
    z_0 = resistances.calc_z_0M(h_C=hc)


    # calculate shortwave radiation of soil and canopy
    difvis, difnir, fvis, fnir = net_radiation.calc_difuse_ratio(S_dn = ssrd_20, sza = np.nanmean(szenith_20))

    skyl = difvis * fvis + difnir * fnir
    S_dn_dir = ssrd_20 * (1.0 - skyl)
    S_dn_dif = ssrd_20 * skyl

    # Leaf spectral properties:{rho_vis_C: visible reflectance, tau_vis_C: visible transmittance, rho_nir_C: NIR reflectance, tau_nir_C: NIR transmittance}
    rho_vis_C=np.full(LAI_pos.shape, 0.05, np.float32)
    tau_vis_C=np.full(LAI_pos.shape, 0.08, np.float32)
    rho_nir_C=np.full(LAI_pos.shape, 0.32, np.float32)
    tau_nir_C=np.full(LAI_pos.shape, 0.33, np.float32) 

    # Soil spectral properties:{rho_vis_S: visible reflectance, rho_nir_S: NIR reflectance}
    rho_vis_S=np.full(LAI_pos.shape, 0.07, np.float32)
    rho_nir_S=np.full(LAI_pos.shape, 0.25, np.float32)

    # F = local LAI
    F = LAI_pos / 1
    # calculate clumping index
    Omega0 = clumping_index.calc_omega0_Kustas(LAI = LAI_np, f_C = 1, x_LAD=1)
    Omega = clumping_index.calc_omega_Kustas(Omega0, np.nanmean(szenith_20))
    LAI_eff = F * Omega

    Sn_C, Sn_S = net_radiation.calc_Sn_Campbell(lai = LAI_pos, sza = np.mean(szenith_20), S_dn_dir = S_dn_dir, S_dn_dif = S_dn_dif, fvis = fvis,
                                        fnir = fnir, rho_leaf_vis = rho_vis_C, tau_leaf_vis = tau_vis_C, rho_leaf_nir = rho_nir_C, 
                                        tau_leaf_nir = tau_nir_C, rsoilv = rho_vis_S, rsoiln = rho_nir_S, x_LAD=1, LAI_eff=LAI_eff)

    # calculate other roughness stuff
    z_0M, d = resistances.calc_roughness(LAI=np.nanmean(LAI_pos), h_C=hc, w_C=1, landcover=11, f_c=None)
    fg = calc_fg_gutman(ndvi = ndvi_20_ma, ndvi_min = np.nanmin(ndvi_20), ndvi_max = np.nanmax(ndvi_20))

    emis_C = 0.98
    emis_S = 0.95
    h_C = hc 
    z_u = 100
    z_T = 100

    output = TSEB.TSEB_PT(lst_20, vza_20, air_temp_20, wind_speed_20, ea, sp_20, Sn_C, Sn_S, L_dn, LAI_pos, h_C, emis_C, emis_S, 
    z_0M, d, z_u, z_T, resistance_form=None, calcG_params=None, const_L=None, 
    kB=0.0, massman_profile=None, verbose=True)

    for stori, ssrd_ras in zip([[storPath_c, storPath_s],[storPath_c_f, storPath_s_f]], [ssrd_mean_calc_20, ssrd_mean_func_20]):
        
        le_c = output[6]/ssrd_20
        heat_latent_scaled_c = ssrd_ras * le_c
        et_daily_c = TSEB.met.flux_2_evaporation(heat_latent_scaled_c, t_k=air_temp_20, time_domain=24)

        le_s = output[8]/ssrd_20
        heat_latent_scaled_s = ssrd_ras * le_s
        et_daily_s = TSEB.met.flux_2_evaporation(heat_latent_scaled_s, t_k=air_temp_20, time_domain=24)

        npTOdisk(et_daily_c, LST_file, stori[0])
        npTOdisk(et_daily_s, LST_file, stori[1])


def Sharp_Evap(tile_to_process, storFolder, path_to_slope, path_to_aspect, path_to_agro, path_to_force, time_start, time_end, compList, predList, S2mask):

    temp_dump_fold = f"{storFolder}temp/{tile_to_process.replace('_', '')}/"
    sharp_outFolder = f'{storFolder}sharpened/{tile_to_process.replace('_', '')}/'
    evap_outFolder = f'{storFolder}evap/{tile_to_process.replace('_', '')}/'

    for foldi in [temp_dump_fold, sharp_outFolder, evap_outFolder]:
        if not os.path.exists(foldi):
            os.makedirs(foldi,exist_ok=False)

    year = time_start[:4]

    # ############## make vrts for slope, aspect and agromask
    slopes = [file for file in getFilelist(path_to_slope, '.tif') if tile_to_process in file] # if any tile name is in file
    # aspect-tiles
    aspects = [file for file in getFilelist(path_to_aspect, '.tif') if tile_to_process in file] # if any tile name is in file
    # thuenen-tiles
    thuenen = [file for file in getFilelist(path_to_agro, '.tif') if tile_to_process in file] # if any tile name is in file

    # get those tiles (and composite if more than one tile is provided)
    slope_path = f'{temp_dump_fold}SLOPE.vrt'
    gdal.BuildVRT(slope_path, slopes)

    aspect_path = f'{temp_dump_fold}ASPECT.vrt'
    gdal.BuildVRT(aspect_path, aspects)

    thuenen_path = f'{temp_dump_fold}THUENEN.vrt'
    gdal.BuildVRT(thuenen_path, thuenen)

    # ################ load force and vrt
    path_to_S2_tiles = f'{path_to_force}/{year}/'
    # get a list with all available tiles
    files = getFilelist(f'{path_to_S2_tiles}/tiles', '.tif', deep=True) 
    files = [file for file in files if any(tile in file for tile in tiles_to_process)]
    date_list = check_forceTSI_compositionDates(files)

    # make the mask ready for S2 masking
    th_ds = gdal.Open(thuenen_path)
    th_arr = th_ds.GetRasterBand(1).ReadAsArray()
    mask = np.where(th_arr == -9999, 0, 1)

    colors = ['BLU', 'GRN', 'RED', 'NIR', 'RE1', 'RE2', 'RE3',  'SW1', 'SW2']

    # this should be the masterloop within the sharpend, evaping and deleting of all files but the above created vrts takes place
    for date in date_list:
        if int(time_start) <= int(date) <= int(time_end):

            # needed lists
            lowRes_files = []
            highRes_files = []
            highRes_names = []
            
            tilesS2 = [file for file in getFilelist(path_to_S2_tiles, '.tif', deep=True) if tile_to_process in file and f'{date}.tif' in file]
            tilesS2 = [t2 for col in colors for t2 in tilesS2 if col in t2]
            S2_path = f'{temp_dump_fold}S2_{date}.vrt'
            print(S2_path)
            vrt = gdal.BuildVRT(S2_path, tilesS2, separate=True)
            vrt = None
            vrt = gdal.Open(S2_path, gdal.GA_Update)  # VRT must be writable
            for idz, bname in enumerate(colors): 
                band = vrt.GetRasterBand(1+idz)
                band.SetDescription(bname)
            vrt = None

            # determine LST and incidence files associated with respective S2 composite
            band_dict = transform_compositeDate_into_LSTbands(date, 4)


            # stat used for compositing
            for comp_stat in compList: #  
                path_to_incident = f'{path_to_inci}{comp_stat}/{year}/'
                path_to_LST = f'{path_to_lst}{comp_stat}/{year}/'

                # get all LST bands that can be sharped with the S2 composite at this date (and sun angle incidence files as well, as they are dependent on that date
                LSTs = []

                for k, v in band_dict.items():
                    month = v['month']
                    band = int(v['band'])
                    v_path = f'{path_to_LST}Daily_LST_{comp_stat}_{year}_{month}.tif'
                    ds = gdal.Open(v_path, 0)

                    # export the LST for that day
                    LST_arr = ds.GetRasterBand(band).ReadAsArray() # store as single Tiff in temp
                    daily_lst_path = f'{temp_dump_fold}Daily_LST_{comp_stat}_{year}_{month}_{band:02d}.tif'
                    makeTif_np_to_matching_tif(LST_arr, v_path, daily_lst_path)
                
                    # store the paths for selecting incidence for corresponding LST
                    incid_date = f'{year}_{month}_{band:02d}.tif'

                    # incidence-tiles
                    incids = [file for file in getFilelist(path_to_incident, '.tif', deep=True) if tile_to_process in file] 
                    incid_path = incids[0]

                    # create highRes file through exapnding the vrt of S2
                    highRes_path = f'{temp_dump_fold}HIGHRES_{comp_stat}_{incid_date.split('.')[0]}.vrt'
                    gdal.BuildVRT(highRes_path, [S2_path, slope_path, aspect_path, incid_path], separate=True)

                    for predi in predList:
                        if predi == 'allpred':
                            maskVRT_water(highRes_path)
                        else:
                            maskVRT_water_and_drop_aux(highRes_path)

                        if S2mask == 1:
                            highRes_files.append(f'{highRes_path.split('.')[0]}_watermask.tif')
                            highRes_names.append(f'S2notMasked_{predi}')
                            lowRes_files.append(daily_lst_path)

                        elif S2mask == 2:
                            maskVRT(f'{highRes_path.split('.')[0]}_watermask.tif', mask, suffix=f'_S2_agromask_{predi}')
                            os.remove(f'{highRes_path.split('.')[0]}_watermask.tif')
                            highRes_files.append(f'{highRes_path.split('.')[0]}_watermask_S2_agromask_{predi}.tif')
                            lowRes_files.append(daily_lst_path)
                            highRes_names.append(f'S2Masked_{predi}')

                        elif S2mask == 3:
                            highRes_files.append(f'{highRes_path.split('.')[0]}_watermask.tif')
                            highRes_names.append(f'S2notMasked_{predi}')
                            lowRes_files.append(daily_lst_path)
                            maskVRT(f'{highRes_path.split('.')[0]}_watermask.tif', mask, suffix=f'_S2_agromask_{predi}')
                            highRes_files.append(f'{highRes_path.split('.')[0]}_watermask_S2_agromask_{predi}.tif')
                            lowRes_files.append(daily_lst_path)
                            highRes_names.append(f'S2Masked_{predi}')
            sharpList = []

            for idx, highResFilename in enumerate(highRes_files):
                    lowResFilename = lowRes_files[idx]
                    # f1 = f'{sharp_outFolder}{'/'.join(highResFilename.split('.')[0].split('_')[2:5])}/'
                    for maskname, mask_lowRes in zip(['withoutLSTmask'], ['']): # 'withLSTmask'  lowmask_bin_path
                    # for maskname, mask_lowRes in zip(['withoutLSTmask', 'withLSTmask'], ['', lowmask_bin_path]):
                        lowResMaskFilename = mask_lowRes
                        # f2 = f'{f1}{maskname}/'
                        for movWin in [15]:
                            for cv in [0]:
                                for regrat in [0.25]:
                                    kombi = f'mvwin{movWin}_cv{cv}_regrat{int(regrat*100):02d}_{highRes_names[idx]}_{maskname}'
                                    # f3 = f'{f2}{highRes_names[idx]}/'
                                    # os.makedirs(f'{f3}Residuals/', exist_ok=True)
                                    # os.makedirs(f'{f3}Values/', exist_ok=True)

                                    os.makedirs(f'{sharp_outFolder}Residuals/', exist_ok=True)
                                    os.makedirs(f'{sharp_outFolder}Values/', exist_ok=True)
                                    
                                    # sharpened_file = f'{f3}{'_'.join(highResFilename.split('.')[0].split('_')[2:6])}_{kombi}_{tile_to_process}.tif'
                                    sharpened_file = f'{sharp_outFolder}{'_'.join(highResFilename.split('.')[0].split('_')[1:5])}_{kombi}_{tile_to_process}.tif'
                                    
                                    runSharpi(highResFilename, lowResFilename, lowResMaskFilename, cv, movWin, regrat, sharpened_file, useDecisionTree = True)

                                    sharpList.append(sharpened_file)
           
            for sharped in sharpList:
                comp = sharped.split('/')[-1].split('_')[0]
                year = sharped.split('/')[-1].split('_')[1]
                month = sharped.split('/')[-1].split('_')[2]
                day = int(sharped.split('/')[-1].split('_')[3])
                mvwin = sharped.split('/')[-1].split('_')[4]
                cv = sharped.split('/')[-1].split('_')[5]
                regrat = sharped.split('/')[-1].split('_')[6]
                sharp = sharped.split('/')[-1].split('_')[8]
                s2Mask = sharped.split('/')[-1].split('_')[7]
                lstMask = sharped.split('/')[-1].split('_')[9]
                tile = '_'.join(sharped.split('/')[-1].split('.')[0].split('_')[-2:])

                runEvapi(year=year, month=month, day=day, comp=comp, sharp=sharp, s2Mask=s2Mask, lstMask=lstMask, tile=tile,
                        tempDir=trash_path, path_to_temp=temp_dump_fold, path_to_sharp=sharp_outFolder,
                        mvwin=mvwin, cv=cv, regrat=regrat, evap_outFolder=evap_outFolder)
            
                
            # at the end of the date loop --> clean up temp folder and sharp
            temp_files_vrt = [file for file in getFilelist(temp_dump_fold, '.vrt')]
            temp_files_tif = [file for file in getFilelist(temp_dump_fold, '.tif')]
            [temp_files_vrt.remove(path) for path in [slope_path, aspect_path, thuenen_path]]
            sharp_files = [file for file in getFilelist(sharp_outFolder, '.tif', deep=True)]
            [os.remove(file) for file in temp_files_vrt]
            [os.remove(file) for file in temp_files_tif]
            [os.remove(file) for file in sharp_files]


# paths
lowmask_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/THUENEN_GER_LST_WARP.tif'
lowmask_bin_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/THUENEN_GER_LST_WARP_BINARY.tif'
temp_dump = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/SharpEvap/Brandenburg/'
trash_path = f'{temp_dump}trash/'

path_to_slope = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/'
path_to_aspect = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/'
path_to_agro = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/THUENEN_2021/'

path_to_force = '/data/Aldhani/eoagritwin/force/output/Guzinski'
path_to_inci = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE2/'
path_to_lst = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/LST_composites/'

# get the Tiles for Brandenburg
bran = pd.read_csv('/data/Aldhani/eoagritwin/misc/state_tile_csv/clipped_grid_bran_tiles.csv')
tiles_to_process = createFORCEtileLIST(list(bran['Tile_X']),
                                        list(bran['Tile_Y']), True)


Sharp_Evap(tile_to_process=tiles_to_process[0], storFolder=temp_dump, path_to_slope=path_to_slope,
           path_to_aspect=path_to_aspect, path_to_agro=path_to_agro, path_to_force=path_to_force,
           time_start='20190501', time_end='20190701', compList=['maxLST'], predList=['S2only'], S2mask=2)