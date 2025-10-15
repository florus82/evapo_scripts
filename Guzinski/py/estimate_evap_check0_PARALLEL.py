import sys

from matplotlib.pylab import dtype
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
from helperToolz.mirmazloumi import *
from pympler import asizeof
from other_repos.pyTSEB.pyTSEB import meteo_utils
from other_repos.pyTSEB.pyTSEB import resistances
from other_repos.pyTSEB.pyTSEB import net_radiation
from other_repos.pyTSEB.pyTSEB import clumping_index 
from other_repos.pyTSEB.pyTSEB import TSEB
import time
from joblib import Parallel, delayed

def run_evapi(year, month, day, comp, sharp, s2Mask, lstMask, tileID, uniqueFolder):

       # set output path
    tempOut = '/data/Aldhani/eoagritwin/et/Auxiliary/trash/'
    
    # set the parameters from sharpener
    mvwin = 15
    cv = 0
    regrat = 25


    outPath = f'{tempOut}{tileID}_{sharp}_{comp}_{year}_{month}_{day:02d}_{lstMask}_{s2Mask}/'
    storPath_c = f'{outPath}{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET_Canopy.tif'
    storPath_s = f'{outPath}{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET_Soil.tif'
    if os.path.isfile(storPath_c):
        return print('file already there')
    else:
        os.makedirs(outPath, exist_ok=True)



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
    if tileID == 'X67_Y42':
        sharp_pathbase = f"/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/sharpened3/{sharp}/{uniqueFolder.split('/')[-1]}/{comp}/{year}/{month}/{day:02d}/{lstMask}/{s2Mask}/Values/"
    else:
        sharp_pathbase = f"/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/sharpened2/{sharp}/{uniqueFolder.split('/')[-1]}/{comp}/{year}/{month}/{day:02d}/{lstMask}/{s2Mask}/Values/"
    s2_pathbase = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/{uniqueFolder}/'


    # the LST acquisition time should determine which sharpened LST files are associatedto be processed (as they are associated with it)
    LST_acq_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/int_format/{year}/Daily_AcqTime_{comp}_{year}_{month}.tif' # epsg 4326

    # the VZA at the time of LST acquisition is need
    VZA_at_acq_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/{comp}/{year}/Daily_VZA_{comp}_{year}_{month}.tif' # epsg 4326

    # sharpened LST
    LST_file = f'{sharp_pathbase}{comp}_{year}_{month}_{day:02d}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{s2Mask}_{lstMask}.tif' # epsg:3035
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

    LST_acq_spatial_arr = LST_acq_spatial_sub.GetRasterBand(day).ReadAsArray()
    VZA_at_acq_file_arr = VZA_at_acq_file_sub.GetRasterBand(day).ReadAsArray()
    dem_arr= dem_sub.GetRasterBand(1).ReadAsArray()
    slope_arr = slope_sub.GetRasterBand(1).ReadAsArray()
    aspect_arr = aspect_sub.GetRasterBand(1).ReadAsArray()
    lat_arr = lat_sub.GetRasterBand(1).ReadAsArray()
    lon_arr = lon_sub.GetRasterBand(1).ReadAsArray()
    geopot_arr = geopot_sub.GetRasterBand(1).ReadAsArray()


    sufix = ''

    npTOdisk(LST_acq_spatial_arr, S2_file, f'{outPath}LST_acq_spatial_sub{sufix}.tif', bands = 1, d_type=gdal.GDT_Int64)
    npTOdisk(VZA_at_acq_file_arr, S2_file, f'{outPath}VZA_at_acq_file_sub{sufix}.tif', bands = 1)
    npTOdisk(dem_arr, S2_file, f'{outPath}dem_sub{sufix}.tif', bands = 1)
    npTOdisk(slope_arr, S2_file, f'{outPath}slope_sub{sufix}.tif', bands = 1)
    npTOdisk(aspect_arr, S2_file, f'{outPath}aspect_sub{sufix}.tif', bands = 1)
    npTOdisk(lat_arr, S2_file, f'{outPath}lat_sub{sufix}.tif', bands = 1)
    npTOdisk(lon_arr, S2_file, f'{outPath}lon_sub{sufix}.tif', bands = 1)
    npTOdisk(geopot_arr, S2_file, f'{outPath}geopot_sub{sufix}.tif', bands = 1)


    # load the era5 variable into cache at LST resolution and read-in the modelled times (one time step per band)
    for path in era5_path_list:
        print(f'processing {path}')
        # check if DEM sharpener needs to be applied
        if '100m_u_component_of_wind' in path:
            # do the warping without sharpening
            try:
                wind100_u = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, lst_acq_file=LST_acq_spatial_sub, doy=day)
            except Exception as e:
                with open(f'{tempOut}ERROR_{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return

        elif '100m_v_component_of_wind' in path:
                # do the warping without sharpening
            try:
                wind100_v = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, lst_acq_file=LST_acq_spatial_sub, doy=day)
            except Exception as e:
                with open(f'{tempOut}ERROR_{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET.log', 'a') as f:
                    f.write(f'{e}')
                return
        # elif 'geopotential' in path:
        #     # do the warping without sharpening
        #     geopot = get_warped_ERA5_at_doy(path_to_era_grib=path, lst_acq_file=LST_acq_file, doy=day)

        elif 'downward' in path: # terrain correction included
            try:
                ssrd, szenith, sazimuth, ssrd_nc, ssrd_mean = get_ssrdsc_warped_and_corrected_at_doy(path_to_ssrdsc_grib=path, reference_path=LST_acq_spatial_sub, 
                                                                                lst_acq_file=LST_acq_spatial_sub, doy=day, 
                                                                                slope_path=slope_sub,
                                                                                aspect_path=aspect_sub,
                                                                                dem_path=dem_sub,
                                                                                lat_path=lat_sub,
                                                                                lon_path=lon_sub)
            except Exception as e:
                with open(f'{tempOut}ERROR_{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET.log', 'a') as f:
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
                with open(f'{tempOut}ERROR_{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET.log', 'a') as f:
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
                with open(f'{tempOut}ERROR_{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET.log', 'a') as f:
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
                    with open(f'{tempOut}ERROR_{comp}_{year}_{month}_{day}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{lstMask}_{s2Mask}_{sharp}_{tileID}_ET.log', 'a') as f:
                        f.write(f'{e}')
                    return
                
    wind_speed_20 = calc_wind_speed(wind100_u, wind100_v) # check wind_u
    ssrd_20 = ssrd
    # ssrd_mean_20 = warp_np_to_reference(ssrd_mean, f'{ssrd_mean_path}surface_solar_radiation_downward_clear_sky_{year}_{int(MONTH_TO_02D[month])}', LST_file) # check this too!!!!!!
    ssrd_mean_20 = ssrd_mean.reshape((1500,1500))

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

    del wind100_u, wind100_v, ssrd, air_temp, dew_temp, sp, szenith, sazimuth, ssrd_mean, ssrd_nc

    npTOdisk(ssrd_mean_20, LST_file, f'{outPath}SSRD_mean2.tif', bands = 1)
    npTOdisk(ssrd_20, LST_file, f'{outPath}SSRD.tif', bands = 1)
    npTOdisk(air_temp_20, LST_file, f'{outPath}TEMP.tif', bands = 1)
    npTOdisk(dew_temp_20, LST_file, f'{outPath}DEW.tif', bands = 1)
    npTOdisk(sp_20, LST_file, f'{outPath}SP.tif', bands = 1)
    npTOdisk(szenith_20, LST_file, f'{outPath}ZEN.tif', bands = 1)
    npTOdisk(sazimuth_20, LST_file, f'{outPath}AZI.tif', bands = 1)
    npTOdisk(wind_speed_20, LST_file, f'{outPath}WSPEED.tif', bands = 1)
    npTOdisk(lst_20, LST_file, f'{outPath}LST.tif', bands = 1)
    npTOdisk(vza_20, LST_file, f'{outPath}VZA.tif', bands = 1)

    condition = (air_temp_20 > 0) & (dew_temp_20 > 0)  & (sp_20 > 0) & (szenith_20 > 0) & (sazimuth_20 > 0) & (wind_speed_20 > 0) & (lst_20 > 0) & (vza_20 > 0)
    ssrd_20[~condition] = np.nan
    ssrd_mean_20[~condition] = np.nan
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

    npTOdisk(ndvi_20, LST_file, f'{outPath}ndvi.tif')
    npTOdisk(ndvi_20_ma, LST_file, f'{outPath}ndvi_pos.tif')
    npTOdisk(LAI_np, LST_file, f'{outPath}LAI.tif')
    npTOdisk(LAI_pos, LST_file, f'{outPath}LAI_pos.tif')
    npTOdisk(hc, LST_file, f'{outPath}canopy_height.tif')
    npTOdisk(ea, LST_file, f'{outPath}vapor_pressure.tif')
    npTOdisk(L_dn, LST_file, f'{outPath}longwave_radiation.tif')
    npTOdisk(d_0_0, LST_file, f'{outPath}resistanceD.tif')
    npTOdisk(z_0, LST_file, f'{outPath}resistanceZ.tif')


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

    le_c = output[6]/ssrd_20
    heat_latent_scaled_c = ssrd_mean_20 * le_c
    et_daily_c = TSEB.met.flux_2_evaporation(heat_latent_scaled_c, t_k=air_temp_20, time_domain=24)

    le_s = output[8]/ssrd_20
    heat_latent_scaled_s = ssrd_mean_20 * le_s
    et_daily_s = TSEB.met.flux_2_evaporation(heat_latent_scaled_s, t_k=air_temp_20, time_domain=24)

    
    npTOdisk(et_daily_c, LST_file, storPath_c)
    npTOdisk(et_daily_s, LST_file, storPath_s)

ncores = 10
joblist = []

# set the year, month and day to estimate evapotranspiration for
year = 2019
month = 'July'
for day in range(19, 28, 1):
    for comp in ['maxLST', 'minVZA']:
        for sharp in ['allpred', 'S2only']:
            for s2Mask in ['S2Masked', 'S2notMasked']:
                for lstMask in ['withoutLSTmask']:
                    for tileID, uniqueFolder in zip(['X64_Y42', 'X67_Y42'], ['tempDump2/ff3f2c872c08977466e5a8dc306d2d2aabc77ad995b0716c30d1c57d0004ebfd', 'tempDump3/bccec09fdbc7a94b5e2bd6c125d9d1174c2509360a2ce529f48e3512655454b5']):
                        joblist.append([year, month, day, comp, sharp, s2Mask, lstMask, tileID, uniqueFolder])

print(f'\n{len(joblist)} times will be vaped\n')


if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time:" + starttime)
    print("")

    Parallel(n_jobs=ncores)(delayed(run_evapi)(job[0], job[1], job[2], job[3], job[4], job[5], job[6], job[7], job[8]) for job in joblist)


print("")
endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print("--------------------------------------------------------")
print("--------------------------------------------------------")
print("start : " + starttime)
print("end: " + endtime)
print("")