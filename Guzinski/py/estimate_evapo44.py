import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
from helperToolz.mirmazloumi import *
from other_repos.pyTSEB.pyTSEB import meteo_utils
from other_repos.pyTSEB.pyTSEB import resistances
from other_repos.pyTSEB.pyTSEB import net_radiation
from other_repos.pyTSEB.pyTSEB import clumping_index 
from other_repos.pyTSEB.pyTSEB import TSEB
import time
from joblib import Parallel, delayed

ncores = 22

def run_evapo(LST_file, day, year, comp, month, tempOut, cdx):    
    # path to era5 raw data
    era5_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/grib/'
    ssrd_mean_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/ssrd_mean_calc/'

    # the LST acquisition time should determine which sharpened LST files are associatedto be processed (as they are associated with it)
    LST_acq_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/int_format/{year}/Daily_AcqTime_{comp}_{year}_{month}.tif' # epsg 4326

    # the VZA at the time of LST acquisition is need
    VZA_at_acq_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/{comp}/{year}/Daily_VZA_{comp}_{year}_{month}.tif' # epsg 4326

    # the DEM, SLOPE, ASPECT, LAT, LON will be used to sharpen some of the era5 variables (the the resolution of the DEM)
    dem_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/DEM_GER_LST_WARP.tif' # epsg 4326
    slope_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/SLOPE_GER_LST_WARP.tif' # epsg 4326
    aspect_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/ASPECT_GER_LST_WARP.tif' # epsg 4326
    lat_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/LAT_GER_LST_WARP.tif' # epsg 4326
    lon_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/LON_GER_LST_WARP.tif' # epsg 4326

    # the geopotential is needed for the sharpening as well
    geopot_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/tiff/low_res/geopotential/geopotential_low_res.tif' # epsg 4326

    # sharpened LST
    LST_file = LST_file
    # for NDVI calculation (estimating LAI and others) and warping to S2 resolution, we use the S2 composite used for sharpening

    # S2_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/tempDump/{LST_file.split(f'/{comp}')[0].split('/')[-1]}/HIGHRES_{comp}_{year}_{month}_{day:02d}_watermask.tif'
    S2_file = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/tempDump/{LST_file.split(f'/{comp}')[0].split('/')[-1]}/S2_20190512.vrt'
    # find era5 file that matches the month of LST observation
    valid_variables = sorted(list(dict.fromkeys(file.split('/')[-2] for file in getFilelist(era5_path, '.grib', deep=True) \
                                    if not any(var in file for var in ['geopotential', 'total_column_water_vapour']))))

    # get a list for those era5 files that match the year and month of the provided LST acquisition file
    era5_path_list = find_grib_file(getFilelist(era5_path, '.grib', deep=True), LST_acq_file)
    era5_path_list = [path for path in era5_path_list if any(variable in path for variable in valid_variables)] # era5 are epsg 4326 and still will be after warping to doy
    temp_pressure_checker(era5_path_list)


    # warp datasets needed for calculations to the spatial extent of the sharpened LST
    LST_acq_spatial_sub = warp_raster_to_reference(source_path=LST_acq_file, reference_path=S2_file, output_path='MEM', resampling='near', keepRes=False)
    VZA_at_acq_file_sub = warp_raster_to_reference(source_path=VZA_at_acq_file, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)
    dem_sub = warp_raster_to_reference(source_path=dem_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)
    slope_sub  = warp_raster_to_reference(source_path=slope_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)
    aspect_sub = warp_raster_to_reference(source_path=aspect_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)
    lat_sub = warp_raster_to_reference(source_path=lat_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)
    lon_sub = warp_raster_to_reference(source_path=lon_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)
    geopot_sub = warp_raster_to_reference(source_path=geopot_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=False)

    # load the era5 variable into cache at LST resolution and read-in the modelled times (one time step per band)
    for path in era5_path_list:
        print(f'processing {path}')
        # check if DEM sharpener needs to be applied
        if '100m_u_component_of_wind' in path:
            # do the warping without sharpening
            wind100_u = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, lst_acq_file=LST_acq_spatial_sub, doy=day)

        elif '100m_v_component_of_wind' in path:
                # do the warping without sharpening
            wind100_v = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, lst_acq_file=LST_acq_spatial_sub, doy=day)

        # elif 'geopotential' in path:
        #     # do the warping without sharpening
        #     geopot = get_warped_ERA5_at_doy(path_to_era_grib=path, lst_acq_file=LST_acq_file, doy=day)

        elif 'downward' in path: # terrain correction included
            ssrd, szenith, sazimuth, ssrd_nc = get_ssrdsc_warped_and_corrected_at_doy(path_to_ssrdsc_grib=path, reference_path=LST_acq_spatial_sub, 
                                                                            lst_acq_file=LST_acq_spatial_sub, doy=day, 
                                                                            slope_path=slope_sub,
                                                                            aspect_path=aspect_sub,
                                                                            dem_path=dem_sub,
                                                                            lat_path=lat_sub,
                                                                            lon_path=lon_sub, nodat=0)

        elif '2m_temperature' in path: # DEM and adiabatic sharpening, following Guzinski 2021
            air_temp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                            lst_acq_file=LST_acq_spatial_sub, doy=day,
                                            sharp_blendheight=100,
                                            sharp_DEM=dem_sub,
                                            sharp_geopot=geopot_sub,
                                            sharp_rate=STANDARD_ADIABAT,
                                            sharpener='adiabatic', nodat=0)

        elif '2m_dewpoint_temperature' in path: # DEM and adiabatic sharpening, following Guzinski 2021
            dew_temp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                            lst_acq_file=LST_acq_spatial_sub, doy=day,
                                            sharp_blendheight=100,
                                            sharp_DEM=dem_sub,
                                            sharp_geopot=geopot_sub,
                                            sharp_rate=MOIST_ADIABAT,
                                            sharpener='adiabatic', nodat=0)

        else: 
            # do warping with DEM sharpening only
            # sanity check
            if not 'surface_pressure' in path:
                raise ValueError('There is and unattended ERA5 variable in the loop - CHECK!!!!')
            else:
                sp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                            lst_acq_file=LST_acq_spatial_sub, doy=day,
                                            sharp_DEM=dem_sub,
                                            sharp_blendheight=100,
                                            sharp_geopot=geopot_sub,
                                            sharp_temp=air_temp,
                                            sharpener='barometric', nodat=0)

    # load the mean ssrd
    ds = gdal.Open(f'{ssrd_mean_path}surface_solar_radiation_downward_clear_sky_{year}_{int(MONTH_TO_02D[month])}')
    ssrd_mean = ds.GetRasterBand(day).ReadAsArray() / 3600


    # bring all variables down to 20m
    wind100_u_20 = wind100_u
    wind100_v_20 = wind100_v
    ssrd_20 = ssrd
    ssrd_mean_20 = warp_np_to_reference(ssrd_mean, f'{ssrd_mean_path}surface_solar_radiation_downward_clear_sky_{year}_{int(MONTH_TO_02D[month])}', LST_file) # check this too!!!!!!
    air_temp_20 = air_temp
    dew_temp_20 = dew_temp
    sp_20 = sp
    szenith_20 = szenith
    sazimuth_20 = sazimuth

    # calculate windspeed
    wind_speed_20 = calc_wind_speed(wind100_u_20, wind100_v_20) # check if this really works

    # load vza
    vza_20 = VZA_at_acq_file_sub.GetRasterBand(day).ReadAsArray()
    # load sharpened LST
    lst_ds = gdal.Open(LST_file)
    lst_20 =lst_ds.GetRasterBand(1).ReadAsArray()


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
    ndvi_20_ma = np.ma.masked_invalid(ndvi_20)
    ndvi_20_ma = np.ma.masked_where(ndvi_20_ma < 0, ndvi_20_ma)
    LAI_np = 0.57*np.exp(2.33*ndvi_20)
    LAI_pos = np.ma.masked_where(LAI_np < 0, LAI_np)

    # estimate canopy height from estimated LAI
    hc = hc_from_lai(LAI_pos, hc_max = 1.2, lai_max = np.nanmax(LAI_np), hc_min=0)

    # estimate long wave irradiance
    ea = meteo_utils.calc_vapor_pressure(T_K=dew_temp_20)
    L_dn = calc_longwave_irradiance(ea = ea, t_a_k = air_temp_20, p = sp_20, z_T = 100, h_C = hc) # ## does that make sense with the 100m!!!!!!!!!!!!!!!!!!!
    d_0_0 = resistances.calc_d_0(h_C=hc)
    z_0 = resistances.calc_z_0M(h_C=hc)

    # calculate shortwave radiation of soil and canopy
    difvis, difnir, fvis, fnir = net_radiation.calc_difuse_ratio(S_dn = ssrd_20, sza = np.mean(szenith_20))

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
    Omega = clumping_index.calc_omega_Kustas(Omega0, np.mean(szenith_20))
    LAI_eff = F * Omega

    Sn_C, Sn_S = net_radiation.calc_Sn_Campbell(lai = LAI_pos, sza = np.mean(szenith_20), S_dn_dir = S_dn_dir, S_dn_dif = S_dn_dif, fvis = fvis,
                                        fnir = fnir, rho_leaf_vis = rho_vis_C, tau_leaf_vis = tau_vis_C, rho_leaf_nir = rho_nir_C, 
                                        tau_leaf_nir = tau_nir_C, rsoilv = rho_vis_S, rsoiln = rho_nir_S, x_LAD=1, LAI_eff=LAI_eff)

    # calculate other roughness stuff
    z_0M, d = resistances.calc_roughness(LAI=np.nanmean(LAI_pos), h_C=hc, w_C=1, landcover=11, f_c=None)
    z_0M_array = np.full(LAI_pos.shape, z_0M)
    d_array = np.full(LAI_pos.shape, d)
    fg = calc_fg_gutman(ndvi = ndvi_20_ma, ndvi_min = np.nanmin(ndvi_20), ndvi_max = np.nanmax(ndvi_20))

    emis_C = 0.98
    emis_S = 0.95
    h_C = hc 
    z_u = 100
    z_T = 100

    output = TSEB.TSEB_PT(lst_20, vza_20, air_temp_20, wind_speed_20, ea, sp_20, Sn_C, Sn_S, L_dn, LAI_pos, h_C, emis_C, emis_S, 
                      z_0M, z_0M_array, z_u, z_T, resistance_form=None, calcG_params=None, const_L=None, 
                      kB=0.0, massman_profile=None, verbose=True)

    ld = output[6]/ssrd_20
    heat_latent_scaled = ssrd_mean_20 * ld
    et_daily_p = TSEB.met.flux_2_evaporation(heat_latent_scaled, t_k=air_temp_20, time_domain=24)

    storPath = f"{tempOut}{LST_file.split('/')[9]}_{LST_file.split('/')[-1].split('.')[0]}_ET_{cdx}.tif"
    npTOdisk(et_daily_p, LST_file, storPath)


base_path = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/sharpened2/'

folders = ['S2onlyshort', 'allpredshort']

comps = ['maxLST', 'minVZA']

days = [f'May_{i:02d}' for i in range(8,17,1)]

masks = ['S2Masked_withLSTmask', 'S2Masked_withoutLSTmask', 'S2notMasked_withLSTmask', 'S2notMasked_withoutLSTmask']

file_dict = {}

for folder in folders:
    for comp in comps:
        for day in days:
            for mask in masks:
                key = (folder, comp, day, mask)
                file_list = [
                    file for file in getFilelist(base_path + folder, '.tif', deep=True)
                    if comp in file and day in file and mask in file and 'resid' not in file
                ]
                file_dict[key] = file_list


# set the year, month and day to estimate evapotranspiration for
year = 2019
month = 'May'
day = 14
comp = 'minVZA'
folder = folders[0]

# set the parameters from sharpener
mvwin = 0
cv = 15
regrat = 25

s2_masked = False
lst_masked = False

if s2_masked:
    s2Mask = 'S2Masked'
else:
    s2Mask = 'S2notMasked'

if lst_masked:
    lstMask = 'withLSTmask'
else:
    lstMask = 'withoutLSTmask'


LST_files = file_dict[(folder, comp, f'{month}_{day:02d}', f'{s2Mask}_{lstMask}')]

# set output path
outDir = f'/data/Aldhani/eoagritwin/et/Auxiliary/S2_ETa/44/single_tiles/{folder}_{comp}_{month}_{day:02d}_{s2Mask}_{lstMask}/'
os.makedirs(outDir, exist_ok=True)

joblist = []

for cdxx, LST_file in enumerate(LST_files):
    joblist.append([LST_file, day, year, comp, month, outDir, cdxx])

print(f'\n{len(joblist)} tiles to evap\n')



if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time:" + starttime)
    print("")

    Parallel(n_jobs=ncores)(delayed(run_evapo)(job[0], job[1], job[2], job[3], job[4], job[5], job[6]) for job in joblist)


    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start : " + starttime)
    print("end: " + endtime)
    print("")

