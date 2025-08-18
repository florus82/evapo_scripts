#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
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
import rasterio


# In[ ]:


# this might dock onto the sharpener process??


# In[2]:


# set the year, month and day to estimate evapotranspiration for
year = 2019
month = 'July'
day = 1
comp = 'minVZA'

# set output path
tempOut = '/data/Aldhani/eoagritwin/et/Auxiliary/trash/'

# set the parameters from sharpener
mvwin = 35
cv = 25
regrat = 20
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

# path to era5 raw data
era5_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/grib/'
ssrd_mean_path = '/data/Aldhani/eoagritwin/et/Auxiliary/ERA5/ssrd_mean_calc/'

# path_base to sharpenend folder and S2_comp
sharp_pathbase = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/sharpened/2514b73c638b86ba52047c6c65eea221652398b499d492380209a557c938f1ba/{comp}/{year}/{month}/{day:02d}/{lstMask}/{s2Mask}/Values/'
s2_pathbase = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/tempDump/2514b73c638b86ba52047c6c65eea221652398b499d492380209a557c938f1ba/'

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
LST_file = f'{sharp_pathbase}{comp}_{year}_{month}_{day:02d}_mvwin{mvwin}_cv{cv}_regrat{regrat}_{s2Mask}_{lstMask}.tif' # epsg:3035

# for NDVI calculation (estimating LAI and others) and warping to S2 resolution, we use the S2 composite used for sharpening
S2_file = [file for file in getFilelist(s2_pathbase, 'vrt', deep=True) if '_Cube.' in file][0]

# find era5 file that matches the month of LST observation
valid_variables = sorted(list(dict.fromkeys(file.split('/')[-2] for file in getFilelist(era5_path, '.grib', deep=True) \
                                   if not any(var in file for var in ['geopotential', 'total_column_water_vapour']))))

# get a list for those era5 files that match the year and month of the provided LST acquisition file
era5_path_list = find_grib_file(getFilelist(era5_path, '.grib', deep=True), LST_acq_file)
era5_path_list = [path for path in era5_path_list if any(variable in path for variable in valid_variables)] # era5 are epsg 4326 and still will be after warping to doy
temp_pressure_checker(era5_path_list)


# In[3]:


# warp datasets needed for calculations to the spatial extent of the sharpened LST
LST_acq_spatial_sub = warp_raster_to_reference(source_path=LST_acq_file, reference_path=S2_file, output_path='MEM', resampling='near', keepRes=True)
VZA_at_acq_file_sub = warp_raster_to_reference(source_path=VZA_at_acq_file, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)
dem_sub = warp_raster_to_reference(source_path=dem_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)
slope_sub  = warp_raster_to_reference(source_path=slope_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)
aspect_sub = warp_raster_to_reference(source_path=aspect_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)
lat_sub = warp_raster_to_reference(source_path=lat_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)
lon_sub = warp_raster_to_reference(source_path=lon_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)
geopot_sub = warp_raster_to_reference(source_path=geopot_path, reference_path=S2_file, output_path='MEM', resampling='bilinear', keepRes=True)


# In[ ]:


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
                                                                         lon_path=lon_sub)

    elif '2m_temperature' in path: # DEM and adiabatic sharpening, following Guzinski 2021
        air_temp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                          lst_acq_file=LST_acq_spatial_sub, doy=day,
                                          sharp_blendheight=100,
                                          sharp_DEM=dem_sub,
                                          sharp_geopot=geopot_sub,
                                          sharp_rate=STANDARD_ADIABAT,
                                          sharpener='adiabatic')

    elif '2m_dewpoint_temperature' in path: # DEM and adiabatic sharpening, following Guzinski 2021
        dew_temp = get_warped_ERA5_at_doy(path_to_era_grib=path, reference_path=LST_acq_spatial_sub, 
                                          lst_acq_file=LST_acq_spatial_sub, doy=day,
                                          sharp_blendheight=100,
                                          sharp_DEM=dem_sub,
                                          sharp_geopot=geopot_sub,
                                          sharp_rate=MOIST_ADIABAT,
                                          sharpener='adiabatic')

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
                                        sharpener='barometric')

# load the mean ssrd
ds = gdal.Open(f'{ssrd_mean_path}surface_solar_radiation_downward_clear_sky_{year}_{int(MONTH_TO_02D[month])}')
ssrd_mean = ds.GetRasterBand(day).ReadAsArray() / 3600


# In[ ]:


# bring all variables down to 20m
wind100_u_20 = warp_np_to_reference(wind100_u, LST_acq_spatial_sub, LST_file)
wind100_v_20 = warp_np_to_reference(wind100_v, LST_acq_spatial_sub, LST_file)
ssrd_20 = warp_np_to_reference(ssrd, LST_acq_spatial_sub, LST_file) # check this too!!!!!!
ssrd_mean_20 = warp_np_to_reference(ssrd_mean, LST_acq_spatial_sub, LST_file) # check this too!!!!!!
ssrd_nc_20 = warp_np_to_reference(ssrd_nc, LST_acq_spatial_sub, LST_file) # check this too!!!!!!
air_temp_20 = warp_np_to_reference(air_temp, LST_acq_spatial_sub, LST_file)
dew_temp_20 = warp_np_to_reference(dew_temp, LST_acq_spatial_sub, LST_file)
sp_20 = warp_np_to_reference(sp, LST_acq_spatial_sub, LST_file)
szenith_20 = warp_np_to_reference(szenith, LST_acq_spatial_sub, LST_file)
sazimuth_20 = warp_np_to_reference(sazimuth, LST_acq_spatial_sub, LST_file)

# calculate windspeed
wind_speed_20 = calc_wind_speed(wind100_u_20, wind100_v_20) # check if this really works

# load vza
vza_ds = warp_raster_to_reference(VZA_at_acq_file_sub, LST_file, output_path='MEM', resampling='bilinear')
vza_20 = vza_ds.GetRasterBand(1).ReadAsArray()
# load sharpened LST
lst_ds = gdal.Open(LST_file)
lst_20 =lst_ds.GetRasterBand(1).ReadAsArray()

del wind100_u, wind100_v, ssrd, air_temp, dew_temp, sp, wind100_u_20, wind100_v_20, szenith, sazimuth, ssrd_mean, ssrd_nc


# In[ ]:


# npTOdisk(ssrd_mean_20, LST_file, f'{tempOut}SSRD_mean.tif', bands = 1)
# npTOdisk(ssrd_20, LST_file, f'{tempOut}SSRD.tif', bands = 1)
# npTOdisk(air_temp_20, LST_file, f'{tempOut}TEMP.tif', bands = 1)
# npTOdisk(dew_temp_20, LST_file, f'{tempOut}DEW.tif', bands = 1)
# npTOdisk(sp_20, LST_file, f'{tempOut}SP.tif', bands = 1)
# npTOdisk(szenith_20, LST_file, f'{tempOut}ZEN.tif', bands = 1)
# npTOdisk(sazimuth_20, LST_file, f'{tempOut}AZI.tif', bands = 1)
# npTOdisk(wind_speed_20, LST_file, f'{tempOut}WSPEED.tif', bands = 1)
# npTOdisk(lst_20, LST_file, f'{tempOut}LST.tif', bands = 1)
# npTOdisk(vza_20, LST_file, f'{tempOut}VZA.tif', bands = 1)


# In[ ]:


#(ssrd_20 > 0)
condition = (ssrd_nc_20 > 0) & (air_temp_20 > 0) & (dew_temp_20 > 0)  & (sp_20 > 0) & (szenith_20 > 0) & (sazimuth_20 > 0) & (wind_speed_20 > 0) & (lst_20 > 0) & (vza_20 > 0)


# In[ ]:


ssrd_20[~condition] = np.nan
ssrd_mean_20[~condition] = np.nan
ssrd_nc_20[~condition] = np.nan
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


# In[ ]:


################ only for testing !!!!!!!!!!!!!! delete once the issue is fixed !!!!!!!!!!!!!!!!!!!!!!!!!!!!
ssrd_20 = ssrd_nc_20


# In[ ]:


# npTOdisk(ssrd_mean_20, LST_file, f'{tempOut}SSRD_mean_mask.tif', bands = 1)
# npTOdisk(ssrd_20, LST_file, f'{tempOut}SSRD_mask.tif', bands = 1)
# npTOdisk(air_temp_20, LST_file, f'{tempOut}TEMP_mask.tif', bands = 1)
# npTOdisk(dew_temp_20, LST_file, f'{tempOut}DEW_mask.tif', bands = 1)
# npTOdisk(sp_20, LST_file, f'{tempOut}SP_mask.tif', bands = 1)
# npTOdisk(szenith_20, LST_file, f'{tempOut}ZEN_mask.tif', bands = 1)
# npTOdisk(sazimuth_20, LST_file, f'{tempOut}AZI_mask.tif', bands = 1)
# npTOdisk(wind_speed_20, LST_file, f'{tempOut}WSPEED_mask.tif', bands = 1)
# npTOdisk(lst_20, LST_file, f'{tempOut}LST_mask.tif', bands = 1)
# npTOdisk(vza_20, LST_file, f'{tempOut}VZA_mask.tif', bands = 1)


# In[ ]:


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


# In[ ]:


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


# In[ ]:


emis_C = 0.98
emis_S = 0.95
h_C = hc 
z_u = 100
z_T = 100


# In[ ]:


output = TSEB.TSEB_PT(lst_20, vza_20, air_temp_20, wind_speed_20, ea, sp_20, Sn_C, Sn_S, L_dn, LAI_pos, h_C, emis_C, emis_S, 
                      z_0M, z_0M_array, z_u, z_T, resistance_form=None, calcG_params=None, const_L=None, 
                      kB=0.0, massman_profile=None, verbose=True)

len(output)
ld = output[6]/ssrd_20
heat_latent_scaled = ssrd_mean_20 * ld
et_daily_p = TSEB.met.flux_2_evaporation(heat_latent_scaled, t_k=air_temp_20, time_domain=24)


# In[ ]:


np.save(f'{tempOut}et.npy', et_daily_p)
npTOdisk(et_daily_p, LST_file, f'{tempOut}et.tif')


# In[ ]:


get_ipython().system('jupyter nbconvert --to script estimate_Evapotrans.ipynb --output /home/potzschf/repos/evapo_scripts/Guzinski/py/estimate_Evapotrans')

