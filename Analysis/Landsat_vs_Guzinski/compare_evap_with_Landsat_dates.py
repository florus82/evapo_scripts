import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
import geopandas as gpd
from collections import defaultdict
from joblib import Parallel, delayed
from datetime import datetime, timedelta
import time


np.random.seed(42)
Landsat_pattern = r"\d{4}_\d{2}_\d{1,2}"
Landsat_folder = '/data/Aldhani/eoagritwin/et/Landsat/daily_extracts/Brandenburg/2019/'

guz_pattern = r"20\d{2}_[A-Z][a-z]+_\d{2}"

daily_LS = [file for file in getFilelist(Landsat_folder, '.tif', deep=True)]
all_LS_dates = re.findall(Landsat_pattern, " ".join(daily_LS))
all_LS_dates = [datetime.strptime(str(compDate), '%Y_%m_%d') for compDate in all_LS_dates]

mask_ds = gdal.Open('/data/Aldhani/eoagritwin/fields/Auxiliary/grid_search/Brandenburg/quick_n_dirty/Fields_as_mask_pixel_20.tif')
mask_arr = mask_ds.GetRasterBand(1).ReadAsArray()
mask_arr[mask_arr>0] = 1

# Number of strata and samples per stratum
n_strata = 10
samples_per_stratum = 1000

results = []


all_guz_Dates = []

for comp in ['maxLST', 'minVZA']:
    pathi = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/SharpEvap/Brandenburg/FirstShot/evap/'
    vrt_folder= f'{pathi}{comp}/vrt/'

    daily_Canopy = [file for file in getFilelist(vrt_folder, '.vrt', deep=True) if '_calc_' in file and 'Canopy' in file]
    daily_Soil = [file for file in getFilelist(vrt_folder, '.vrt', deep=True) if '_calc_' in file and 'Soil' in file]

    for daily_C in daily_Canopy:
        for daily_S in daily_Soil:

            # find date matching soil and canopy files
            if daily_C.split('Canopy_')[-1] == daily_S.split('Soil_')[-1]:

                # get the date as time object
                g_date = datetime.strptime(str(daily_C.split('Canopy_calc_')[-1].split('.')[0]), '%Y_%B_%d')

                # check if a landsat image exists for that date
                if g_date in all_LS_dates:
                    print(g_date)

                    # add soil and canopy
                    arr_C = stackReader(checkPath(daily_C))
                    arr_S = stackReader(checkPath(daily_S))
                    guz_arr = arr_C + arr_S
                    guz_arr_masked = guz_arr * mask_arr
                    guz_arr_masked = np.where(guz_arr_masked > 0, guz_arr_masked, np.nan)

                    # load and warp landsat to match estimation
                    LS_Match = daily_LS[all_LS_dates.index(g_date)]
                    LandS_arr = stackReader(checkPath(LS_Match))
                    warped_LandS_arr = warp_np_to_reference(LandS_arr, LS_Match, target_tif_path=daily_C)
                    warped_masked_LandS = warped_LandS_arr * mask_arr

                    # Compute percentile thresholds
                    percentiles = np.nanpercentile(guz_arr_masked, np.linspace(0, 100, n_strata + 1))
                    print(percentiles)
                    for i in range(n_strata):

                        lower, upper = percentiles[i], percentiles[i + 1]
                        stratum_mask = (guz_arr_masked >= lower) & (guz_arr_masked < upper) & (warped_masked_LandS > 0) & (~np.isnan(warped_masked_LandS))

                        idx = np.argwhere(stratum_mask)

                        chosen_idx = idx[np.random.choice(len(idx), min(samples_per_stratum, len(idx)), replace=False)]

                        chosen_landsat = warped_masked_LandS[chosen_idx[:,0], chosen_idx[:,1]]
                        chosen_evapo_est = guz_arr_masked[chosen_idx[:,0], chosen_idx[:,1]]

                        for pix in range(len(chosen_landsat)):

                            results.append({
                                'Date': g_date.strftime('%Y_%m_%d'),
                                'Stratum': f'{i}_{i+1}th',
                                'Landsat': chosen_landsat[pix],
                                'Guzinski': chosen_evapo_est[pix],
                                'row': chosen_idx[pix][0],
                                'col': chosen_idx[pix][1],
                                'comp': comp
                            })

                else:
                    print('no landsat scene at that day :(')

# convert and export
df = pd.DataFrame(results)
df.to_csv(f'daily_extracts.csv', index=False)