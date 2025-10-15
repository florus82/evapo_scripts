import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 

tempOut = '/data/Aldhani/eoagritwin/et/Auxiliary/trash/Compare_ET/'

# get all folders
folders = [os.path.join(tempOut, entry) for entry in os.listdir(tempOut) if os.path.isdir(os.path.join(tempOut, entry))]

results = []
i = 0
for fold in folders:
    files = getFilelist(fold, 'tif') 
    if len(files) != 29:
        i += 1
        continue
    else:
        i += 1
        print(i)
        for file in files[:]:
            if 'ET_Canopy' in file:
                ds = gdal.Open(file)
                et_C = ds.GetRasterBand(1).ReadAsArray()
                files.remove(file)
            elif 'ET_Soil' in file:
                ds = gdal.Open(file)
                et_S = ds.GetRasterBand(1).ReadAsArray()
                files.remove(file)
            else:
                pass

   

    y = fold.split('/')[-1].split('_')[4]
    m = MONTH_TO_02D[fold.split('/')[-1].split('_')[5]]
    d = fold.split('/')[-1].split('_')[6]

    if fold.split('_')[-2] == 'withLSTmask':
        lstmask = 'yes'
    else:
        lstmask = 'no'

    if fold.split('_')[-1] == 'S2notMasked':
        s2mask = 'no'
    else:
        s2mask = 'yes'

    for file in files:

        # load the variable
        var_name = file.split('/')[-1].split('.')[0]
        ds = gdal.Open(file)
        arr = ds.GetRasterBand(1).ReadAsArray()

        # take care of LST acq file unix seconds 
        if var_name == 'LST_acq_spatial':
            arr = pd.to_datetime(arr, unit='s').values.reshape(arr.shape)
        # loop over both et products and extract values
        for et_name, et_var in zip(['Soil', 'Canopy'], [et_S, et_C]):
            # loop over inside/outside mask for easier appending
            mask0 = (et_var <= 0) & (~np.isnan(et_var))
            mask_inv = (et_var > 0) & (~np.isnan(et_var))

            for maskname, master_mask in zip(['Where0', 'WhereNOT0'], [mask0, mask_inv]):

                if  np.sum(master_mask) > 0:
                    v_max = np.nanmax(arr[master_mask])
                    v_min = np.nanmin(arr[master_mask])
                    v_p25, v_p50, v_p75 = np.nanpercentile(arr[master_mask], [25, 50, 75])
                else:
                    v_max = v_min = v_p25 = v_p50 = v_p75 = np.nan

                # Append result
                results.append({
                    'Tile': '_'.join(fold.split('/')[-1].split('_')[:2]),
                    'preds': fold.split('/')[-1].split('_')[2],
                    'composite': fold.split('/')[-1].split('_')[3],
                    'date': f'{d}-{m}-{y}',
                    'LST_Mask': lstmask,
                    'S2_Mask': s2mask, 
                    'ET_var': et_name,
                    'Mask_Name': maskname,
                    'name_var': var_name,
                    'vmax': v_max,
                    'vmin': v_min,
                    'vp25': v_p25,
                    'vp50': v_p50,
                    'vp75': v_p75,
                    'vlenght': arr[master_mask].shape[0],
                })

# convert and export
df = pd.DataFrame(results)
df.to_csv('/home/potzschf/repos/evapo_scripts/Analysis/extract_evap0.csv', index=False)