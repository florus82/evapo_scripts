#!/usr/bin/env python
# coding: utf-8

# In[2]:


import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
from other_repos.pyDMS.pyDMS.pyDMS import *
import time
from joblib import Parallel, delayed

os.environ["GDAL_MAX_DATASET_POOL_SIZE"] = "600"

ncores = 40

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

    print("Training regressor...")
    disaggregator.trainSharpener()
    print("Sharpening...")
    downscaledFile = disaggregator.applySharpener(highResFilename, lowResFilename)
    print("Residual analysis...")
    residualImage, correctedImage = disaggregator.residualAnalysis(downscaledFile, lowResFilename,
                                                                    lowResMaskFilename,
                                                                    doCorrection=True)
    print("Saving output...")

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

    print(time.time() - start_time, "seconds")


# In[3]:


# paths
lowmask_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/THUENEN_GER_LST_WARP.tif'
lowmask_bin_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/THUENEN_GER_LST_WARP_BINARY.tif'
temp_dump = '/data/Aldhani/eoagritwin/et/Sentinel3//LST/LST_values/tempDump/'

# create vrts of slope, aspect and landcover (for masking)
path_to_slope = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/'
path_to_aspect = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/'
path_to_agro = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/THUENEN_2021/'

# get all highRes datasets: S2 composites, aspect, ratio and incedence for the tile/time and stack them (in a composite if more than one tile is processed)
# set the highRes S2 tiles that will be processed
# tiles_to_process = createFORCEtileLIST(list(bran['Tile_X']),
#                                         list(bran['Tile_Y']))

Tiles_X = [[i] for i in range(60,71,1)] * 4
Tiles_Y = [[i] for i in range(40,44,1) for _ in range(int(len(Tiles_X)/4))]

for Tile_X, Tile_Y in zip(Tiles_X, Tiles_Y):
    print(Tile_X)
    print(Tile_Y)

    tiles_to_process = createFORCEtileLIST(Tile_X, Tile_Y)


    # In[4]:


    # make a specific folder for this run and store the info together
    csv_path = f'{temp_dump}folders.csv'
    rand_foldname = getUniqueIDfromTILESXY(Tile_X, Tile_Y)
    temp_dump_fold = f'{temp_dump}{rand_foldname}/'

    if not os.path.exists(temp_dump_fold):
        os.makedirs(temp_dump_fold,exist_ok=False)

        df = pd.DataFrame({
                'Folder': temp_dump_fold,
                'Tile_X': Tile_X,
                'Tile_Y': Tile_Y   
            })

        if not os.path.exists(csv_path):
            df.to_csv(csv_path, index=False)
        else:
            df_exist = pd.read_csv(csv_path)
            df_new = pd.concat([df_exist, df], ignore_index=True)
            df_new.to_csv(csv_path, index=False)
            # add the lines

        # slope-tiles
        slopes = [file for file in getFilelist(path_to_slope, '.tif') if any(tile in file for tile in tiles_to_process)] # if any tile name is in file
        # aspect-tiles
        aspects = [file for file in getFilelist(path_to_aspect, '.tif') if any(tile in file for tile in tiles_to_process)] # if any tile name is in file
        # thuenen-tiles
        thuenen = [file for file in getFilelist(path_to_agro, '.tif') if any(tile in file for tile in tiles_to_process)] # if any tile name is in file

        # get those tiles (and composite if more than one tile is provided)
        slope_path = f'{temp_dump_fold}SLOPE_vrt'
        gdal.BuildVRT(slope_path, slopes)

        aspect_path = f'{temp_dump_fold}ASPECT.vrt'
        gdal.BuildVRT(aspect_path, aspects)

        thuenen_path = f'{temp_dump_fold}THUENEN.vrt'
        gdal.BuildVRT(thuenen_path, thuenen)
    else:
        print('Tile combination already processed before - slope and aspect vrt should already be present')
        slope_path = f'{temp_dump_fold}SLOPE.vrt'
        aspect_path = f'{temp_dump_fold}ASPECT.vrt'
        thuenen_path = f'{temp_dump_fold}THUENEN.vrt'


    # year 
    year = 2019

    # paths
    path_to_S2_tiles = f'/data/Aldhani/eoagritwin/force/output/Guzinski/{year}/'

    ##### which tiles should be processed
    # get a list with all available tiles
    files = getFilelist(f'{path_to_S2_tiles}/tiles', '.tif', deep=True) 
    files = [file for file in files if any(tile in file for tile in tiles_to_process)]
    date_list = check_forceTSI_compositionDates(files)

    th_ds = gdal.Open(thuenen_path)
    th_arr = th_ds.GetRasterBand(1).ReadAsArray()
    mask = np.where(th_arr == -9999, 0, 1)

    lowRes_files = []
    highRes_files = []
    highRes_names = []

    colors = ['BLU', 'GRN', 'RED', 'NIR', 'RE1', 'RE2', 'RE3',  'SW1', 'SW2']

    #### S2 composites are time sensitive (need to be aligned with date of LST observation), so is incidence
    for date in date_list:
        # if not os.path.exists(f'{temp_dump_fold}INCIDENCE_{date}.vrt'):
        if date != '20190723': # '20190705':
            continue
        
        # check if for that date vrts were already processed
        if os.path.exists(f'{temp_dump_fold}S2_{date}'):
            print('S2 data for this date already processed')
            not_masked = [file for file in getFilelist(temp_dump_fold, '.vrt') if 'HIGHRES' in file]
            masked = [file for file in getFilelist(temp_dump_fold, '.tif') if 'HIGHRES' in file]
            highRes_files = [item for pair in zip(not_masked, masked) for item in pair]
            highRes_names = ["S2notMasked" if ".vrt" in file else "S2Masked" for file in highRes_files]
            lowRes_files = [f for f in [file for file in getFilelist(temp_dump_fold, '.tif') if 'Daily_LST' in file] for _ in range(2)]
            continue
        else:
            # get those tiles (and composite if more than one tile is provided)
            if len(tiles_to_process) == 1:
                tilesS2 = [file for file in getFilelist(path_to_S2_tiles, '.tif', deep=True) if tiles_to_process[0] in file and f'{date}.tif' in file]
                S2_path = f'{temp_dump_fold}S2_{date}.vrt'
                print(S2_path)
                vrt = gdal.BuildVRT(S2_path, tilesS2, separate=True)
                vrt = None
                vrt = gdal.Open(S2_path, gdal.GA_Update)  # VRT must be writable
                for idz, bname in enumerate(colors): 
                    band = vrt.GetRasterBand(1+idz)
                    band.SetDescription(bname)
                vrt = None

            else:
                tilesS2 = [file for file in getFilelist(path_to_S2_tiles, '.tif', deep=True) if any(tile in file for tile in tiles_to_process) and f'{date}.tif' in file] 
                force_to_vrt(tilesS2,
                        getCOLORSinOrderFORCELIST(tilesS2, colors, single=False),
                        f'{temp_dump_fold}S2_{date}',
                        False,
                        bandnames= colors)
    
                S2_path = [file for file in getFilelist(f'{temp_dump_fold}S2_{date}', '.vrt', deep=True) if '_Cube' in file][0]

                # determine LST and incidence files associated with respective S2 composite
            band_dict = transform_compositeDate_into_LSTbands(date, 4)


            # stat used for compositing
            for comp_stat in ['minVZA', 'maxLST']:
                path_to_incident = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE/{comp_stat}/{year}/'
                path_to_LST = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/LST_composites/{comp_stat}/{year}/'

                # get all LST bands that can be sharped with the S2 composite at this date (and sun angle incidence files as well, as they are dependent on that date
                LSTs = []

                for k, v in band_dict.items():
                    month = v['month']
                    band = int(v['band'])
                    v_path = f'{path_to_LST}Daily_LST_{comp_stat}_{year}_{month}.tif'
                    ds = gdal.Open(v_path, 0)

                    # export the LST for that day
                    LST_arr = ds.GetRasterBand(band).ReadAsArray() # store as single Tiff in temp
                    makeTif_np_to_matching_tif(LST_arr, v_path, f'{temp_dump_fold}Daily_LST_{comp_stat}_{year}_{month}_{band:02d}.tif')

                    # store the paths for selecting incidence for corresponding LST
                    incid_date = f'{year}_{month}_{band:02d}.tif'
                    lowRes_files.append(f'{temp_dump_fold}Daily_LST_{comp_stat}_{year}_{month}_{band:02d}.tif')

                    # incidence-tiles
                    incids = [file for file in getFilelist(path_to_incident, '.tif', deep=True) if any(tile in file for tile in tiles_to_process) and incid_date in file] 
                    # get those tiles (and composite if more than one tile is provided)
                    if len(tiles_to_process) == 1:
                        incid_path = incids[0]

                    else:
                        incid_path = f'{temp_dump_fold}INCIDENCE_{comp_stat}_{incid_date.split('.')[0]}.vrt'
                        gdal.BuildVRT(incid_path, incids)

                    # create highRes file through exapnding the vrt of S2
                    highRes_path = f'{temp_dump_fold}HIGHRES_{comp_stat}_{incid_date.split('.')[0]}.vrt'
                    gdal.BuildVRT(highRes_path, [S2_path, slope_path, aspect_path, incid_path], separate=True) # 
                    maskVRT_water(highRes_path) # maskVRT_water_and_drop_aux(highRes_path)
                    highRes_files.append(f'{highRes_path.split('.')[0]}_watermask.tif')
                    highRes_names.append('S2notMasked')
                    maskVRT(f'{highRes_path.split('.')[0]}_watermask.tif', mask)
                    highRes_files.append(f'{highRes_path.split('.')[0]}_watermask_S2_agromask.tif')
                    lowRes_files.append(f'{temp_dump_fold}Daily_LST_{comp_stat}_{year}_{month}_{band:02d}.tif')
                    highRes_names.append('S2Masked')



    joblist = []
    outFolder = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/sharpened2/allpred/{rand_foldname}/'
    for idx, highResFilename in enumerate(highRes_files):
        lowResFilename = lowRes_files[idx]
        f1 = f'{outFolder}{'/'.join(highResFilename.split('.')[0].split('_')[2:6])}/'
        for maskname, mask_lowRes in zip(['withoutLSTmask', 'withLSTmask'], ['', lowmask_bin_path]):
            lowResMaskFilename = mask_lowRes
            f2 = f'{f1}{maskname}/'
            for movWin in [15]:
                for cv in [0]:
                    for regrat in [0.25]:
                        kombi = f'mvwin{movWin}_cv{cv}_regrat{int(regrat*100):02d}_{highRes_names[idx]}_{maskname}'
                        f3 = f'{f2}{highRes_names[idx]}/'
                        os.makedirs(f'{f3}Residuals/', exist_ok=True)
                        os.makedirs(f'{f3}Values/', exist_ok=True)
                        joblist.append([highResFilename, 
                                    lowResFilename,
                                    lowResMaskFilename,
                                    cv, movWin, regrat,
                                    f'{f3}{'_'.join(highResFilename.split('.')[0].split('_')[2:6])}_{kombi}.tif'])

    print(f'\n{len(joblist)} times will be sharpened\n')


    if __name__ == '__main__':
        starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        print("--------------------------------------------------------")
        print("Starting process, time:" + starttime)
        print("")

        Parallel(n_jobs=ncores)(delayed(runSharpi)(job[0], job[1], job[2], job[3], job[4], job[5], job[6]) for job in joblist)


    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start : " + starttime)
    print("end: " + endtime)
    print("")