import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *
from helperToolz.dicts_and_lists import *
from helperToolz.guzinski import * 
from other_repos.pyDMS.pyDMS.pyDMS import *
import time

os.environ["GDAL_MAX_DATASET_POOL_SIZE"] = "600"
# get tiles for Brandenburg
bran = pd.read_csv('/data/Aldhani/eoagritwin/misc/state_tile_csv/clipped_grid_bran_tiles.csv')

# create vrts of slope and aspect
path_to_slope = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/'
path_to_aspect = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/'
temp_dump = '/data/Aldhani/eoagritwin/et/Sentinel3//LST/LST_values/tempDump/'


# get all highRes datasets: S2 composites, aspect, ratio and incedence for the tile/time and stack them (in a composite if more than one tile is processed)
# set the highRes S2 tiles that will be processed
tiles_to_process = createFORCEtileLIST(list(bran['Tile_X']),
                                        list(bran['Tile_Y']))

tiles_to_process = createFORCEtileLIST([58, 59, 58, 59],
                                       [33, 34, 34, 34])

# tiles_to_process = get_forceTSI_output_Tiles(getFilelist(path_to_slope, '.tif'))

# slope-tiles
slopes = [file for file in getFilelist(path_to_slope, '.tif') if any(tile in file for tile in tiles_to_process)] # if any tile name is in file
# aspect-tiles
aspects = [file for file in getFilelist(path_to_aspect, '.tif') if any(tile in file for tile in tiles_to_process)] # if any tile name is in file

# get those tiles (and composite if more than one tile is provided)
if len(tiles_to_process) == 1:
    
    slope_path = slopes[0]
    aspect_path = aspects[0]

else:
    slope_path = f'{temp_dump}SLOPE.vrt'
    gdal.BuildVRT(slope_path, slopes)

    aspect_path = f'{temp_dump}ASPECT.vrt'
    gdal.BuildVRT(aspect_path, aspects)

# year 
for year in [2019]:#range(2017,2025):
  
    # paths
    path_to_S2_tiles = f'/data/Aldhani/eoagritwin/force/output/Guzinski/{year}/'
    


    ##### which tiles should be processed
    # get a list with all available tiles
    files = getFilelist(f'{path_to_S2_tiles}/tiles', '.tif', deep=True) 
    files = [file for file in files if any(tile in file for tile in tiles_to_process)]
    date_list = check_forceTSI_compositionDates(files)


    #### S2 composites are time sensitive (need to be aligned with date of LST observation), so is incidence

    for date in date_list:
    # date = date_list[0] # will be replaced through loop
        if date != '20190705':
            continue
        # get those tiles (and composite if more than one tile is provided)
        if len(tiles_to_process) == 1:

            tilesS2 = [file for file in getFilelist(path_to_S2_tiles, '.tif', deep=True) if tiles_to_process[0] in file and f'{date}.tif' in file]
            S2_path = f'{temp_dump}S2_{date}.vrt'
            gdal.BuildVRT(S2_path, tilesS2)

        else:
            tilesS2 = [file for file in getFilelist(path_to_S2_tiles, '.tif', deep=True) if any(tile in file for tile in tiles_to_process) and f'{date}.tif' in file] 
            force_to_vrt(tilesS2,
                    getCOLORSinOrderFORCELIST(tilesS2, list(dict.fromkeys(tile.split('SEN2L_')[-1].split('_TSI')[0] for tile in tilesS2)), single=False),
                    f'{temp_dump}S2_{date}',
                    False,
                    bandnames= list(dict.fromkeys(tile.split('SEN2L_')[-1].split('_TSI')[0] for tile in tilesS2)))
            print(list(dict.fromkeys(tile.split('SEN2L_')[-1].split('_TSI')[0] for tile in tilesS2)))
            S2_path = [file for file in getFilelist(f'{temp_dump}S2_{date}', '.vrt', deep=True) if '_Cube' in file][0]
            
            # determine LST and incidence files associated with respective S2 composite
        band_dict = transform_compositeDate_into_LSTbands(date, 4)


        # stat used for compositing
        for comp_stat in ['minVZA', 'maxLST']:
            path_to_incident = f'/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE/{comp_stat}/{year}/'
            path_to_LST = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/{comp_stat}/{year}/'

            # get all LST bands that can be sharped with the S2 composite at this date (and sun angle incidence files as well, as they are dependent on that date
            incid_dates = []
            LSTs = []

            for k, v in band_dict.items():
                month = v['month']
                band = int(v['band'])
                v_path = f'{path_to_LST}Daily_LST_{comp_stat}_{year}_{month}.tif'
                ds = gdal.Open(v_path, 0)
                LST_arr = ds.GetRasterBand(band).ReadAsArray() # store as single Tiff in temp
                makeTif_np_to_matching_tif(LST_arr, v_path, f'{temp_dump}Daily_LST_{comp_stat}_{year}_{month}_{band:02d}.tif')

                # store the paths for selecting incidence for corresponding LST
                incid_dates.append(f'{year}_{month}_{band:02d}.tif')
                LSTs.append(f'{temp_dump}Daily_LST_{comp_stat}_{year}_{month}_{band:02d}.tif')
                ##### loop over the LST files and go

            for idx, LST_file in enumerate(LSTs):
                # incidence-tiles
                incids = [file for file in getFilelist(path_to_incident, '.tif', deep=True) if any(tile in file for tile in tiles_to_process) and incid_dates[idx] in file] 
                
                # get those tiles (and composite if more than one tile is provided)
                if len(tiles_to_process) == 1:
                    incid_path = incids[0]

                else:
                    incid_path = f'{temp_dump}INCIDENCE.vrt'
                    gdal.BuildVRT(incid_path, incids)

                # sanity check for incidence and LST date
                if (LSTs[idx].split(f'{year}')[-1] == incids[0].split(f'{year}')[-1]):
                    
                    # get LST file
                    lowRes_path = LSTs[idx]
                    # create highRes file through exapnding the vrt of S2
                    highRes_path = f'{temp_dump}HIGHRES.vrt'
                    gdal.BuildVRT(highRes_path, [S2_path, slope_path, aspect_path, incid_path], separate=True)
                else:
                    raise ValueError('Something is seriously wrong with the alignment of LST and incidence dates!!!!!')
                
                movWin = 40
                cv = 25
                highResFilename = highRes_path
                lowResFilename = lowRes_path
                outputFilename = f'/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/sharpened/germany/{comp_stat}_mvwin{movWin}_cv{cv}_{lowRes_path.split('_LST_')[-1]}'


                useDecisionTree = True

                commonOpts = {"highResFiles":               [highResFilename],
                                "lowResFiles":              [lowResFilename],
                                "lowResQualityFiles":         [],# [lowResMaskFilename],
                                "lowResGoodQualityFlags":     [],#[255],
                                "cvHomogeneityThreshold":     cv,
                                "movingWindowSize":           movWin,
                                "disaggregatingTemperature":  True}
                dtOpts =     {"perLeafLinearRegression":    True,
                                "linearRegressionExtrapolationRatio": 0.25}
                sknnOpts =   {'hidden_layer_sizes':         (10,),
                                'activation':                 'tanh'}
                nnOpts =     {"regressionType":             REG_sklearn_ann,
                                "regressorOpt":               sknnOpts}

                start_time = time.time()

                if useDecisionTree:
                    opts = commonOpts.copy()
                    opts.update(dtOpts)
                    disaggregator = DecisionTreeSharpener(**opts)
                else:
                    opts = commonOpts.copy()
                    opts.update(nnOpts)
                    disaggregator = NeuralNetworkSharpener(**opts)

                print("Training regressor...")
                disaggregator.trainSharpener()
                print("Sharpening...")
                downscaledFile = disaggregator.applySharpener(highResFilename, lowResFilename)
                print("Residual analysis...")
                residualImage, correctedImage = disaggregator.residualAnalysis(downscaledFile, lowResFilename,
                                                                            # lowResMaskFilename,
                                                                                doCorrection=True)
                print("Saving output...")
                highResFile = gdal.Open(highResFilename)
                if correctedImage is not None:
                    outImage = correctedImage
                else:
                    outImage = downscaledFile
                # outData = utils.binomialSmoother(outData)
                outFile = utils.saveImg(outImage.GetRasterBand(1).ReadAsArray(),
                                        outImage.GetGeoTransform(),
                                        outImage.GetProjection(),
                                        outputFilename)
                residualFile = utils.saveImg(residualImage.GetRasterBand(1).ReadAsArray(),
                                            residualImage.GetGeoTransform(),
                                            residualImage.GetProjection(),
                                            os.path.splitext(outputFilename)[0] + "_residual" +
                                            os.path.splitext(outputFilename)[1])

                outFile = None
                residualFile = None
                downsaceldFile = None
                highResFile = None

                print(time.time() - start_time, "seconds")

