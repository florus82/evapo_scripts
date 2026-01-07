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


# paths
lowmask_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/THUENEN_GER_LST_WARP.tif'
lowmask_bin_path = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/reprojected/THUENEN_GER_LST_WARP_BINARY.tif'
temp_dump = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/SharpEvap/Brandenburg/SecondShot_2ndyear/'
trash_path = f'{temp_dump}trash/'

path_to_slope = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/SLOPE/'
path_to_aspect = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/ASPECT/'
path_to_agro = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/THUENEN_2021/'

path_to_dem = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/DEM/'
path_to_lat = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LAT/'
path_to_lon = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/LON/'
path_to_acq = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/Acq_time/'
path_to_vza = '/data/Aldhani/eoagritwin/et/Sentinel3/VZA/comp/'
path_to_vaa = '/data/Aldhani/eoagritwin/et/Sentinel3/VAA/comp/'

path_to_force = '/data/Aldhani/eoagritwin/force/output/Guzinski'
path_to_inci = '/data/Aldhani/eoagritwin/et/Auxiliary/DEM/Force_Tiles/INCIDENCE2/'
path_to_lst = '/data/Aldhani/eoagritwin/et/Sentinel3/LST/LST_values/LST_composites/'

# get the Tiles for Brandenburg
bran = pd.read_csv('/data/Aldhani/eoagritwin/misc/state_tile_csv/clipped_grid_bran_tiles.csv')
tiles_to_process = createFORCEtileLIST(list(bran['Tile_X']),
                                        list(bran['Tile_Y']), True)

tiles_to_process.sort()
joblist = []

for idx, tile_to_process in enumerate(tiles_to_process):
    joblist.append([tile_to_process, temp_dump, path_to_slope, path_to_aspect, path_to_agro, path_to_force,
                    path_to_inci, path_to_lst, '20180401', '20180910', ['maxLST'], ['S2only'], 2, False,
                    path_to_dem, path_to_lat, path_to_lon, path_to_acq, path_to_vaa, path_to_vza])


# joblist.append([tiles_to_process[-37], temp_dump, path_to_slope, path_to_aspect, path_to_agro, path_to_force,
#                 path_to_inci, path_to_lst, '20190628', '20190707', ['maxLST'], ['S2only'], 2, True])

# joblist.append([tiles_to_process[-43], temp_dump, path_to_slope, path_to_aspect, path_to_agro, path_to_force,
#                 path_to_inci, path_to_lst, '20190628', '20190707', ['maxLST'], ['S2only'], 2, True])

print(f'\n{len(joblist)} tiles will be processed\n')


ncores = 26


if __name__ == '__main__':
    starttime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("Starting process, time:" + starttime)
    print("")

    Parallel(n_jobs=ncores)(delayed(Sharp_Evap)(job[0], job[1], job[2], job[3], job[4], job[5], job[6], job[7], job[8], job[9], job[10],
                                                job[11], job[12], job[13], job[14], job[15], job[16],job[17], job[18], job[19]) for job in joblist)

    print("")
    endtime = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    print("--------------------------------------------------------")
    print("--------------------------------------------------------")
    print("start : " + starttime)
    print("end: " + endtime)
    print("")


# Sharp_Evap(tile_to_process=tiles_to_process[-7],
#            storFolder=temp_dump,
#            path_to_slope=path_to_slope,
#            path_to_aspect=path_to_aspect,
#            path_to_agro=path_to_agro,
#            path_to_force=path_to_force,
#            path_to_inci=path_to_inci,
#            path_to_lst=path_to_lst,
#            time_start='20190401',
#            time_end='20190408',
#            compList=['maxLST'],
#            predList=['S2only'],
#            S2mask=3,
#            printEvapInter=True)