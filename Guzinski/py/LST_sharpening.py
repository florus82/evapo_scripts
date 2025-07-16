import sys
sys.path.append('/home/potzschf/repos/')
from helperToolz.guzinski import *
workhorse = True

if workhorse:
    origin = 'Aldhani/eoagritwin/'
else:
    origin = ''

# the sharpening is done in the follwoing order: Year, Tile, Day


############## define paths and parameters
base_path = f'/data/{origin}force/output/Guzinski/'
year = 2019

############## sanity checks
# first, all files for months outside of studay aim will be deleted to save storgage
reduce_forceTSI_output_to_validmonths(f'{base_path}{year}/tiles', 4, 10)
print('All FORCE files outside April-October deleted')
# get a list with all available tiles
files = getFilelist(f'{base_path}{year}/tiles', '.tif', deep=True) 
unique_tiles = get_forceTSI_output_Tiles(files)
print(f'There are {len(unique_tiles)} tiles available for processing for the year {year}')


# define the tiles to process and check if they contain composites for the same dates so that they can be stacked in the mosaiked
# (there is no sanity check for mosaiking)
tiles_to_process = unique_tiles#createFORCEtileLIST([67], [41])
date_list = []
for tile in tiles_to_process:
    date_list.append((get_forceTSI_output_DOYS([file for file in files if tile in file])))
for i in range(0,len(date_list)-1):
    if date_list[i] == date_list[i + 1]:
        continue
    else:
        print('the doys of composites across tiles is not equal - Better check!!!!!')
date_list = date_list[0]