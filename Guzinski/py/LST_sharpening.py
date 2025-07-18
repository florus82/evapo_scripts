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
comp_stat = 'max' # the statistic the daily LST images are based on 
LST_path = '/data/Aldhani/eoagritwin/et/Sentinel3/tiffs/LST/daily_observations_all/'

############## files clean up & sanity checks
# first, all files for months outside of studay aim will be deleted to save storgage
reduce_forceTSI_output_to_validmonths(f'{base_path}{year}/tiles', 4, 10)
print('All FORCE files outside April-October deleted')

# get a list with all available tiles
files = getFilelist(f'{base_path}{year}/tiles', '.tif', deep=True) 
unique_tiles = get_forceTSI_output_Tiles(files)
print(f'There are {len(unique_tiles)} tiles available for processing for the year {year}')

# check if they contain composites for the same dates so that they can be stacked in the mosaiked
date_list = check_forceTSI_compositionDates(files)


############### define the tiles to process and check if they contain composites for the same dates so that they can be stacked in the mosaiked
# (there is no sanity check for mosaiking --> are tiles neighbouting each other)
tiles_to_process = createFORCEtileLIST([67], [41]) # unique_tiles

for tile_processing in tiles_to_process:
    for date in date_list:
        