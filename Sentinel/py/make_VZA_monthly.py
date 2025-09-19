import sys
import time
sys.path.append('/home/potzschf/repos/')
from helperToolz.helpsters import *

# the raw nc files for VZA were downloaded for 15 days, not per month as LST
# for easier managing, files will be merged here

files = sorted(getFilelist('/data/Aldhani/eoagritwin/et/Sentinel3/VZA/raw', '.nc'))

for year in [i for i in range(2017,2025,1)]:
    outPath = f'/data/Aldhani/eoagritwin/et/Sentinel3/VZA/monthly_tiff_values/{year}/'
    os.makedirs(outPath, exist_ok=True) 

    # get a subset of files for that year
    yearFiles = [file for file in files if int(file.split('/')[-1].split('_')[-1][0:4]) == year]

    for month in [f'{i:02d}' for i in range(1,13)]:
        filename = f'VZA_{year}_{month}.tif'
        if os.path.exists(f'{outPath}{filename}'):
            t = time.localtime()
            ti = time.strftime("%H:%M:%S", t)
            print(f"already exists - next one at {ti}")
        else:
            print(f'working on month {month} in {year}')
            yfiles = [yF for yF in yearFiles if month == yF.split('-')[1]]

            monthL = []
            bnames = []
            for file in yfiles:
                monthL.append(getDataFromNC_VZA(file))
                accDateTimes = getAllDatesS3(file)
                for l in range(len(accDateTimes)):
                    bnames.append(str(accDateTimes[l].astype('datetime64[s]')))

            block = np.dstack(monthL)
            # block[np.isnan(block)] = 50
            # block2 = np.where(block > 45, 0, 1)
            # block3 = block2.astype(np.int8)

            exportNCarrayDerivatesInt(file, outPath,
                                    filename,
                                    bnames,
                                    block, #block3,
                                    #datType=gdal.GDT_Int8,
                                    numberOfBands=block.shape[2])