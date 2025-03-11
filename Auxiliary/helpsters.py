import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import math
import os
import time
import xarray as xr 
from osgeo import ogr, osr
import random

# getFilelist returns a list with all files of a certain type in a path
def getFilelist(originpath, ftyp, deep = False):
    out   = []
    if deep == False:
        files = os.listdir(originpath)
        for i in files:
            if i.split('.')[-1] in ftyp:
                if originpath.endswith('/'):
                    out.append(originpath + i)
                else:
                    out.append(originpath + '/' + i)
            # else:
            #     print("non-matching file - {} - found".format(i.split('.')[-1]))
    else:
        for path, subdirs, files in os.walk(originpath):
            for i in files:
                if i.split('.')[-1] in ftyp:
                    out.append(os.path.join(path, i))
    return out

def plotter(array, row=1, col=1, names=False, title=False):

    # Plot the slices
    fig, axes = plt.subplots(row, col, figsize=(col*5, row*5), constrained_layout=False)  # 4 slices
    # Create a colormap
    cmap = plt.cm.viridis
    
    if col != 1:
        slice_indices = np.linspace(0, (row * col) -1, col * row, dtype=int)
        #print(slice_indices)
        for ax, idx in zip(axes.ravel(), slice_indices):
            im = ax.imshow(array[:, :, idx], cmap=cmap)
            if names == False:
                ax.set_title(f"Slice {idx}")
            else:
                ax.set_title(names[idx], fontsize=10)     
            # ax.set_xticks([0, 32, 64, 96, 127])
            # ax.set_yticks([0, 32, 64, 96, 127])
            # ax.set_xticklabels(['X0', 'X32', 'X64', 'X96', 'X127'])
            # ax.set_yticklabels(['Y0', 'Y32', 'Y64', 'Y96', 'Y127'])

            cbar_ax = ax.inset_axes([0.1, -0.2, 0.8, 0.05])  # [x, y, width, height]
            cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
            cbar.set_label('Value Scale')
    else:
            im = axes.imshow(array[:, :], cmap=cmap)
            #ax.set_xticks([0, 32, 64, 96, 127])
            #ax.set_yticks([0, 32, 64, 96, 127])
            #ax.set_xticklabels(['X0', 'X32', 'X64', 'X96', 'X127'])
            #ax.set_yticklabels(['Y0', 'Y32', 'Y64', 'Y96', 'Y127'])

            cbar_ax = axes.inset_axes([0.1, -0.2, 0.8, 0.05])  # [x, y, width, height]
            cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
            cbar.set_label('Value Scale')
 
    fig.subplots_adjust(hspace=0.2, wspace=0.2)
    if title != False:
        fig.suptitle(f'Date is {title[0]} at canals {title[1]}', fontsize=12)
    plt.tight_layout()
    plt.show()

def getNestedListMinLengthIndex(nestedList):
    res = [index for index, band in enumerate(nestedList) if len(band) == min([len(i) for i in nestedList])]
    return res[0]

def getBandNames(rasterstack):
    bands = []
    ds = gdal.Open(rasterstack)
    numberBands = ds.RasterCount
    for i in range(numberBands):
        bands.append(ds.GetRasterBand(i+1).GetDescription())
    return bands

def makeZeroNAN(arr):
    arr[arr == 0] = np.nan
    return arr


#####################################################################################
#####################################################################################
################# S3 downloads to numpy/tiff files ##################################
#####################################################################################
##################################################################################### 

def getAllDatesS3(listOfFiles, year='all'):
    '''Takes a list of paths of .nc files for Sentinel-3 if year == all, all paths are considered. 
    If a year is provided, the dates are only extracted for the corresponding yearand is returned 
    as well as a pathlist subsetted to this year. Expected filenaming convention: Germany_2017-01-01_2017-01-31.nc.
    The second 2017 is important here'''

    if year != 'all':
        listOfFiles = [file for file in listOfFiles if int(file.split('_')[-1][0:4]) == year]
    if type(listOfFiles) != list:
        listOfFiles = [listOfFiles]
    for e, file in enumerate(listOfFiles):
        dat = xr.open_dataset(file)
        if e == 0:
            tim = dat['t'].to_numpy()
        else:
            tim = np.concatenate((tim, dat['t'].to_numpy()))
    if year != 'all':
        return np.sort(tim), listOfFiles
    else:
        return np.sort(tim)

def getGeoTransFromNC(ncfile):
    '''Takes a path to an ncfile or an xarray_dataset and returns a tupel that can be used for gdal's SetGeotransform()'''
    
    if type(ncfile) == str:
        ncfile = xr.open_dataset(ncfile)
    upperLeft_X = np.min(ncfile.coords['x'])
    upperLeft_Y = np.max(ncfile.coords['y'])
    rotation = 0
    pixelWidth = ncfile.coords['x'][1] - ncfile.coords['x'][0]
    pixelHeight = -pixelWidth
    return (upperLeft_X, pixelWidth, rotation, upperLeft_Y, rotation, pixelHeight)

def getShapeFromNC(ncfile):
    '''Takes a path to an ncfile or an xarray_dataset and returns shape[1], shape[0], shape[2] ,comparable to np.array.shape()'''
    
    if type(ncfile) == str:
        ncfile = xr.open_dataset(ncfile)
    return len(ncfile.coords['x'].values), len(ncfile.coords['y'].values), ncfile[','.join(ncfile.data_vars.keys()).split(',')[-1]].shape[0]

def getDataFromNC(ncfile):
    '''Takes a path to an ncfile or an xarray_dataset and returns a 3D numpy array of the data'''
    
    if type(ncfile) == str:
        ncfile = xr.open_dataset(ncfile)
    arr = ncfile[','.join(ncfile.data_vars.keys()).split(',')[-1]].to_numpy()
    return np.swapaxes(np.swapaxes(arr, 0, 1), 1, 2)
    
def getCRS_WKTfromNC(ncfile):
    '''Takes a path to an ncfile or an xarray_dataset and returns coordinate sys as wkt'''
    
    if type(ncfile) == str:
        ncfile = xr.open_dataset(ncfile)
    return ncfile['crs'].attrs['crs_wkt']

def convertNCtoTIF(ncfile, storPath, fileName, accDT, make_uint16 = False, explode=False):
    '''Converts a filepath to an nc file or a .nc file to a .tif with option to store it UINT16 (Kelvin values are multiplied by 100 before decimals are cut off)'''
    
    gtiff_driver = gdal.GetDriverByName('GTiff')
    geoTrans = getGeoTransFromNC(ncfile)
    geoWKT = getCRS_WKTfromNC(ncfile)
    typi = gdal.GDT_Float64
    numberOfXpixels, numberOfYpixels, numberofbands = getShapeFromNC(ncfile)
    dat = getDataFromNC(ncfile)
    noDataVal = np.nan

    if make_uint16 == True:
        dat = dat * 100
        dat.astype(np.uint16)
        typi = gdal.GDT_UInt16
        fileName = fileName.split('.tif')[0] + '_UINT16.tif'
        noDataVal = 0
    
    if explode == False:
        out_ds = gtiff_driver.Create(storPath + fileName, numberOfXpixels, numberOfYpixels, numberofbands, typi)
        out_ds.SetGeoTransform(geoTrans)
        out_ds.SetProjection(geoWKT)

        for band in range(numberofbands):
            out_ds.GetRasterBand(band + 1).WriteArray(dat[:, :, band])
            out_ds.GetRasterBand(band + 1).SetNoDataValue(noDataVal)
            out_ds.GetRasterBand(band + 1).SetDescription(str(accDT[band]).split('.')[0])
        del out_ds
    # explode=True --> write single raster tifs for each accDateTime
    else:
        for band in range(numberofbands):
            if make_uint16 == False:
                out_ds = gtiff_driver.Create(storPath + str(accDT[band]).split('.')[0].replace(':', '_') + '.tif', numberOfXpixels, numberOfYpixels, 1, typi)
            else:
                out_ds = gtiff_driver.Create(storPath + str(accDT[band]).split('.')[0].replace(':', '_') + '_UINT16.tif', numberOfXpixels, numberOfYpixels, 1, typi)
            out_ds.SetGeoTransform(geoTrans)
            out_ds.SetProjection(geoWKT)
            out_ds.GetRasterBand(1).WriteArray(dat[:, :, band])
            out_ds.GetRasterBand(1).SetNoDataValue(noDataVal)
            out_ds.GetRasterBand(1).SetDescription(str(accDT[band]).split('.')[0])
            del out_ds

def getAccDateTimesByfilename(dicti, filename):
    """Finds the index of a filename in lookUp['filename'] and retrieves and return corresponding accDateTimes"""
    
    if filename in dicti["filename"]:
        index = dicti["filename"].index(filename)  # Get index
        accDateTimes = dicti["accDateTimes"][index]  # Get corresponding times
        return accDateTimes
    else:
        return print(f'Filename {filename} not found!')

def exportNCarrayDerivatesInt(ncfile, storPath, fileName, bandname, arr, make_uint16 = False):

    gtiff_driver = gdal.GetDriverByName('GTiff')
    numberOfXpixels, numberOfYpixels, numberofbands = getShapeFromNC(ncfile)

    typi = gdal.GDT_Float32
    if make_uint16 == True:
        typi = gdal.GDT_Int16

    out_ds = gtiff_driver.Create(storPath + fileName, numberOfXpixels, numberOfYpixels, 1, typi)
    out_ds.SetGeoTransform(getGeoTransFromNC(ncfile))
    out_ds.SetProjection(getCRS_WKTfromNC(ncfile))
    out_ds.GetRasterBand(1).WriteArray(arr)

    out_ds.GetRasterBand(1).SetDescription(bandname	)
    del out_ds
