download LST
download VZA (could only be downloaded at 14 day intervals; make_VZA_monthly.ipynb merges them to monthly .nc files for easier comparison with monthly LST)
download ERA5
download DEM

preprocess DEM -> slope & aspect # there is an one-pixel too much issue with the slope/aspect/incidence per tile, because of the stupid shapefile 
                     --> report on github and just cut off last column
create lat/lon raster
airTemp -> DEM sharpen at LST spatial resolution 
ssrdcs -> correct for terrain etc. at LST spatial resolution # maybe do on-the fly only for bands(timesteps) needed???
preprocess LST (airTemp, VZA, > 0°C) minVZA vs maxLST

calculate Incidence (slope, aspect, acquisition time LST, lat/lon)
prepare FORCE S2
sharpen LST (LST ~ S2,slope,aspect,incidence)


## there is really an issue with the compositing of LST values:
    - missing pixels in LST lead to holes in acquisiton time
    - this leads to holes in incidence
    - this leads to holes in sharpened images
    - also, it leads era5 holes as raw era5 is warped to lst acquisition times for value extraction at satellite overpass (LST scene time)