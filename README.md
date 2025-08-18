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
sharpen LST (LST ~ S2,slope,aspect,incidence) # WE NEED MASKS for LC that we dont WANT!!!!!!!


## there is really an issue with the compositing of LST values:
    - missing pixels in LST lead to holes in acquisiton time
    - this leads to holes in incidence
    - this leads to holes in sharpened images
    - also, it leads era5 holes as raw era5 is warped to lst acquisition times for value extraction at satellite overpass (LST scene time)


## Efficiency considerations
  --> the entire data set? makes sense, if more than one LST day is processed in a loop (warp_at_doy needs to implement that loop)
  --> only the bands needed for that specific day could be warped if only one (or a few days) of a month are processed (then warp to reference needed     
                    implementation, as well as warp_at_doy needed adaptation to a) get the bands b) loop differently over bands???)
  --> not load all the tifs (slope, aspect, lat, lon, dem) for germany, but use the .vrt files created during the sharpening process
          --> does not save a lot, as they are only needed for the correction of ssrd



# get_warped_ERA5_at_doy
 - output only for warped, not doy so far
 - when sharping to 


# ein model pro lst observation und dann mit allen (?) den Tag vorhersagen: immer noch das Problem der Uhrzeiten --> vielleicht ein min/max Uhrzeit threshold?
# incidence raus und testen: um die LST Pixel hoffentlich im sharpened image um eliminieren, bzw. zu minimieren (da ja das compositing Problem noch besteht)
# anderes compositing --> bis jetzt minVZA und maxLST (+o°C Threshold, und) --> das vorhandene compositing nochmal richtig prüfen!!!
# ssrd tuning --> zum einen das warpen in den Griff bekommen, zum anderen die Values checken --> formeln nochmals richtig durchgehen
