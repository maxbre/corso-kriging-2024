library(terra)

# import vector RV
# it imports all geographic features 
# that are going to be transferred in the next objects

v_rv <- vect('./data/veneto.gpkg')

# create a void raster object
# equatore circa 40000 km, 
# 1 grado long = 40000/360 circa 111 km, 0.1 gradi circa 11 km
# alla latitudine italia 45 molto meno quasi la metà
r_void <- rast(v_rv, res=0.2) 
r_void

library(geosphere)
# returns distance in meters
distm(c(12.5, 45.5), c(12.6, 45.5), fun = distHaversine)
# shortest distance (geodesic, great circle distance) based on WGS84 ellipsoid
# https://en.wikipedia.org/wiki/Great-circle_distance
distGeo(c(12.5, 45.5), c(12.6, 45.5))


r_rv <- rasterize(v_rv, r_void)
values(r_rv)<- -999
names(r_rv) <- "blank"
plot(r_rv)
r_rv
# resample the temp object, imported from tif file
# to get the same spatial resolution
# resample from terra

# check if that format from terra is causing the problem with kriging

# turning a SpatRaster object to a SpatialGridDataFrame
# this is the trick
# https:/raster#https://stackoverflow.com/questions/70756463/turning-a-spatraster-object-to-a-spatialgriddataframe
# gstat do not understand SpatRaster object from terra package
# need to convert to SpatialGridDataFrame (or SpatialPixelsDataFrame) 
# by passing through the use of the raster package
# r_temp_rv <- as(r_temp_rv, "SpatialGridDataFrame")
