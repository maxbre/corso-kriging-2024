# load packages ----
library(terra)

# naming objects, for sure this is the most difficult part! -----
# general convention adopted here for naming objects
# template: "o_vars_geo"
# where:
# o = the object, e.g.: v for vectors, r for rasters, and so on...
# vars = the content of object, e.g.: temp, dem
# geo = the geographical boundaries, e.g.: rv for regione veneto, s for samplings


# import raw data ----
# read raster temperature
r_temp <- rast("./data_raw/VenetoCorrectedMODIS_LST_Avg2017.tif")
# read raster elevation (dem)
r_dem <- rast("data_raw/VenetoDEM.tif")

# check for different resolution
res(r_temp) == res(r_dem)

# rasters are with different resolution
# and this can cause a lot of headheaches
# in the further manipulation of the dataset
# and that's why of the next quite elaborated manipulation of data

# pivot raster -----
# trick: use a pivot void raster with the shape of RV
v_shp_rv <- vect("./data_raw/veneto.gpkg")
plot(v_shp_rv)

# next code is about modfying the resolution of the tow raster temp and dem
# by using a pivot model to transfer spatial information (by resample)
# it's a bit tricky (take some time to ponder!)

# create a void raster object
# pay attention about the resolution lon,lat
# here define the granularity of dataset
r_shp_rv <- rast(v_shp_rv, res = 0.02) # it is in decimal degrees because of lon,lat
r_shp_rv


# equator is about 40000 km,
# 1 degree long = 40000/360 ~ 111 km, 0.1 gradi ~ 11 km
# at latitude of IT a lot less, halve it, roughly!
# to know what does it mean 0.1 degree in longitude, here in rv (about)
# library(geosphere)
# distm(c(12.5, 45.5), c(12.6, 45.5), fun = distHaversine)
# shortest distance (geodesic, great circle distance) based on WGS84 ellipsoid
# https://en.wikipedia.org/wiki/Great-circle_distance
# distGeo(c(12.5, 45.5), c(12.6, 45.5))

# rasterize the vector shp rv
# this is going to be a sort of pivot raster to align all other raster objects
# transfer values associated with the geometries of vector data to a raster
r_shp_rv <- rasterize(v_shp_rv, r_shp_rv)
r_shp_rv <- classify(r_shp_rv, cbind(1, -999))
names(r_shp_rv) <- "shape_rv"
plot(r_shp_rv)

# export as void raster rv in the form of tif (to be eventually used later)
writeRaster(r_shp_rv, './data_input/r_shp_rv.tif', overwrite = TRUE)

# export as void raster rv in the form of csv
write.csv(as.data.frame(r_shp_rv, xy = TRUE),
          './data_input/r_shp_rv.csv',
          row.names = FALSE)

# resampling -------
# resample the temp object, imported from tif file
# to get the same spatial resolution by using the pivot raster shape of rv
r_temp_rv <- terra::resample(r_temp, r_shp_rv, method = "bilinear")
plot(r_temp_rv)

# change name to something more meaningful
names(r_temp_rv) <- "temp"

# crop and mask RV raster temp
r_temp_rv <- crop(r_temp_rv, r_shp_rv, mask = TRUE)
plot(r_temp_rv)

# export raster temp rv (to be later used again)
writeRaster(r_temp_rv, './data_input/r_temp_rv.tif', overwrite = TRUE)

# resample the dem object, imported from tif file
# to get the same spatial resolution
r_dem_rv <- terra::resample(r_dem, r_shp_rv, method = "bilinear")
plot(r_dem_rv)

# change name to something more meaningful
names(r_dem_rv) <- "elev_m"

# crop and mask RV raster dem
r_dem_rv <- crop(r_dem_rv, r_shp_rv, mask = TRUE)
plot(r_dem_rv)

# export raster dem rv (to be later used again)
writeRaster(r_dem_rv, './data_input/r_dem_rv.tif', overwrite = TRUE)

# export rasters as csv -----

library(tidyverse)

# export raster as df (to be later used again)
r_temp_rv %>%
  as.data.frame(xy = TRUE) %>%
  write_csv('./data_input/r_temp_rv.csv')

# export raster as df (to be later used again)
r_dem_rv %>%
  as.data.frame(xy = TRUE) %>%
  write_csv('./data_input/r_dem_rv.csv')

# get some samples from rasters -----

# sample both temp and dem rasters
# draw some sample points of temp
# it is a spatial vector
set.seed(123)
v_temp_rv_s <- spatSample(
  r_temp_rv,
  method = "random",
  replace = FALSE,
  # default with random
  size = 1000,
  na.rm = TRUE,
  as.points = TRUE,
  xy = TRUE
  )

plot(r_temp_rv)
plot(
  v_temp_rv_s,
  add = TRUE,
  col = "red",
  cex = 0.3,
  pch = 3
  )

# export temp samp vector as gpkg
# to be later reused
writeVector(v_temp_rv_s,
            './data_input/v_temp_rv_sampling.gpkg',
            overwrite = TRUE)

# extract elev sample from raster dem RV through temp points
# various way to accomplish this, see online help
# output is a df
df_dem_rv_s <- terra::extract(r_dem_rv,
                              v_temp_rv_s,
                              #df_temp_rv_s[c("x","y")],
                              xy = TRUE,
                              exact = TRUE,
                              ID = FALSE)

# transform spatial vector of temp to df
df_temp_rv_s <- as.data.frame(v_temp_rv_s, geom = "XY")

# inner join temp and dem sampling as df
# this seems a more safer approach because is checking for same exact coordinates
df_temp_dem_rv_s <- inner_join(df_temp_rv_s, df_dem_rv_s) %>%
  select(x, y, temp, elev_m)

# check for any NA
sum(is.na(df_temp_dem_rv_s))

# where are NAs
which(is.na(df_temp_dem_rv_s), arr.ind = TRUE)

# filter NAs
# https://stackoverflow.com/questions/74525000/how-to-filter-for-rows-containing-na
# tidyverse solution (it seems slightly more convoluted than base approach)
# df_temp_elev <- df_temp_dem_s %>%
# filter(!if_any(everything(), is.na))

# base solution (better, more concise and straight, I sort of prefer it!)
# not needed
df_temp_dem_rv_s <- df_temp_dem_rv_s[complete.cases(df_temp_dem_rv_s), ]

# export dataset temp and dem for sampling data
write_csv(df_temp_dem_rv_s,
          './data_input/df_temp_dem_rv_sampling.csv')

# export rasters for rv temp and dem as csv -------

df_temp_rv <- r_temp_rv %>%
  as.data.frame(xy = TRUE)

df_dem_rv <- r_dem_rv %>%
  as.data.frame(xy = TRUE)

# check resolution
res(r_temp_rv) == res(r_dem_rv)

# here I use left join
# because of different number of rows
# between temp and dem (not comletely clear why...)
df_temp_dem_rv <- left_join(df_temp_rv, df_dem_rv)

# check for any NA
sum(is.na(df_temp_dem_rv))

# where are NAs?
i <- which(is.na(df_temp_dem_rv), arr.ind = TRUE)
# 12.3881 44.79701

# select rows
r <- i[, 1]

# this is a bit rough, but no time now to clean up too much things...
lon <- df_temp_dem_rv[r, 1]
lat <- df_temp_dem_rv[r, 2]

# so now the big question is:
#.... but where the heck are located those missing pixels?

pixel <- data.frame(lon = lon, lat = lat)

# convert to an sf object
pixel <- sf::st_as_sf(pixel, coords = c("lon", "lat"), crs = 4326) # lat/long coordinate reference system
                      
# just take a look at it by using mapview
mapview::mapview(pixel)

# missing data in a raster, how to impute the value?
# considering the position zero elevation is almost fine

# now fix it by hand and close this sort of nightmare!
df_temp_dem_rv[r, 4] <- 0

# for raster
# https://stackoverflow.com/questions/57214163/extract-raster-value-based-on-list-of-coordinates-sptransform

# export dataset for sampling data
write_csv(df_temp_dem_rv, './data_input/r_temp_dem_rv.csv')

# closing up and exporting as rda -----
#save objects to rda file
#save(xx, xxx, xxxxx, file = "./data_input/00_objects.rda")

#library(tidyverse)
#write_rds()
