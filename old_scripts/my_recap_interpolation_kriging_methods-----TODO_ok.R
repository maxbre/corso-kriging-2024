library(tidyverse)
library(sp)
library(gstat)
library(terra)
library(raster)

# read df temp and dem elevation
sp_temp_dem_s <- read_csv('./data_input/df_temp_dem_sampling.csv')

# create spatial objects for analysis with gstat
# points
coordinates(sp_temp_dem_s) <- ~x+y
proj4string(sp_temp_dem_s) <- CRS("EPSG:4326")
is.projected(sp_temp_dem_s)

# read raster temp RV
# why using raster package see the next note
r_temp_rv <- raster('./data_input/r_temp_rv.tif')

# read raster dem RV
# why using raster package see the next note
r_dem_rv <- raster('./data_input/r_dem_rv.tif')

# create spatial object for analysis with gstat
#
# turning a SpatRaster object to a SpatialGridDataFrame
# this is the trick!
# https:/raster#https://stackoverflow.com/questions/70756463/turning-a-spatraster-object-to-a-spatialgriddataframe
# gstat do not understand SpatRaster object from terra package
# you need to convert to SpatialGridDataFrame (or SpatialPixelsDataFrame) 
# by passing through the use of the raster package

#sp_temp_rv <- as(r_temp_rv, "SpatialGridDataFrame")
#sp_temp_rv <- as(r_temp_rv, "SpatialPixelsDataFrame")
# here about the difference between SpatialGridDataFrame vs SpatialPixelDataFrame
# https://gis.stackexchange.com/questions/363109/sf-spatialgriddataframe-vs-spatialpixelsdataframe-processing-visualization-an
#gridded(sp_temp_rv) <- TRUE
#proj4string(sp_temp_rv) <- CRS("EPSG:4326")
# is.projected(sp_temp_rv)
# but this is the problem! name of columns x and y are different s1,s2
# names(as.data.frame(sp_temp_rv))
# leave out this approach, too much problems with it!


####### alternatively deali with grid csv (a much safer option, in a sense!)
## uncomment next lines eventually
## from here 
# read raster RV in the form of a df
sp_temp_rv <- read_csv('./data_input/r_temp_rv.csv')

# turn it to a spatial grid
# For raster data, gstat only understands formats from the older sp package and the newer stars package.
#coordinates(sp_temp_rv) <- ~x+y
#gridded(sp_temp_rv) <- TRUE
#proj4string(sp_temp_rv) <- CRS("EPSG:4326")
#is.projected(sp_temp_rv)
## up to here 

# eventually you can play with the resolution of tif
# and the results change, obviously!
# for reference see the specific part of the code 
# in the file about the "interpolation_deterministic_methods"

# read raster RV in the form of a df
sp_temp_rv <- read_csv('./data_input/r_temp_rv.csv')

# turn it to a spatial grid
# For raster data, gstat only understands formats from the older sp package and the newer stars package.
coordinates(sp_temp_rv) <- ~x+y
gridded(sp_temp_rv) <- TRUE
proj4string(sp_temp_rv) <- CRS("EPSG:4326")
is.projected(sp_temp_rv)

# ordinary kriging -----

# analysis of variograms
# experimental sample variogram
# here using the same parametrization by pirotti exercise

# here two important parameters to define
#width: subsequent distance intervals into which data point pairs are grouped for semivariance estimates
#cutoff: spatial separation distance up to which point pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three.

width <- 0.01 # pay attention: lon,lat because of original coordinates
cutoff <- 15

# define the experimental variogram
exp_v_ok <- variogram(temp~1,
               data = sp_temp_dem_s,
               width = width, 
               cutoff = cutoff)
plot(exp_v_ok)
#plot(exp_v_ok, plot.numbers = TRUE)

# show models
vgm()
show.vgms()

# set model variogram
# intial parameter are set by eye esitmation
# model, sill, range, nugget 

model <- "Gau"
psill <- 5
range <- 10
nugget <- NA

# model variogram
model_v_ok <- vgm(psill = psill, 
               model = model, 
               range = range, 
               nugget = nugget)


# fit the model to the experimental variogram
# least square fit
fit_v_ok <- fit.variogram(exp_v_ok, model_v_ok)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
plot(exp_v_ok, fit_v_ok)
#plot(exp_v_ok, fit_v_ok, plot.numbers = TRUE)

# The range is the distance after which the variogram levels off. 
# The physical meaning of the range is that pairs of points that are this distance or greater apart 
# are not spatially correlated !
# The sill is the total variance contribution, or the maximum variability between pairs of points.

# strange enough, here it is something weird about the use of the function krige()
# by not naming the arguments it works correctly
# otherwise sometimes it fails (not clear to me why)
# https://github.com/r-spatial/gstat/issues/69

ok_temp <- krige(temp~1,
                 sp_temp_dem_s,     #data
                 sp_temp_rv,    #newdata
                 model = fit_v_ok)


plot(ok_temp)

# how to export to tif?
# https://gis.stackexchange.com/questions/419758/export-kriging-plot-as-geotiff
# writeRaster(ok_temp$var1.pred, 
#             "./data_input/r_temp_rv_ok.tif",
#             overwrite = FALSE)

# here about how to create a tif directly from gsta object, easy!
# https://gis.stackexchange.com/questions/453042/exporting-interpolation-idw-ok-uk-prediction-from-gstat-in-ascii-format

# now transform SpatialGridDataFrame to raster 
# because you need to calculate the delta with 
# the "original" raster of temperature for rv
# here just selecting the predicted values by ordinary kriging
# not dealing with the variance! 
r_temp_rv_ok <- raster::raster(ok_temp['var1.pred'])
names(r_temp_rv_ok) <- "temp_pred_ok"
plot(r_temp_rv_ok)

# export to raster
raster::writeRaster(r_temp_rv_ok, "./data_input/r_temp_rv_ok.tif", overwrite=TRUE)

# calculate the difference in temp between rasters
delta_r_temp_rv_ok <- r_temp_rv - r_temp_rv_ok

# rename to something more meaningful
names(delta_r_temp_rv_ok) <- "delta_temp"

plot(delta_r_temp_rv_ok)

hist(delta_r_temp_rv_ok, breaks=50)
summary(delta_r_temp_rv_ok)
quantile(delta_r_temp_rv_ok$delta_temp[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_ok <- as.data.frame(delta_r_temp_rv_ok, xy=TRUE) %>% cbind(method="krig_ord")
# here eventually pay attention to NAs and omit them

# universal kriging -----

width <- 0.01 # pay attention: lon,lat because of original coordinates
cutoff <- 15

# define the experimental variogram
exp_v_uk <- variogram(temp~x+y,
                   data = sp_temp_dem_s,
                   width = width, 
                   cutoff = cutoff)
plot(exp_v_uk)
#plot(exp_v_uk, plot.numbers = TRUE)

# set model variogram
# intial parameter are set by eye esitmation
# model, sill, range, nugget 

model <- "Gau"
psill <- 5
range <- 10
nugget <- NA

# model variogram
model_v_uk <- vgm(psill = psill, 
               model = model, 
               range = range, 
               nugget = nugget)


# fit the model to the experimental variogram
# least square fit
fit_v_uk <- fit.variogram(exp_v_uk, model_v_uk)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
plot(exp_v_uk, fit_v_uk)
#plot(exp_v_uk, fit_v_uk, plot.numbers = TRUE)

uk_temp <- krige(temp~x+y,
                 sp_temp_dem_s, #data
                 sp_temp_rv,    #newdata
                 model = fit_v_uk)

plot(uk_temp)

# now transform SpatialPixelDataFrame to raster 
# because you need to calculate the delta with 
# the "original" raster of temperature for rv
# here just selecting the predicted values by ordinary kriging
# not dealing with the variance! 
r_temp_rv_uk <- raster::raster(uk_temp['var1.pred'])
names(r_temp_rv_uk) <- "temp_pred_uk"
plot(r_temp_rv_uk)

# export to raster
raster::writeRaster(r_temp_rv_uk, "./data_input/r_temp_rv_uk.tif", overwrite=TRUE)

# calculate the difference in temp between rasters
delta_r_temp_rv_uk <- r_temp_rv - r_temp_rv_uk

# rename to something more meaningful
names(delta_r_temp_rv_uk) <- "delta_temp"

plot(delta_r_temp_rv_uk)

hist(delta_r_temp_rv_uk, breaks=50)
summary(delta_r_temp_rv_uk)
quantile(delta_r_temp_rv_uk$delta_temp[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_uk <- as.data.frame(delta_r_temp_rv_uk, xy=TRUE) %>% cbind(method="krig_univ")
# here eventually pay attention to NAs and omit them

# kriging with external drift -----

# now because we are considering the elevation
# as an external drift for the temperature
# we need to deal with a raster with both variable variables temp and elevation

#join data from raster
sp_temp_dem_rv <- r_temp_rv %>% 
  as.data.frame(xy=TRUE) %>% 
  inner_join( r_dem_rv %>% as.data.frame(xy=TRUE)) 

# export to csv
sp_temp_dem_rv %>%  write_csv('./data_input/r_temp_dem_rv.csv')

# #### attention here, check again!
sp_temp_dem_rv <- na.omit(sp_temp_dem_rv)

# turn it to a spatial grid
# For raster data, gstat only understands formats from the older sp package and the newer stars package.
coordinates(sp_temp_dem_rv) <- ~x+y
gridded(sp_temp_dem_rv) <- TRUE
proj4string(sp_temp_dem_rv) <- CRS("EPSG:4326")
is.projected(sp_temp_dem_rv)

width <- 0.01 # pay attention: lon,lat because of original coordinates
cutoff <- 15

# define the experimental variogram
exp_v_ked <- variogram(temp~elev_m+x+y,
                      data = sp_temp_dem_rv,
                      width = width, 
                      cutoff = cutoff)

plot(exp_v_ked)
#plot(exp_v_ked, plot.numbers = TRUE)

# set model variogram
# intial parameter are set by eye esitmation
# model, sill, range, nugget 

model <- "Gau"
psill <- 5
range <- 10
nugget <- NA

# model variogram
model_v_ked <- vgm(psill = psill, 
                  model = model, 
                  range = range, 
                  nugget = nugget)


# fit the model to the experimental variogram
# least square fit
fit_v_ked <- fit.variogram(exp_v_ked, model_v_ked)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
plot(exp_v_ked, fit_v_ked)
#plot(exp_v_ked, fit_v_ked, plot.numbers = TRUE)

########## pay attention
# here elevation is missing on sp_temp_rv
# go back and import it!
ked_temp <- krige(temp~elev_m+x+y,
                  df_temp_elev, #data
                  sp_temp_rv,    #newdata
                  model = fit_v_ked)

plot(ked_temp)
