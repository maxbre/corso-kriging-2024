library(terra)
library(tidyverse)
library(viridis)

# read gpkg as vector of RV
v_rv <- vect("./data/veneto.gpkg")
plot(v_rv)

# read raster temperature
r_temp <- rast("./data/VenetoCorrectedMODIS_LST_Avg2017.tif")
# change name to something more meaningful
names(r_temp) <- "temp"
# crop and mask RV raster temp
r_temp_rv <- crop(r_temp, v_rv, mask=TRUE)
plot(r_temp_rv)

# read raster elevation (dem)
r_dem <- rast("data/VenetoDEM.tif")
# change name to something more meaningful
names(r_dem) <- "elev_m"
# crop and mask RV raster temp
r_dem_rv <- crop(r_dem, v_rv, mask=TRUE)
plot(r_dem_rv)

# check for resolution of rasters
res(r_temp_rv) == res(r_dem_rv)

# then resample
# transfer values of a SpatRaster to another one with a different geometry
# here the geometry is resampled from r_dem_rv to r_temp_rv
# it can be also the other way round
r_temp_rv <- resample(r_temp_rv, r_dem_rv)

# now check gain for resolution of rasters
res(r_temp_rv) == res(r_dem_rv)

# about the problem of aligning rasters with different resolution
# you have to to resample! check for answer by Hijmans (the author of terra package)
# https://stackoverflow.com/questions/71925095/convert-raster-from-coarse-to-finer-resolution
# https://stackoverflow.com/questions/32278825/how-to-change-the-resolution-of-a-raster-layer-in-r/70115798#70115798

# export raster temp rv (to be later used again)
writeRaster(r_temp_rv, './data_input/r_temp_rv.tif', overwrite = TRUE)

# export raster as df (to be later used again)
r_temp_rv %>% 
  as.data.frame(xy=TRUE) %>% 
  write_csv('./data_input/r_temp_rv.csv')

# export raster dem rv (to be later used again)
writeRaster(r_dem_rv, './data_input/r_dem_rv.tif', overwrite = TRUE)

# export raster as df (to be later used again)
r_dem_rv %>% 
  as.data.frame(xy=TRUE) %>% 
  write_csv('./data_input/r_dem_rv.csv')

# sample both temp and dem rasters
# draw some sample points of temp
# it is a spatial vector
set.seed(123)
v_temp_rv_s <- spatSample(r_temp_rv,
                          method = "random",
                          size = 1000,
                          na.rm = TRUE,
                          as.points = TRUE,
                          xy = TRUE)
plot(r_temp_rv)
plot(v_temp_rv_s, add=TRUE, col="red", cex=0.3, pch=3)

# transform spatial vector to df
df_temp_s <- as.data.frame(v_temp_rv_s, geom="XY")

# extract elev sample from raster dem RV through temp points
# various way to accomplish this, see online help
df_dem_s <- terra::extract(r_dem_rv, 
                           df_temp_s[c("x","y")], 
                           xy=TRUE, 
                           exact=TRUE, 
                           ID=FALSE)

# inner join
# this seems a more safer approach because is checking for same exact coordinates
df_temp_elev <- inner_join(df_temp_s, df_dem_s) %>% 
  select(x, y, temp, elev_m)

# check for any NA
sum(is.na(df_temp_elev))

# where are NAs
which(is.na(df_temp_elev), arr.ind=TRUE)

# filter NAs
# https://stackoverflow.com/questions/74525000/how-to-filter-for-rows-containing-na

# tidyverse solution (slightly more convoluted)
#df_temp_elev <- df_temp_elev %>%
#  filter(!if_any(everything(), is.na))

# base solution (better, more concise, I prefer it!)
df_temp_elev <- df_temp_elev[complete.cases(df_temp_elev),]

# export dataset
write_csv(df_temp_elev, './data_input/df_temp_elev.csv')

##################### end of data munging

# interpolation deterministic methods ----

# nn ------
# nearest neighbour (plain, no linear interpolation between points)
# interpolation between points is FALSE by defualt  
r_temp_rv_nn <- interpNear(r_temp_rv,
                           v_temp_rv_s,
                           field="temp",
                           radius=1,
                           interpolate=FALSE)


# mask
r_temp_rv_nn <- mask(r_temp_rv_nn, r_temp_rv)

# plot masked
plot(r_temp_rv_nn)

# calculate delta
delta_r_temp_rv_nn <- r_temp_rv_nn - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_nn) <- "delta_temp"

plot(delta_r_temp_rv_nn)

hist(delta_r_temp_rv_nn, breaks=50)
summary(delta_r_temp_rv_nn)
quantile(delta_r_temp_rv_nn$delta_temp[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_nn <- as.data.frame(delta_r_temp_rv_nn, xy=TRUE) %>% cbind(method="near_neigh")

# nn_ iln ------
# nearest neighbour with linear interpolation between points
# interpolation between points is FALSE by defualt  
# note here interpolate = TRUE
r_temp_rv_nn_iln <- interpNear(r_temp_rv,
                               v_temp_rv_s,
                               field = "temp",
                               radius = 1,
                               interpolate = TRUE)

# mask
r_temp_rv_nn_iln <- mask(r_temp_rv_nn_iln, r_temp_rv)

# plot masked
plot(r_temp_rv_nn_iln)

# calculate delta
delta_r_temp_rv_nn_iln <- r_temp_rv_nn_iln - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_nn_iln) <- "delta_temp"

plot(delta_r_temp_rv_nn_iln)

hist(delta_r_temp_rv_nn_iln, breaks=50)
summary(delta_r_temp_rv_nn_iln)
quantile(delta_r_temp_rv_nn_iln$delta_temp[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_nn_iln <- as.data.frame(delta_r_temp_rv_nn_iln, xy=TRUE) %>% cbind(method="near_neigh_iln")

# idw pw1 -----
# inverse distance weighting 
# power 1

r_temp_rv_idwp1 <- interpIDW(r_temp_rv, 
                      v_temp_rv_s, 
                      field="temp", 
                      radius=1, 
                      power=1)

# mask
r_temp_rv_idwp1 <- mask(r_temp_rv_idwp1, r_temp_rv)

# plot mask
plot(r_temp_rv_idwp1)

# calculate delta
delta_r_temp_rv_idwp1 <- r_temp_rv_idwp1 - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_idwp1) <- "delta_temp"

plot(delta_r_temp_rv_idwp1)

hist(delta_r_temp_rv_idwp1, breaks=50)
summary(delta_r_temp_rv_idwp1)
quantile(delta_r_temp_rv_idwp1$delta_temp[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_idw_p1 <- as.data.frame(delta_r_temp_rv_idwp1, xy=TRUE) %>% cbind(method="idw_pw1")

# idw pw2 -----
# inverse distance weighting 
# power 2

r_temp_rv_idwp2 <- interpIDW(r_temp_rv, 
                             v_temp_rv_s, 
                             field = "temp", 
                             radius = 1, 
                             power = 2)

# mask
r_temp_rv_idwp2 <- mask(r_temp_rv_idwp2, r_temp_rv)

# plot mask
plot(r_temp_rv_idwp2)

# calculate delta
delta_r_temp_rv_idwp2 <- r_temp_rv_idwp2 - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_idwp2) <- "delta_temp"

plot(delta_r_temp_rv_idwp2)

hist(delta_r_temp_rv_idwp2, breaks=50)
summary(delta_r_temp_rv_idwp2)
quantile(delta_r_temp_rv_idwp2$delta_temp[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_idw_p2 <- as.data.frame(delta_r_temp_rv_idwp2, xy=TRUE) %>% cbind(method="idw_pw2")

### closing up...
# bind rows all data frames ----
df_delta_temp <- bind_rows(df_dt_nn,
                           df_dt_nn_iln,
                           df_dt_idw_p1,
                           df_dt_idw_p2)


# plot of deltas for different interpolation approaches

# boxplots
df_delta_temp %>% 
  ggplot() + 
  geom_boxplot(aes(x=method, y=delta_temp, colour = method))+
  theme_minimal()

# maps of deltas
df_delta_temp %>% 
  ggplot()+
  geom_raster(aes(x=x, y=y, fill = delta_temp))+
  scale_fill_viridis(option = "plasma")+
  facet_wrap(vars(method))+
  coord_sf(default_crs = sf::st_crs(4326))+
  theme_void()

