library(terra)

# read raster temperature RV
temp_r <- rast("data/VenetoCorrectedMODIS_LST_Avg2017.tif")

plot(temp_r)

# read raster elevation (dem) RV
dem_r <- rast("data/VenetoDEM.tif")

## temp sample 1000 points from raster temp RV
temp_s <- spatSample(temp_r,
                     method="random",
                     size=1000,
                     na.rm=T,
                     as.points=T)
# rename variable
names(temp_s) <- "temp"

# extract elev sample from raster dem RV through temp points
elev_s <- extract(dem_r, temp_s, ID= FALSE)

# rename variable
names(elev_s) <- c("elev_m")

# extract coordinates
xy <- geom(temp_s)[,3:4]

# transform to each object to a df
temp <- as.data.frame(temp_s)
elev <- as.data.frame(elev_s)

#compose final df to work with by leaving out na
mydf <- na.omit(data.frame(xy, temp, elev))


########## grid interpolazione

# read gpkg as vector 
RV_v <- vect("./data/veneto.gpkg")

# also possible this approach
#library(sf)
#RV_sf <- read_sf("./data/veneto.gpkg")

# transform vector to raster
# data are in lat long!
RV_r <- rast(RV_v, resolution = 0.01)

# add some values to the raster
values(RV_r) <- -999

# set the name of field

names(RV_r) <- 'blank'
# crop the raster to the shape of region
RV_r <- crop(RV_r, RV_v, mask=TRUE)


####################################################################
# this

RV_r_temp <- crop(temp_r, RV_v, mask=TRUE)
plot(RV_r_temp)


temp_samp <- spatSample(RV_r_temp,
                     method="random",
                     size=500,
                     na.rm=T,
                     as.points=T)

plot(RV_r_temp)
plot(temp_samp, add=TRUE, col="white", cex=0.4, pch=3)

####################################################################

# plot 
plot(RV_r)

############### interpolazione metodo deterministico

# interpolate nearest neighbour
RV_r_nn <- interpNear(RV_r, temp_s, field="temp", radius=1)
# mask
RV_r_nn <- mask(RV_r_nn, RV_r)

temp_s_mask <- mask(temp_s, RV_v)

plot(RV_r_nn)
plot(temp_s_mask, add=T, col="white", cex=0.4, pch=3)

plot(RV_r-RV_r_nn)
