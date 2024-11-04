###### ------ prepare dataset as data.frame ----

# read raster temperature RV
temp_r <- terra::rast("data/VenetoCorrectedMODIS_LST_Avg2017.tif")

# read raster elevation (dem) RV
dem_r <- terra::rast("data/VenetoDEM.tif")

## temp sample 1000 points from raster temp RV
temp_s <- terra::spatSample(temp_r,
                            method="random",
                            size=1000,
                            na.rm=T,
                            as.points=T)
# rename variable
names(temp_s) <- "temp"

# extract elev sample from raster dem RV through temp points
elev_s <- terra::extract(dem_r, temp_s, ID= FALSE)

# rename variable
names(elev_s) <- c("elev_m")

# extract coordinates
xy <- terra::geom(temp_s)[,3:4]

# transform to each object to a df
temp <- as.data.frame(temp_s)
elev <- as.data.frame(elev_s)

#compose final df to work with by leaving out na
mydf <- na.omit(data.frame(xy, temp, elev))

##############----- variogram ----

library(sp)
library(gstat)

coordinates(mydf) <- ~x+y
proj4string(mydf)<-CRS("epsg:4326")

# lagged scatterplot for first 300 rows
hscat(temp~1, data= mydf[1:300,], breaks=c(0, 0.05, 0.1, 0.5, 1, 2, 3))

############# case 1 ----
# variogram cloud for temp, selection of first 300 rows
plot(variogram(temp~1, data= mydf[1:300,], cloud=TRUE))

# experimental variogram for temp, with number of obs
p1 <- plot(variogram(temp~1, data= mydf), plot.numbers=TRUE)

# this is to just to update the size of labels
# I know it's a strange way, this is just a workaround....
# because the size of labels is actually hardcoded in the function
# but do not focus too much on this issue, it's not important now...
opts <- lattice::trellis.par.get()
opts$add.text$cex <- 0.6
update(p1, par.settings = opts)

############# case 2 ----
# experimental variogram for temp and covariate elev, with numbers of obs
p2 <- plot(variogram(temp~elev_m, data= mydf), plot.numbers=TRUE)

# same as before for about the labels
update(p2, par.settings = opts)

############# case 3 ----
# experimental variogram for temp and covariates x+y, with numbers of obs
p3 <- plot(variogram(temp~x+y, data= mydf), plot.numbers=TRUE)

# same as before for about the labels
update(p3, par.settings = opts)
