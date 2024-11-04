# load packages ----
library(tidyverse)
library(gstat)
library(sp)
library(gridExtra) # arrange plots
#library(terra)
#library(raster)

# naming objects, for sure the most difficult part! -----
# general convention adopted here for naming objects
# "o_vars_geo"
# o = the object, e.g.: v for r vector, r for rasters, and so on...
# vars = the content of object, e.g.: temp, dem
# geo = the geographical boundaries, e.g.: rv for regione veneto, s for sampling points


# import data_input ----

# read df temp and dem elevation
sp_temp_dem_s <- read_csv('./data_input/df_temp_dem_rv_sampling.csv')

# create spatial objects for analysis with gstat
# sampling points
coordinates(sp_temp_dem_s) <- ~x+y
proj4string(sp_temp_dem_s) <- CRS("EPSG:4326")
is.projected(sp_temp_dem_s)

# read temp and dem raster for RV 
sp_temp_dem_rv <- read_csv('./data_input/r_temp_dem_rv.csv')

# turn it to a spatial grid
# For raster data, gstat only understands formats from the older sp package and the newer stars package.
coordinates(sp_temp_dem_rv) <- ~x+y
gridded(sp_temp_dem_rv) <- TRUE
proj4string(sp_temp_dem_rv) <- CRS("EPSG:4326")
is.projected(sp_temp_dem_rv)

# co-kriging -----

# For modeling of Cross-Variogram 
# we have to build the gstat model sequentially
# by using the gstat method. 
# First we have to build a gstat structure for target variable (temp) and covariates (elev).
# Then, we will add fit variogram models to the gstat object.

## direct variogram target variable temp ---

variogram(temp ~ 1,
          data = sp_temp_dem_s,
          width = 0.1,
          cutoff = 10) -> exp_v_temp

plot(exp_v_temp)

automap::autofitVariogram(temp~1, sp_temp_dem_s, 
                          width = 0.1, 
                          cutoff = 10
                          )$var_model

model <- "Ste"
psill <- 2.5
range <- 8
kappa <- 0.6

# model variogram
model_v_temp <- vgm(psill = psill, 
                    model = model, 
                    range = range,
                    kappa = kappa)

# fit the model to the experimental variogram
# least square fit
fit_v_temp <- fit.variogram(exp_v_temp, model_v_temp)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
p_temp <- plot(exp_v_temp, fit_v_temp,
               main = "variogram temp")

p_temp

## direct variogram co-variable elev ----

exp_v_dem <- variogram(elev_m~1,
                       data = sp_temp_dem_s,
                       width = 0.1, 
                       cutoff = 10)
plot(exp_v_dem)

automap::autofitVariogram(elev_m~1, 
                          sp_temp_dem_s,
                          width = 0.1, 
                          cutoff = 10)$var_model

model <- "Ste"
psill <- 5*10^4
range <- 5
kappa <- 0.6

# model variogram
model_v_dem <- vgm(psill = psill, 
                   model = model, 
                   range = range, 
                   kappa = kappa)


# fit the model to the experimental variogram
# least square fit
fit_v_dem <- fit.variogram(exp_v_dem, model_v_dem)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
p_dem <- plot(exp_v_dem, fit_v_dem,
              main = "variogram elevation m")

p_dem

## compare two variograms
grid.arrange(p_temp, p_dem, ncol = 2)  # Multiplot 

# co-variogram or cross variogram for both variables ----
g_temp_dem <- gstat(NULL, 
                    id = "temp", 
                    formula = temp ~ 1, 
                    data = sp_temp_dem_s)

g_temp_dem <- gstat(g_temp_dem,
                    id = "dem",
                    formula = elev_m ~ 1,
                    data = sp_temp_dem_s)

# now display the two direct variograms and one cross-variogram
v_cross_temp_dem <- variogram(g_temp_dem, cutoff=10)

plot(v_cross_temp_dem, pl=F)

# add variogram models to the gstat object and
# fit a them using the linear model of co-regionalisation.
## By filling all the frames with one model (using the fill.all = T argument),
# these conditions are automatically met.

g_temp_dem <- gstat(g_temp_dem, 
                    id = "temp", 
                    model = model_v_temp, 
                    fill.all=T)

g_temp_dem

# Fit a Linear Model of Coregionalization 
# to a Multivariable Sample Variogram
g_temp_dem <- fit.lmc(v_cross_temp_dem, g_temp_dem)

g_temp_dem

plot(v_cross_temp_dem, model=g_temp_dem$model)

# coKriging prediction at grid locations -----
ck_temp_dem <- predict(g_temp_dem, sp_temp_dem_rv)

summary(ck_temp_dem)

plot(ck_temp_dem, main = "temp ck dem")
spplot(ck_temp_dem["temp.pred"], main = "temp ck dem")


####### UP TO HERE - OK - CODE CHECK! #################
#######################################################

# block kriging ----
