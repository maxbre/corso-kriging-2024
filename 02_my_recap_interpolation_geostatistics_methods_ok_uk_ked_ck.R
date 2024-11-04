# load packages ----
library(tidyverse)
library(gstat)
library(sp)
library(automap)
library(gridExtra) # arrange plots
library(terra)
#library(raster)

# naming objects, for sure this is the most difficult part! -----
# general convention adopted here for naming objects
# template: "o_vars_geo"
# where:
# o = the object, e.g.: v for vectors, r for rasters, and so on...
# vars = the content of object, e.g.: temp, dem
# geo = the geographical boundaries, e.g.: rv for regione veneto, s for samplings


# import data_input ----

# read df temp and dem elevation for samplig points
sp_temp_dem_s <- read_csv('./data_input/df_temp_dem_rv_sampling.csv')

# create spatial objects for analysis with gstat
# sampling points
coordinates(sp_temp_dem_s) <- ~ x + y
proj4string(sp_temp_dem_s) <- CRS("EPSG:4326")
is.projected(sp_temp_dem_s)

# read temp and dem raster for RV
sp_temp_dem_rv <- read_csv('./data_input/r_temp_dem_rv.csv')

# turn it to a spatial grid
# for raster data, gstat only understands formats from the older sp package and the newer stars package.
coordinates(sp_temp_dem_rv) <- ~ x + y
gridded(sp_temp_dem_rv) <- TRUE
proj4string(sp_temp_dem_rv) <- CRS("EPSG:4326")
is.projected(sp_temp_dem_rv)

# to note it can be used also the void raster
# to be interpolated by the kriging methods
sp_shp_rv <- read_csv('./data_input/r_shp_rv.csv')

# the only precaution is that you need the same exact variables
# to be interpolated with

sp_shp_rv <- sp_shp_rv %>%
  rename(temp = shape_rv) %>%
  mutate(elev_m = -999)

# turn it to a spatial grid
# for raster data, gstat only understands formats from the older sp package and the newer stars package.
coordinates(sp_shp_rv) <- ~ x + y
gridded(sp_shp_rv) <- TRUE
proj4string(sp_shp_rv) <- CRS("EPSG:4326")
is.projected(sp_shp_rv)

# ordinary kriging -----

# analysis of variograms
# experimental sample variogram
# here using the same parametrization by pirotti exercise

# here two important parameters to define
#width: subsequent distance intervals into which data point pairs are grouped for semivariance estimates
#cutoff: spatial separation distance up to which point pairs are included in semivariance estimates; as a default, the length of the diagonal of the box spanning the data is divided by three.

plot(variogram(temp ~ 1, data = sp_temp_dem_s))

width <- 0.01 # pay attention: lon,lat because of original coordinates
cutoff <- 10

# define the experimental variogram
exp_v_ok <- variogram(temp ~ 1,
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

automap::autofitVariogram(temp ~ 1, sp_temp_dem_s, width = 0.01, cutoff = 10)$var_model

model <- "Ste"
psill <- 2
range <- 10
kappa <- 0.6

# model variogram
model_v_ok <- vgm(
  psill = psill,
  model = model,
  range = range,
  kappa = kappa
)


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

ok_temp <- krige(temp ~ 1, sp_temp_dem_s, #data
                 #sp_temp_dem_rv,   #newdata
                 sp_shp_rv, #newdata
                 model = fit_v_ok)

plot(ok_temp['var1.pred'], main = "temp ok temp~1")
plot(ok_temp['var1.var'], main = "var temp ok temp~1")

# transform to raster grid rv
r_temp_dem_rv <- terra::rast(sp_temp_dem_rv["temp"])
# transfor to vect sample points rv
v_temp_dem_s <- terra::vect(sp_temp_dem_s)

# transform to raster ordinary kriging for rv
r_temp_ok_rv <- terra::rast(ok_temp['var1.pred'])

# calculate delta
delta_r_temp_ok_rv <- r_temp_dem_rv - r_temp_ok_rv
# rename to something more meaningful
names(delta_r_temp_ok_rv) <- "delta_temp"

dev.off() # this is necessary to avoid an error of the plot device
plot(delta_r_temp_ok_rv)
plot(
  v_temp_dem_s,
  add = TRUE,
  col = 'red',
  pch = 3,
  cex = 0.3
)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_ok <- as.data.frame(delta_r_temp_ok_rv, xy = TRUE) %>% cbind(method =
                                                                     "krig_ord")

# universal kriging -----

width <- 0.01 # pay attention: lon,lat because of original coordinates
cutoff <- 10

# define the experimental variogram
exp_v_uk <- variogram(temp ~ x + y,
                      data = sp_temp_dem_s,
                      width = width,
                      cutoff = cutoff)
plot(exp_v_uk)
#plot(exp_v_uk, plot.numbers = TRUE)

# set model variogram
# intial parameter are set by eye esitmation
# model, sill, range, nugget

automap::autofitVariogram(temp ~ x + y, sp_temp_dem_s, width = 0.01, cutoff = 10)$var_model
model <- "Ste"
psill <- 6
range <- 10
kappa <- 0.6

# model variogram
model_v_uk <- vgm(
  psill = psill,
  model = model,
  range = range,
  kappa = kappa
)


# fit the model to the experimental variogram
# least square fit
fit_v_uk <- fit.variogram(exp_v_uk, model_v_uk)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
plot(exp_v_uk, fit_v_uk)
#plot(exp_v_uk, fit_v_uk, plot.numbers = TRUE)

uk_temp <- krige(temp ~ x + y, sp_temp_dem_s, #data
                 #sp_temp_dem_rv,    #newdata
                 sp_shp_rv, #newdata
                 model = fit_v_uk)

plot(uk_temp['var1.pred'], main = "temp uk temp~x+y")
plot(uk_temp['var1.var'], main = "var temp uk temp~x+y")

# transform to raster ordinary kriging for rv
r_temp_uk_rv <- terra::rast(uk_temp["var1.pred"])

# calculate delta
delta_r_temp_uk_rv <- r_temp_dem_rv - r_temp_uk_rv
# rename to something more meaningful
names(delta_r_temp_uk_rv) <- "delta_temp"

dev.off()
plot(delta_r_temp_uk_rv)
plot(
  v_temp_dem_s,
  add = TRUE,
  col = 'red',
  pch = 3,
  cex = 0.3
)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_uk <- as.data.frame(delta_r_temp_uk_rv, xy = TRUE) %>% cbind(method =
                                                                     "krig_univ")

# kriging with external drift -----

# now because we are considering the elevation
# as an external drift for the temperature
# we need to deal with a raster with both variable variables temp and elevation

width <- 0.01 # pay attention: lon,lat because of original coordinates
cutoff <- 10

# define the experimental variogram
exp_v_ked <- variogram(temp ~ elev_m + x + y,
                       data = sp_temp_dem_s,
                       width = width,
                       cutoff = cutoff)

plot(exp_v_ked)
#plot(exp_v_ked, plot.numbers = TRUE)

# set model variogram
# intial parameter are set by eye esitmation
# model, sill, range, nugget


automap::autofitVariogram(temp ~ elev_m + x + y,
                          sp_temp_dem_s,
                          width = 0.01,
                          cutoff = 10)$var_model

model <- "Ste"
psill <- 0.6
range <- 10
kappa <- 0.4

# model variogram
model_v_ked <- vgm(
  psill = psill,
  model = model,
  range = range,
  kappa = kappa
)


# fit the model to the experimental variogram
# least square fit
fit_v_ked <- fit.variogram(exp_v_ked, model_v_ked)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
plot(exp_v_ked, fit_v_ked)
#plot(exp_v_ked, fit_v_ked, plot.numbers = TRUE)

# kriging with external drift
# in fact it is an universal kriging
# with a variable other than coordinates x,y
ked_temp <- krige(temp ~ elev_m + x + y, sp_temp_dem_s, #data
                  sp_temp_dem_rv, #newdata
                  #sp_shp_rv,         #newdata
                  model = fit_v_ked)

plot(ked_temp['var1.pred'], main = "temp ked temp~elev_m+x+y")
plot(ked_temp['var1.var'], main = "var temp ked temp~elev_m+x+y")

# transform to raster ordinary kriging for rv
r_temp_ked_rv <- terra::rast(ked_temp["var1.pred"])

# calculate delta
delta_r_temp_ked_rv <- r_temp_dem_rv - r_temp_ked_rv
# rename to something more meaningful
names(delta_r_temp_ked_rv) <- "delta_temp"

dev.off() # this is necessary to avoid an error of the plot device
plot(delta_r_temp_ked_rv)
plot(
  v_temp_dem_s,
  add = TRUE,
  col = 'red',
  pch = 3,
  cex = 0.3
)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_ked <- as.data.frame(delta_r_temp_ked_rv, xy = TRUE) %>% cbind(method =
                                                                       "krig_extdr")


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

automap::autofitVariogram(temp ~ 1, sp_temp_dem_s, width = 0.1, cutoff = 10)$var_model

model <- "Ste"
psill <- 2.5
range <- 8
kappa <- 0.6

# model variogram
model_v_temp <- vgm(
  psill = psill,
  model = model,
  range = range,
  kappa = kappa
)

# fit the model to the experimental variogram
# least square fit
fit_v_temp <- fit.variogram(exp_v_temp, model_v_temp)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
p_temp <- plot(exp_v_temp, fit_v_temp, main = "variogram temp")

p_temp

## direct variogram co-variable elev ----

exp_v_dem <- variogram(elev_m ~ 1,
                       data = sp_temp_dem_s,
                       width = 0.1,
                       cutoff = 10)
plot(exp_v_dem)

automap::autofitVariogram(elev_m ~ 1, sp_temp_dem_s, width = 0.1, cutoff = 10)$var_model

model <- "Ste"
psill <- 5 * 10 ^ 4
range <- 5
kappa <- 0.6

# model variogram
model_v_dem <- vgm(
  psill = psill,
  model = model,
  range = range,
  kappa = kappa
)

# fit the model to the experimental variogram
# least square fit
fit_v_dem <- fit.variogram(exp_v_dem, model_v_dem)

# plot experimental variogram and fitted variogram
# against data (number of pairs)
p_dem <- plot(exp_v_dem, fit_v_dem, main = "variogram elevation m")

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
v_cross_temp_dem <- variogram(g_temp_dem, cutoff = 10)

plot(v_cross_temp_dem, pl = F)

# add variogram models to the gstat object and
# fit a them using the linear model of co-regionalisation.
## By filling all the frames with one model (using the fill.all = T argument),
# these conditions are automatically met.

g_temp_dem <- gstat(g_temp_dem,
                    id = "temp",
                    model = model_v_temp,
                    fill.all = T)

g_temp_dem

# Fit a Linear Model of Coregionalization
# to a Multivariable Sample Variogram
g_temp_dem <- fit.lmc(v_cross_temp_dem, g_temp_dem)

g_temp_dem

plot(v_cross_temp_dem, model = g_temp_dem$model)

# coKriging prediction at grid locations -----
ck_temp_dem <- predict(g_temp_dem, 
                       # sp_shp_rv     # newdata
                       sp_temp_dem_rv # newdata
                       )

summary(ck_temp_dem)

plot(ck_temp_dem['temp.pred'], main = "temp ck temp~elev_m+x+y")
plot(ck_temp_dem['temp.var'], main = "var temp ck temp~elev_m+x+y")

# transform to raster ordinary kriging for rv
r_temp_dem_ck_rv <- terra::rast(ck_temp_dem["temp.pred"])

# calculate delta
delta_r_temp_dem_ck_rv <- r_temp_dem_rv - r_temp_dem_ck_rv
# rename to something more meaningful
names(delta_r_temp_dem_ck_rv) <- "delta_temp"

dev.off() # this is necessary to avoid an error of the plot device
plot(delta_r_temp_dem_ck_rv)
plot(
  v_temp_dem_s,
  add = TRUE,
  col = 'red',
  pch = 3,
  cex = 0.3
)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_ck <- as.data.frame(delta_r_temp_dem_ck_rv, xy = TRUE) %>%
  cbind(method = "krig_co_dem")

# closing up ----

# bind rows all data frames ----
df_delta_temp <- bind_rows(df_dt_ok, df_dt_uk, df_dt_ked, df_dt_ck)

# export delta of temp for different kriiging methods
# to be later used in the analysis
write_csv(df_delta_temp,
          './data_input/df_delta_temp_methods_ok_uk_ked_ck.csv')

# plot of deltas for different interpolation approaches

# boxplots
df_delta_temp %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = delta_temp, colour = method)) +
  theme_minimal()

library(viridis)
library(ggthemes)

# maps of deltas
df_delta_temp %>%
  ggplot() +
  geom_raster(aes(x = x, y = y, fill = delta_temp)) +
  scale_fill_viridis(option = "plasma") +
  facet_wrap(vars(method)) +
  coord_sf(default_crs = sf::st_crs(4326)) +
  theme_map() +
  theme(legend.position = "right")


#save(....)
