# load packages ----
library(terra)
library(tidyverse)

# naming objects, for sure this is the most difficult part! -----
# general convention adopted here for naming objects
# template: "o_vars_geo"
# where:
# o = the object, e.g.: v for vectors, r for rasters, and so on...
# vars = the content of object, e.g.: temp, dem
# geo = the geographical boundaries, e.g.: rv for regione veneto, s for samplings


# import data_input ----

# read raster of temp for rv
r_temp_rv <- rast('./data_input/r_temp_rv.tif')

# read vector points for temp sampling
v_temp_rv_s <- vect('./data_input/v_temp_rv_sampling.gpkg')

# interpolation deterministic methods ----

# nn ------
# nearest neighbour (plain, no linear interpolation between points)
# interpolation between points is FALSE by defualt
r_temp_rv_nn <- interpNear(
  r_temp_rv,
  v_temp_rv_s,
  field = "temp",
  radius = 1,
  interpolate = FALSE
)

# mask
r_temp_rv_nn <- mask(r_temp_rv_nn, r_temp_rv)

# plot masked
plot(r_temp_rv_nn)

# calculate delta
delta_r_temp_rv_nn <- r_temp_rv_nn - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_nn) <- "delta_temp"

plot(delta_r_temp_rv_nn)

hist(delta_r_temp_rv_nn, breaks = 50)
summary(delta_r_temp_rv_nn)
quantile(delta_r_temp_rv_nn$delta_temp[],
         c(0.1, 0.25, 0.5, 0.75, 0.9),
         na.rm = TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_nn <- as.data.frame(delta_r_temp_rv_nn, xy = TRUE) %>% cbind(method =
                                                                     "near_neigh")

# nn_ iln ------
# nearest neighbour with linear interpolation between points
# interpolation between points is FALSE by defualt
# note here interpolate = TRUE
r_temp_rv_nn_iln <- interpNear(
  r_temp_rv,
  v_temp_rv_s,
  field = "temp",
  radius = 1,
  interpolate = TRUE
)

# mask
r_temp_rv_nn_iln <- mask(r_temp_rv_nn_iln, r_temp_rv)

# plot masked
plot(r_temp_rv_nn_iln)

# calculate delta
delta_r_temp_rv_nn_iln <- r_temp_rv_nn_iln - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_nn_iln) <- "delta_temp"

plot(delta_r_temp_rv_nn_iln)

hist(delta_r_temp_rv_nn_iln, breaks = 50)
summary(delta_r_temp_rv_nn_iln)
quantile(delta_r_temp_rv_nn_iln$delta_temp[],
         c(0.1, 0.25, 0.5, 0.75, 0.9),
         na.rm = TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_nn_iln <- as.data.frame(delta_r_temp_rv_nn_iln, xy = TRUE) %>% cbind(method =
                                                                             "near_neigh_iln")

# idw pw1 -----
# inverse distance weighting
# power 1

r_temp_rv_idwp1 <- interpIDW(
  r_temp_rv,
  v_temp_rv_s,
  field = "temp",
  radius = 1,
  power = 1
)

# mask
r_temp_rv_idwp1 <- mask(r_temp_rv_idwp1, r_temp_rv)

# plot mask
plot(r_temp_rv_idwp1)

# calculate delta
delta_r_temp_rv_idwp1 <- r_temp_rv_idwp1 - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_idwp1) <- "delta_temp"

plot(delta_r_temp_rv_idwp1)

hist(delta_r_temp_rv_idwp1, breaks = 50)
summary(delta_r_temp_rv_idwp1)
quantile(delta_r_temp_rv_idwp1$delta_temp[],
         c(0.1, 0.25, 0.5, 0.75, 0.9),
         na.rm = TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_idw_p1 <- as.data.frame(delta_r_temp_rv_idwp1, xy = TRUE) %>% cbind(method =
                                                                            "idw_pw1")

# idw pw2 -----
# inverse distance weighting
# power 2

r_temp_rv_idwp2 <- interpIDW(
  r_temp_rv,
  v_temp_rv_s,
  field = "temp",
  radius = 1,
  power = 2
)

# mask
r_temp_rv_idwp2 <- mask(r_temp_rv_idwp2, r_temp_rv)

# plot mask
plot(r_temp_rv_idwp2)

# calculate delta
delta_r_temp_rv_idwp2 <- r_temp_rv_idwp2 - r_temp_rv

# rename to something more meaningful
names(delta_r_temp_rv_idwp2) <- "delta_temp"

plot(delta_r_temp_rv_idwp2)

hist(delta_r_temp_rv_idwp2, breaks = 50)
summary(delta_r_temp_rv_idwp2)
quantile(delta_r_temp_rv_idwp2$delta_temp[],
         c(0.1, 0.25, 0.5, 0.75, 0.9),
         na.rm = TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_idw_p2 <- as.data.frame(delta_r_temp_rv_idwp2, xy = TRUE) %>% cbind(method =
                                                                            "idw_pw2")

# closing up ----

# bind rows all data frames ----
df_delta_temp <- bind_rows(df_dt_nn, df_dt_nn_iln, df_dt_idw_p1, df_dt_idw_p2)

# export delta of temp for different deterministic methods
# to be later used in the analysis
write_csv(df_delta_temp,
          './data_input/df_delta_temp_deterministiv_methods.csv')

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
