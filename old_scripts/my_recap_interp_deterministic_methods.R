library(terra)
library(tidyverse)
library(viridis)

# read raster temperature
r_t <- rast("data/VenetoCorrectedMODIS_LST_Avg2017.tif")

# read RV gpkg as vector 
RVv <- vect("./data/veneto.gpkg")

# crop and mask RV raster temp
RVr_t <- crop(r_t, RVv, mask=TRUE)

plot(RVr_t)

# sample of 500 points within RV
RVv_t_s <- spatSample(RVr_t,
                        method="random",
                        size=500,
                        na.rm=T,
                        as.points=T)

# change the name to something more meaningful
names(RVv_t_s)<- "temp"

plot(RVr_t)
plot(RVv_t_s, add=TRUE, col="red", cex=0.3, pch=3)

############### interpolazione con metodi deterministici

# nearest neighbour (plain, no linear interpolation between points)
# interpolation between points is FALSE by defualt  
RVr_t_nn <- interpNear(RVr_t, RVv_t_s, field="temp", radius=1, interpolate=FALSE)

# mask
RVr_t_nn <- mask(RVr_t_nn, RVr_t)

plot(RVr_t_nn)

# delta
d_RVr_t_nn <- RVr_t_nn - RVr_t

# rename to something more meaningful
names(d_RVr_t_nn) <- "delta_t"

plot(d_RVr_t_nn)

hist(d_RVr_t_nn, breaks=100)
summary(d_RVr_t_nn)
quantile(d_RVr_t_nn$delta_t[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_nn <- as.data.frame(d_RVr_t_nn, xy=TRUE) %>% cbind(method="near_neigh")
#########

# nearest neighbour with linear interpolation between points
# interpolation between points is FALSE by defualt  

RVr_t_lnn <- interpNear(RVr_t, RVv_t_s, field="temp", radius=1, interpolate=TRUE)

# mask
RVr_t_lnn <- mask(RVr_t_lnn, RVr_t)

plot(RVr_t_lnn)

# delta
d_RVr_t_lnn <- RVr_t_lnn - RVr_t

# rename to something more meaningful
names(d_RVr_t_lnn) <- "delta_t"

plot(d_RVr_t_lnn)

hist(d_RVr_t_lnn, breaks=100)
summary(d_RVr_t_lnn)
quantile(d_RVr_t_lnn$delta_t[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_lnn <- as.data.frame(d_RVr_t_lnn, xy=TRUE) %>% cbind(method="ln_near_neigh")

##  IDW power 1 -------

RVr_t_p1 <- interpIDW(RVr_t, RVv_t_s, field="temp", radius=1, power=1)

# mask
RVr_t_p1 <- mask(RVr_t_p1, RVr_t)

plot(RVr_t_p1)

# delta
d_RVr_t_p1 <- RVr_t_p1 - RVr_t

# rename to something more meaningful
names(d_RVr_t_p1) <- "delta_t"

plot(d_RVr_t_p1)

hist(d_RVr_t_p1, breaks=100)
summary(d_RVr_t_p1)
quantile(d_RVr_t_p1$delta_t[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_p1 <- as.data.frame(d_RVr_t_p1, xy=TRUE) %>% cbind(method="idw_pw_1")


##  IDW power 2 -------

RVr_t_p2 <- interpIDW(RVr_t, RVv_t_s, field="temp", radius=1, power=2)

# mask
RVr_t_p2 <- mask(RVr_t_p2, RVr_t)

plot(RVr_t_p2)

# delta
d_RVr_t_p2 <- RVr_t_p2 - RVr_t

# rename to something more meaningful
names(d_RVr_t_p2) <- "delta_t"

plot(d_RVr_t_p2)

hist(d_RVr_t_p2, breaks=100)
summary(d_RVr_t_p2)
quantile(d_RVr_t_p2$delta_t[], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE)

# transform raster to a dataframe, add a new column for the method used for computation
df_dt_p2 <- as.data.frame(d_RVr_t_p2, xy=TRUE) %>% cbind(method="idw_pw_2")

##################

# bond rows all data frames
df_dt <- bind_rows(df_dt_nn, df_dt_lnn, df_dt_p1, df_dt_p2)

# plot delta

df_dt %>% 
  ggplot() + 
  geom_boxplot(aes(x=method, y=delta_t, colour = method))+
  theme_minimal()

df_dt %>% 
  ggplot()+
  geom_raster(aes(x=x, y=y, fill = delta_t))+
  scale_fill_viridis(option = "plasma")+
  facet_wrap(vars(method))+
  coord_sf(default_crs = sf::st_crs(4326))+
  theme_void()
