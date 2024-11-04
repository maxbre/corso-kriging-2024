library(sf)

nc <- st_read(system.file("shape/nc.shp", package="sf"))
nc <- nc[1,]
set.seed(31234)
p1 <- st_sample(nc, 15)
p2 <- st_sample(nc, 1)
plot(st_geometry(nc))
plot(p1, add = TRUE, col="blue")
plot(p2, add = TRUE, pch = 2, col ="red")
