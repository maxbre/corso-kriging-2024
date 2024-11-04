load('./ok_temp.rda')

raster::raster(ok_temp)
temp_ok<-terra::rast(ok_temp)
temp_ok

raster::raster(sp_temp_dem_rv)
rv<-terra::rast(sp_temp_dem_rv)
rv

samp <- terra::vect(sp_temp_dem_s)

delta <-rv["temp"] - temp_ok["var1.pred"]

plot(delta)
plot(samp, add=TRUE, col='red', pch=3, cex=0.3)

plot(temp_ok["var1.var"], 
     box=FALSE, 
     plg=list(x="top", at=c(0,0.5,1)),
     alpha=0.8)
plot(samp, add=TRUE, col='red', pch=3, cex=0.3)

par(mfrow=c(1,2))
plot(temp_ok["var1.pred"], box=FALSE, axes=FALSE, main = "estimated temp")
plot(temp_ok["var1.var"], box=FALSE, axes=FALSE, main = "estimated variance of temp") 

