###Application of kriging when the rainfall intensities are normally distributed ###

library(sp)
library(rgdal)
library(gstat)
library(rgeos)
library(zoo) 

setwd('..')
source("./2_Scripts/kriging_genutil.R")

num.pts=NULL

#inputs
npoint=8 #number of stations
avr.int=5 #averaging interval of rainfall intensity (min)
lowthres=10
upthres=200
vgbound=c(25,50,100,150,200,250,300,350,400)
#vgbound=c(75,105,125,175,250,305,400)
#vgbound=c(25,50,100,200,300,400)
#vgbound=c(1:5*100,1200,1500,2000,2500,3000,4000)


#set threshold
rain.data=read.csv(paste0("./1_Data/int_brad_",avr.int,".csv"))
coord=read.csv("./1_Data/loc_brad.csv")
geo_data=geostat_data(rain.data,coord,lowthres,upthres)
# read and explore data
d <- data.frame(geo_data)
names(d); dim(d)
max(d$x) - min(d$x); max(d$y) - min(d$y)
hist(d$int, col="Lightblue")

# small shift of x-coordinate to avoid predicting at observation locations
d$x <- d$x + 0.01

# 1. standardisation steps
mean.tp=rollapply(d$int,width=npoint,by=npoint, FUN=mean)
sd.tp=rollapply(d$int,width=npoint,by=npoint, FUN=sd)
plot(mean.tp, type="b"); plot(sd.tp)

d2 <- d
d2$mean.tp=rep(mean.tp,each=npoint)
d2$sd.tp=rep(sd.tp,each=npoint)
d2$std.int=(d2$int-d2$mean.tp)/d2$sd.tp
hist(d2$std.int, col="Lightblue")

# 2. variogram estimation

# large shift x-coordinate each time instant to enable variogram estimation using all data
# to achieve that only paired comparisons of observations at the same time are used
d2$x <- d2$x + 1000*(rep(1:(nrow(d2)/npoint),each=npoint)-1)

coordinates(d2) <- ~x+y
g <- gstat(formula=std.int~1, data=d2)
vg <- variogram(g, boundaries=vgbound)
vgm <- vgm(psill=1.2, "Exp", range=250, nugget=0.6)
vgm <- fit.variogram(vg, vgm)
plot(vg, vgm, plot.nu=T,main="Avr. Int : 15 min, Range : 5-20 mm/hr")

# 3. Kriging

# choose a time instant
i <- 1
d <- d[i:(i+7),]
coordinates(d) <- ~x+y

# scale by multiplying with the variance
vgm
vgm[,2] <- vgm[,2]*sd.tp[i]^2  
vgm

# 3.1 point kriging

# define grid
xy <- expand.grid(seq(415300,415735,1), seq(432715, 432890,1))
names(xy) <- c('x','y')
gridded(xy) <- ~x+y

pk <- krige(formula=int~1, locations=d, newdata=xy, model=vgm)
pk$var1.sd <- sqrt(pk$var1.var); spplot(pk, zcol="var1.pred"); spplot(pk, zcol="var1.sd")

# 3.2 point kriging

# define block
bdry <- as.data.frame(matrix(c(415300, 432715,415735,432715,415735,
  432890, 415300,432890), ncol=2, byrow=TRUE))
names(bdry) <- c("x","y")

# create SpatialPolygons from boundary
p <- Polygon(bdry)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps), proj4string=CRS(as.character(NA)))

# block kriging, must use predict.gstat instead of krige to control block discretization
g <- gstat(formula=int~1, locations=d, model=vgm)
bk <- predict(g, newdata=sps, sps.args=list(n=10000, type="regular"))
bk$var1.sd <- sqrt(bk$var1.var)

# 4. Comparison of results
# check that block prediction equals average of point predictions in area
mean(pk$var1.pred); bk$var1.pred
mean(pk$var1.pred) - bk$var1.pred
