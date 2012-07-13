
## @knitr unnamed-chunk-1
library(ggplot2)
library(vcd)
gg.opts <- opts(panel.grid.major=theme_blank(), 
panel.grid.minor=theme_blank(), 
panel.background=theme_rect())
# maps
library(maps)
# make the results reproducable
set.seed(11123)


## @knitr unnamed-chunk-2
library(maptools)
survey.area<-readShapeSpatial("data/Study_ar")
survey.area<-data.frame(survey.area@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area)<-c("longitude","latitude")


## @knitr unnamed-chunk-3
load("data/dolphins.RData")

## @knitr unnamed-chunk-4
# centroid
lon0 <- -88.31951
lat0 <- 27.01594

# the dsm package has the function that we need
library(dsm)

sa.tmp <- latlong2km(survey.area$longitude, survey.area$latitude, 
lon0=lon0, lat0=lat0)

survey.area <- data.frame(x=1000*sa.tmp$km.e, y=1000*sa.tmp$km.n)

rm(sa.tmp)

### @knitr unnamed-chunk-7
#suppressPackageStartupMessages(library(Distance))
#
#
### @knitr unnamed-chunk-8
#hn.model<-ds(distdata,max(distdata$distance),monotonicity="strict")
#
### @knitr unnamed-chunk-9
#mod1<-dsm.fit(hn.model$ddf, response="indiv", formula=~s(x,y),
#obsdata=obsdata,segdata=segdata)
#
#
#off.set <- fitted(hn.model$ddf)[1]*preddata$width*preddata$height
#
#mod1.pred <- dsm.predict(mod1, preddata,off=off.set)
#pp <- cbind(preddata,mod1.pred)
#
### @knitr unnamed-chunk-15
#offset <- fitted(hn.model$ddf)[1]*pp$width*pp$height
#block.size<-4
#n.boot<-1000
#mod1.movblk <- dsm.var.movblk(mod1, preddata, n.boot,
#                              block.size=block.size,
#                              off.set=offset,bar=FALSE,
#                              bs.file="mexico-bs.csv",
#                              ds.uncertainty=TRUE)
#
#summary(mod1.movblk)

load("cvplot.RData")

## @knitr unnamed-chunk-17
plot(mod1.movblk,xlab="Metres from -88.31 longitude",ylab="Metres from 27.01 latitude",limits=c(0,10))
ggsave("cvplot-movblk.pdf",height=5,width=9)

quartz()

### @knitr unnamed-chunk-18
#preddata.varprop <- list()
#offset.varprop <- list()
#for(i in 1:nrow(preddata)){
#  preddata.varprop[[i]] <- preddata[i,]
#  offset.varprop[[i]] <- offset[i]
#}
#mod1.varprop<-dsm.var.prop(mod1,pred.data=preddata.varprop,off.set=offset.varprop)
#
#
### @knitr unnamed-chunk-19
#summary(mod1.varprop)


## @knitr unnamed-chunk-20
plot(mod1.varprop,xlab="Metres from -88.31 longitude",ylab="Metres from 27.01 latitude",limits=c(0,10))
ggsave("cvplot-varprop.pdf",height=5,width=9)


#save.image("cvplot.RData")
