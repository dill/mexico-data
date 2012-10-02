
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

library(maptools)
survey.area<-readShapeSpatial("../data/Study_ar")
survey.area<-data.frame(survey.area@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area)<-c("longitude","latitude")

load("../data/dolphins.RData")
# centroid
lon0 <- -88.31951
lat0 <- 27.01594

# the dsm package has the function that we need
library(dsm)

sa.tmp <- latlong2km(survey.area$longitude, survey.area$latitude, 
lon0=lon0, lat0=lat0)

survey.area <- data.frame(x=1000*sa.tmp$km.e, y=1000*sa.tmp$km.n)

rm(sa.tmp)

suppressPackageStartupMessages(library(Distance))
hn.model<-ds(distdata,max(distdata$distance),monotonicity="strict")

mod1<-dsm.fit(hn.model$ddf, response="indiv", formula=~s(x,y),
obsdata=obsdata,segdata=segdata)

off.set <- fitted(hn.model$ddf)[1]*preddata$width*preddata$height

mod1.pred <- dsm.predict(mod1, preddata,off=off.set)
pp <- cbind(preddata,mod1.pred)

## @knitr unnamed-chunk-15
offset <- fitted(hn.model$ddf)[1]*pp$width*pp$height
#block.size<-4
#n.boot<-1000
#mod1.movblk <- dsm.var.movblk(mod1, preddata, n.boot,
#                              block.size=block.size,
#                              off.set=offset,bar=FALSE,
#                              bs.file="cvplot-mexico-bs.csv",
#                              ds.uncertainty=TRUE)
#
#summary(mod1.movblk)



## @knitr unnamed-chunk-18
preddata.varprop <- list()
offset.varprop <- list()
for(i in 1:nrow(preddata)){
  preddata.varprop[[i]] <- preddata[i,]
  offset.varprop[[i]] <- offset[i]
}
mod1.varprop<-dsm.var.prop(mod1,pred.data=preddata.varprop,
                           off.set=offset.varprop)

#save.image("cvplot.RData")
#load("cvplot.RData")

summary(mod1.varprop)


# @knitr unnamed-chunk-20
#plot.limits <- c(0,100)
#plot.breaks <- c(seq(0,1,len=100),seq(2,10,1),seq(11,100,len=10),seq(100,1000,100))
#legend.breaks <- c(seq(0,1,len=5),10,1000)
plot.limits <- c(0,5)
plot.breaks <- seq(0,5,len=100)
legend.breaks <- round(seq(0,5,len=8),2)

#plot(mod1.movblk,xlab="Easting",ylab="Northing",
#     limits=plot.limits,breaks=plot.breaks,legend.breaks=legend.breaks)
##ggsave("cvplot-movblk.pdf",height=5,width=9)
#quartz()

# annoying switch to km
mod1.varprop$pred.data <- lapply(mod1.varprop$pred.data,
                                 function(x){x$x<-x$x/1000
                                             x$y<-x$y/1000
                                             x$width <- x$width/1000
                                             x$height <- x$height/1000
                                             return(x)})
mod1.varprop$dsm.object$data$x <- mod1.varprop$dsm.object$data$x/1000
mod1.varprop$dsm.object$data$y <- mod1.varprop$dsm.object$data$y/1000
mod1.varprop$dsm.object$ddf$data$x <- mod1.varprop$dsm.object$ddf$data$x/1000
mod1.varprop$dsm.object$ddf$data$y <- mod1.varprop$dsm.object$ddf$data$y/1000

plot(mod1.varprop,xlab="Easting",ylab="Northing",
     limits=plot.limits,breaks=plot.breaks,legend.breaks=legend.breaks)
ggsave("cvplot-varprop.pdf",height=5,width=9)

