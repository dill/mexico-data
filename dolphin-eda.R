
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
survey.area<-readShapeSpatial("data/Study_ar")
survey.area<-data.frame(survey.area@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area)<-c("longitude","latitude")
load("data/dolphins.RData")
# centroid
lon0 <- -88.31951
lat0 <- 27.01594
library(dsm)
sa.tmp <- latlong2km(survey.area$longitude, survey.area$latitude, 
lon0=lon0, lat0=lat0)
survey.area <- data.frame(x=1000*sa.tmp$km.e, y=1000*sa.tmp$km.n)
rm(sa.tmp)



### actual plotting code here


p <- ggplot(preddata)
p <- p + gg.opts
p <- p + coord_equal()
p <- p + labs(fill="Depth (m)",x="Metres from -88.31 longitude",y="Metres from 27.01 latitude",size="Group size")
p <- p + geom_tile(aes(x=x,y=y,fill=depth, width = width, height = height))
p <- p + geom_line(aes(x, y, group=Transect.Label),data=segdata)
p <- p + geom_point(aes(x, y, size=size), data=distdata, colour="red",alpha=I(0.7))
print(p)
ggsave("depth-transects.pdf",width=9,height=5)

## @knitr EDA-plots
# save graphics options
#o<-par("mfrow")

pdf("distances-groups.pdf",width=9,height=5)
par(mfrow=c(1,2))

# histograms
hist(distdata$distance,main="",xlab="Distance (m)")

# plots of distance vs. size
plot(distdata$distance,distdata$size, main="",xlab="Distance (m)",ylab="Group size",pch=19,cex=0.5,col=rgb(0.74,0.74,0.74,0.7))

# lm fit
l.dat<-data.frame(distance=seq(0,8000,len=1000))
lo<-lm(size~distance, data=distdata)
lines(l.dat$distance,as.vector(predict(lo,l.dat)))

#par(o)

dev.off()

### @knitr count-spatplot
#p<-qplot(data=survey.area,x=x,y=y,geom="polygon" , ylab="y", xlab="x", alpha=I(0.7),fill=I("lightblue"))
#p<-p+gg.opts
#p <- p + coord_equal()
#p <- p + labs(size="Group size")
#p <- p + geom_line(aes(x, y, group=Transect.Label),data=segdata)
#p <- p + geom_point(aes(x, y, size=size), data=distdata, colour="red",alpha=I(0.7))
#print(p)
#
#
### @knitr unnamed-chunk-5
#p <- ggplot(preddata)
#p <- p + gg.opts
#p <- p + coord_equal()
#p <- p + labs(fill="Depth",x="x",y="y")
#p <- p + geom_tile(aes(x=x,y=y,fill=depth, width = width, height = height))
#print(p)





