
## @knitr unnamed-chunk-1
library(ggplot2)
library(vcd)
gg.opts <- opts(panel.grid.major=theme_blank(), 
panel.grid.minor=theme_blank(), 
panel.background=theme_rect(),
legend.key=theme_blank())
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
p <- p + labs(fill="Depth (m)",x="Metres from centre point",y="Metres from centre point",size="Group size")
p <- p + geom_tile(aes(x=x,y=y,fill=depth, width = width, height = height))
p <- p + geom_line(aes(x, y, group=Transect.Label),data=segdata)
p <- p + geom_point(aes(x, y, size=size), data=distdata, colour="blue",alpha=I(0.7))
p <- p + scale_fill_gradientn(colours = heat_hcl(10), 
                              limits=c(min(preddata$depth),max(preddata$depth)),
                              breaks=c(0,200,300,400,500,750,1000,2000,3000,3500))
print(p)
ggsave("depth-transects.pdf",width=9,height=5)

## @knitr EDA-plots
# save graphics options
#o<-par("mfrow")

pdf("distances-groups.pdf",width=9,height=5)
par(mfrow=c(1,2))

# histograms
#hist(distdata$distance,main="",xlab="Distance (m)")
library(Distance)

## @knitr unnamed-chunk-8
hn.model<-ds(distdata,max(distdata$distance),monotonicity="strict")
plot(hn.model,pl.den=0,showpoints=FALSE,main="")

# plots of distance vs. size
plot(distdata$distance,distdata$size, main="",xlab="Distance (m)",ylab="Group size",pch=19,cex=0.5,col=rgb(0.74,0.74,0.74,0.7))

# lm fit
l.dat<-data.frame(distance=seq(0,8000,len=1000))
lo<-lm(size~distance, data=distdata)
lines(l.dat$distance,as.vector(predict(lo,l.dat)))

dev.off()

