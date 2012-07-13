
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


## @knitr unnamed-chunk-7
suppressPackageStartupMessages(library(Distance))


## @knitr unnamed-chunk-8
hn.model<-ds(distdata,max(distdata$distance),monotonicity="strict")

off.set <- fitted(hn.model$ddf)[1]*preddata$width*preddata$height




## @knitr unnamed-chunk-21
mod2<-dsm.fit(hn.model$ddf, response="indiv", formula=~s(depth), obsdata=obsdata,segdata=segdata)
mod2.pred<-dsm.predict(mod2, preddata,off=off.set)
cat("depth abundance=",sum(mod2.pred),"\n")

mod2.xy<-dsm.fit(hn.model$ddf, response="indiv", formula=~s(x,y)+s(depth), obsdata=obsdata,segdata=segdata)

mod2.xy.pred<-dsm.predict(mod2.xy, preddata,off=off.set)
cat("depth+x-y abundance=",sum(mod2.xy.pred),"\n")


plot.lims<-c(min(mod2.pred,mod2.xy.pred),max(mod2.pred,mod2.xy.pred))


## @knitr mod2-preds
pp<-cbind(preddata,mod2.pred)
p<-ggplot(pp)+gg.opts
p <- p + labs(fill="Abundance",x="Metres from -88.31 longitude",y="Metres from 27.01 latitude")
p<-p+geom_tile(aes(x=x,y=y,fill=mod2.pred,width=width,height=height))
p <- p + coord_equal()
p <- p + scale_fill_gradientn(colours=heat_hcl(20),limits=plot.lims)
p<-p+geom_path(aes(x=x, y=y),data=survey.area)
print(p)
ggsave("fit-depth.pdf", width=9,height=5)

quartz()
pdf("fit-depth-gam.pdf", width=5,height=5)
plot(mod2$result,select=1)
dev.off()

quartz()

## @knitr mod2-preds
pp<-cbind(preddata,mod2.xy.pred)
p<-ggplot(pp)+gg.opts
p <- p + labs(fill="Abundance",x="Metres from -88.31 longitude",y="Metres from 27.01 latitude")
p<-p+geom_tile(aes(x=x,y=y,fill=mod2.xy.pred,width=width,height=height))
p <- p + coord_equal()
p <- p + scale_fill_gradientn(colours=heat_hcl(20),limits=plot.lims)
p<-p+geom_path(aes(x=x, y=y),data=survey.area)
print(p)
ggsave("fit-depth-xy.pdf", width=9,height=5)


