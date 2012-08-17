
## @knitr unnamed-chunk-1
library(ggplot2)
library(vcd)
gg.opts <- opts(panel.grid.major=theme_blank(), 
panel.grid.minor=theme_blank(), 
panel.background=theme_rect())
# maps
library(maps)
# make the results reproducable
set.seed(1123)


## @knitr unnamed-chunk-2
library(maptools)
survey.area<-readShapeSpatial("data/Study_ar")
survey.area<-data.frame(survey.area@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area)<-c("longitude","latitude")


## @knitr unnamed-chunk-3
load("data/dolphins.RData")



## @knitr unnamed-chunk-4
lon0 <- -88.31951
lat0 <- 27.01594

# the dsm package has the function that we need
library(dsm)

sa.tmp <- latlong2km(survey.area$longitude, survey.area$latitude, 
lon0=lon0, lat0=lat0)

survey.area <- data.frame(x=1000*sa.tmp$km.e, y=1000*sa.tmp$km.n)

rm(sa.tmp)


## @knitr unnamed-chunk-5
segdata$Effort<-segdata$Effort*1000

seg.tmp <- latlong2km(segdata$longitude, segdata$latitude, lon0=lon0, lat0=lat0)

segdata <- cbind(segdata, x=1000*seg.tmp$km.e, y=1000*seg.tmp$km.n)
rm(seg.tmp)


## @knitr unnamed-chunk-6
ds.tmp <- latlong2km(distdata$longitude, distdata$latitude, lon0=lon0, lat0=lat0)

distdata <- cbind(distdata, x=1000*ds.tmp$km.e, y=1000*ds.tmp$km.n)
rm(ds.tmp)


## @knitr unnamed-chunk-7
obsdata$Effort<-obsdata$Effort*1000


## @knitr unnamed-chunk-8
pred.tmp <- latlong2km(preddata$longitude, preddata$latitude, lon0=lon0, lat0=lat0)

preddata <- cbind(preddata, x=1000*pred.tmp$km.e, y=1000*pred.tmp$km.n)
rm(pred.tmp)


## @knitr unnamed-chunk-9
suppressPackageStartupMessages(library(Distance))


## @knitr unnamed-chunk-10
hn.model<-ds(distdata,max(distdata$distance),monotonicity="strict")
summary(hn.model)


obsdata[,"size"][obsdata[,"size"]>0] <- 1

responsedata<-aggregate(obsdata[,"size"],
                            list(obsdata[,"Sample.Label"]), sum)
names(responsedata) <- c("seglab","count")


ss<-data.frame(x=segdata$x,
               y=segdata$y,
               n=responsedata$count,
               seglab=responsedata$seglab,
               height=segdata$Effort,
               width=fitted(hn.model$ddf)[1]*max(distdata$distance))


p<-ggplot(ss)+gg.opts
p <- p + labs(fill="Groups")
p<-p+geom_tile(aes(x=x,y=y,fill=n,width=width,height=height))
p <- p + coord_equal()
p <- p + scale_fill_gradientn(colours=heat_hcl(20),limits=c(0,4))
p<-p+geom_path(aes(x=x, y=y),data=survey.area)
print(p)

## @knitr hn-detfct
#plot(hn.model)


### @knitr unnamed-chunk-11
#mod1<-dsm.fit(hn.model$ddf, response="group", formula=~s(x,y),
#obsdata=obsdata,segdata=segdata)
#summary(mod1)
#
#dsm.check(mod1)

## @knitr unnamed-chunk-12
# need to faff around a bit here for plotting since we need the width and height of each cell
#lr <- c(preddata$longitude-1/6, preddata$longitude+1/6)
#tb <- c(preddata$latitude-1/6, preddata$latitude+1/6)
#
#lr.tmp <- latlong2km(lr, rep(preddata$latitude,2), lon0=lon0, lat0=lat0)
#tb.tmp <- latlong2km(rep(preddata$longitude,2), tb, lon0=lon0, lat0=lat0)
#
#preddata$width <- 1000*
#(lr.tmp$km.e[(length(preddata$latitude)+1):length(lr.tmp$km.e)]-
#lr.tmp$km.e[1:length(preddata$latitude)])
#preddata$height <- 1000*
#(tb.tmp$km.n[(length(preddata$latitude)+1):length(tb.tmp$km.n)]-
#tb.tmp$km.n[1:length(preddata$latitude)])
#
#rm(lr, tb, lr.tmp, tb.tmp)
#
#off.set <- fitted(hn.model$ddf)[1]*preddata$width*preddata$height
#
#
### @knitr unnamed-chunk-13
#mod1.pred <- dsm.predict(mod1, preddata,off=off.set)
#
#
### @knitr unnamed-chunk-14
#pp <- cbind(preddata,mod1.pred)
#
#
### @knitr unnamed-chunk-15
#sum(mod1.pred)
#
#
### @knitr unnamed-chunk-16
#quantile(pp$mod1.pred)


## @knitr unnamed-chunk-17
# median # segs per transect is 9
#for(trans in unique(segdata$Transect.Label)){
#cat(sum(segdata$Transect.Label==trans),"\n")
#}


### @knitr unnamed-chunk-18
#offset <- fitted(hn.model$ddf)[1]*pp$width*pp$height
#block.size<-4
#nbs<-20
#mod1.movblk <- dsm.var.movblk(mod1, preddata, nbs,
#                              off.set=offset,
#                              block.size=block.size,
#                              ds.uncertainty=FALSE)
#plot(mod1.movblk,limits=c(0,3))
#quartz()
#
#mod1.movblk.2 <- dsm.var.movblk(mod1, preddata, nbs,
#                              off.set=offset,
#                              block.size=block.size,
#                              ds.uncertainty=FALSE,
#                              bs.file="bp.test.csv")
#
#plot(mod1.movblk.2,limits=c(0,3))

### @knitr unnamed-chunk-19
#summary(mod1.movblk)
#
#
### @knitr unnamed-chunk-20
#mod1.movblk.dsu <- dsm.var.movblk(mod1, preddata, 10,
#                                  off.set=offset,
#                                  block.size=block.size,
#                                  ds.uncertainty=TRUE)
#
#
### @knitr unnamed-chunk-21
#summary(mod1.movblk.dsu)
#
#
### @knitr unnamed-chunk-22

#pd <- list()
#of <- list()
#for(i in 1:nrow(preddata)){
#  pd[[i]] <- preddata[i,]
#  of[[i]] <- offset[i]
#}
#preddata<-pd
#offset<-of
#
#
#mod1.varprop<-dsm.var.prop(mod1,pred.data=preddata,off.set=offset)
#
#
#plot(mod1.varprop)
#
#
### @knitr unnamed-chunk-23
#summary(mod1.varprop)
#
#
### @knitr unnamed-chunk-24
#mod2<-dsm.fit(hn.model$ddf, response="group", formula=~s(x,y)+s(depth), obsdata=obsdata,segdata=segdata)
#summary(mod2)
#
#
### @knitr mod2-preds
#mod2.pred<-dsm.predict(mod2, preddata,off=off.set)
#pp<-cbind(preddata,mod2.pred)
#p<-ggplot(pp)+gg.opts
#p <- p + labs(fill="Groups")
#p<-p+geom_tile(aes(x=x,y=y,fill=mod2.pred,width=width,height=height))
#p <- p + coord_equal()
#p <- p + scale_fill_gradientn(colours=heat_hcl(20),limits=c(0,4))
#p<-p+geom_path(aes(x=x, y=y),data=survey.area)
#print(p)
#
#
### @knitr unnamed-chunk-25
#quantile(pp$mod2.pred)
#
#
### @knitr unnamed-chunk-26
#mod3<-dsm.fit(hn.model$ddf, response="group", formula=~s(x,y), 
#		obsdata=obsdata,segdata=segdata,
#		model.defn=list(fn="gam",family="Tweedie",family.pars=list(p=1.2)))
#summary(mod3)
#
#
### @knitr unnamed-chunk-27
#soap.knots <- make_soap_grid(survey.area,c(11,6))
#
#
### @knitr unnamed-chunk-28
#onoff <- inSide(x=segdata$x, y=segdata$y, bnd=survey.area)
#segdata <- segdata[onoff,]
#
#
### @knitr unnamed-chunk-29
#mod4<-dsm.fit(hn.model$ddf, response="group", 
#formula=~s(x,y, bs="so",k=10,xt=list(bnd=list(survey.area)))+s(depth), 
#obsdata=obsdata, segdata=segdata, 
#model.defn=list(fn="gam",family="quasipoisson",knots=soap.knots))
#summary(mod4)
#
#
### @knitr mod4-pred
#mod4.pred<-dsm.predict(mod4,preddata, off=off.set)
#pp<-cbind(preddata,mod4.pred)
#
#p<-ggplot(pp)+gg.opts
#p<-p+geom_tile(aes(x=x,y=y,fill=mod4.pred,width=width,height=height))
#p <- p + coord_equal()
#p <- p + scale_fill_gradientn(colours=heat_hcl(20),limit=c(0,4))
#p<-p+geom_path(aes(x=x, y=y),data=survey.area)
#p <- p + labs(fill="Groups")
#print(p)
#
#
### @knitr unnamed-chunk-30
#hr.beau.model<-ds(distdata,max(distdata$distance),formula=~as.factor(beaufort),
#key="hr", adjustment=NULL)
#summary(hr.beau.model)
#hn.beau.model<-ds(distdata,max(distdata$distance),formula=~as.factor(beaufort),
#key="hn", adjustment=NULL)
#summary(hn.beau.model)
#
#
### @knitr unnamed-chunk-31
#mod5<-dsm.fit(hr.beau.model$ddf, response="group.est", formula=~s(x,y),
#obsdata=obsdata,segdata=segdata)
#summary(mod5)
#
#
### @knitr mod5-pred
#mod5.pred<-dsm.predict(mod5,preddata, off=off.set)
#pp<-cbind(preddata,mod5.pred)
#
#p <- ggplot(pp)+gg.opts
#p <- p + geom_tile(aes(x=x,y=y,fill=mod5.pred,width=width,height=height))
#p <- p + coord_equal()
#p <- p + scale_fill_gradientn(colours=heat_hcl(20),limit=c(0,4))
#p <- p + geom_path(aes(x=x, y=y),data=survey.area)
#p <- p + labs(fill="Groups")
#print(p)
#
#
### @knitr unnamed-chunk-32
#quantile(pp$mod5.pred)
#
#
