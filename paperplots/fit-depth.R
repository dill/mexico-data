# generate figure 2 and 3


## @knitr unnamed-chunk-1
library(ggplot2)
library(vcd)
gg.opts <- opts(panel.grid.major=theme_blank(),
                panel.grid.minor=theme_blank(),
                panel.background=theme_blank(),
                legend.key=theme_blank())
# maps
library(maps)
# make the results reproducable
set.seed(11123)
library(gridExtra)

library(maptools)
survey.area<-readShapeSpatial("../data/Study_ar")
survey.area<-data.frame(survey.area@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area)<-c("longitude","latitude")
load("../data/dolphins.RData")
lon0 <- -88.31951
lat0 <- 27.01594
library(dsm)
sa.tmp <- latlong2km(survey.area$longitude, survey.area$latitude,
                     lon0=lon0, lat0=lat0)
#survey.area <- data.frame(x=1000*sa.tmp$km.e, y=1000*sa.tmp$km.n)
survey.area <- data.frame(x=sa.tmp$km.e, y=sa.tmp$km.n)
rm(sa.tmp)

suppressPackageStartupMessages(library(Distance))
hn.model<-ds(distdata,max(distdata$distance),monotonicity="strict")

#off.set <- preddata$width*preddata$height
off.set <- 444*1000^2 
mod2<-dsm(N~s(depth), hn.model, segdata, obsdata)
mod2.pred<-predict(mod2, preddata, off.set)
cat("depth abundance=",sum(mod2.pred),"\n")

mod2.xy<-dsm(N~s(x,y)+s(depth), hn.model, segdata, obsdata)
mod2.xy.pred<-predict(mod2.xy, preddata, off.set)
cat("depth+x-y abundance=",sum(mod2.xy.pred),"\n")

## plot the check
#pdf("dsm-check.pdf",width=7,height=7)
#gam.check(mod2.xy)
#dev.off()

# plot stuff

plot.breaks <- c(0,10,50,100,200,500,1000,2000,3000)#,5000,7000,7500)
plot.lims<-c(0,max(mod2.pred,mod2.xy.pred))
grad.obj <- scale_fill_gradient(high="black",low="white",limits=plot.lims,trans="sqrt")


## plot smooth of depth
#quartz()
pdf("fit-depth-gam.pdf", width=5,height=5)
plot(mod2,select=1,ylab="Smooth of depth")
dev.off()

# from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## @knitr mod2-preds
pp<-cbind(preddata,mod2.pred)
pp$x<-pp$x/1000
pp$y<-pp$y/1000
pp$width <- pp$width/1000
pp$height <- pp$height/1000

p.depth<-ggplot(pp)+gg.opts
p.depth <- p.depth + labs(fill="Abundance",x="Easting",y="Northing")
p.depth<-p.depth+geom_tile(aes(x=x,y=y,fill=mod2.pred,width=width,height=height))
p.depth <- p.depth + coord_equal()
p.depth <- p.depth + grad.obj
p.depth<-p.depth+geom_path(aes(x=x, y=y),data=survey.area)
#print(p.depth)
#ggsave("fit-depth.pdf", width=9,height=5)


#quartz()

## predictons for xy and depth
pp<-cbind(preddata,mod2.xy.pred)
# to get northings and eastings
pp$x<-pp$x/1000
pp$y<-pp$y/1000
pp$width <- pp$width/1000
pp$height <- pp$height/1000
p.depth.xy<-ggplot(pp)+gg.opts
p.depth.xy <- p.depth.xy + labs(fill="Abundance",x="Easting",y="Northing")
p.depth.xy<-p.depth.xy+geom_tile(aes(x=x,y=y,fill=mod2.xy.pred,width=width,height=height))
p.depth.xy <- p.depth.xy + coord_equal()
p.depth.xy <- p.depth.xy + grad.obj
p.depth.xy<-p.depth.xy+geom_path(aes(x=x, y=y),data=survey.area)
#print(p.depth.xy)
#ggsave("fit-depth-xy.pdf", width=9,height=5)

multiplot(p.depth,p.depth.xy)
ggsave("fit-depths.pdf",width=9,height=10)


## uncertainty plotting

plot.limits <- c(0, 1000)
plot.breaks <- c(seq(0, 1, len = 100), seq(2, 10, 1), seq(11, 100, len = 10),
    seq(100, 1000, 100))
legend.breaks <- c(seq(0, 1, len = 5), 10, 1000)

mod2.xy.varprop <- dsm.var.prop(mod2.xy, split(preddata,1:nrow(preddata)), off.set)

fixit <- function(pp){
  pp$x <- pp$x/1000
  pp$y <- pp$y/1000
  pp$width <- pp$width/1000
  pp$height <- pp$height/1000
  return(pp)
}
mod2.xy.varprop$pred.data <- lapply(mod2.xy.varprop$pred.data,fixit)

plot(mod2.xy.varprop, xlab = "Easting", ylab = "Northing", limits = plot.limits,
    breaks = plot.breaks, legend.breaks = legend.breaks,observations=FALSE,poly=survey.area)

ggsave("cvplot-varprop.pdf", width=9,height=5)


