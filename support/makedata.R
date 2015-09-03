library(dsm)
library(rgdal)
library(maptools)
library(plyr)

options(stringsAsFactors=FALSE)


obs<-read.csv("../data/csv/mexico-observation.csv",header=FALSE)
seg<-read.csv("../data/csv/mexico-segment.csv",header=FALSE)
#tran<-read.csv("../data/csv/mexico-transect.csv",header=FALSE)
pred <- readShapeSpatial("../data/DProject/Dolphins.dat/Predict")
pred.d <- read.table("../data/csv/mexico-pred.txt",skip=5)
survey.area <- readShapeSpatial("../data/DProject/Dolphins.dat/Study_Ar")

# proj 4 string
# using http://spatialreference.org/ref/esri/north-america-lambert-conformal-conic/
lcc_proj4 <- CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")

# centroid of the survey area
#lon0 <- -88.31951
#lat0 <- 27.01594


# create obsdata
# * ``object`` - object id
# * ``Segment.Label`` - the segment the observation occurred in
# * ``size`` - group size for the observation
# * ``distance`` - perpendicular/radial distance to observation

# get rid of the column titles put there by Distance
obs <- obs[6:nrow(obs),]



obsdata <- data.frame(
                      object=1:nrow(obs),
                      Sample.Label=obs[,13],
                      size=as.numeric(obs[,21]),
                      distance=as.numeric(obs[,19]),
                      Effort=1000*as.numeric(obs[,14])
                     )
obsdata$size[is.na(obsdata$size)] <- 0
obsdata$distance[is.na(obsdata$distance)] <- 0

# create segdata
# * ``x`` - centreline of the transect (i.e. "across the transect")
# * ``y`` - centre in the direction of the transect (i.e. "up/down")
# * ``Effort`` - the effort expended
# * ``Transect.Label`` - identifier for the transect this segment is in
# * ``Segment.Label`` - identifier for the segment (unique!)
# * ``esw`` - effective strip width... ?

# top 5 rows are the column titles from Distance
seg <- seg[6:nrow(seg),]

# make an sp object to project
segsp <- SpatialPoints(cbind(as.numeric(seg[,16]),
                             as.numeric(seg[,15])))

# give the sp object a projection
proj4string(segsp) <-CRS("+proj=longlat +datum=WGS84")
# re-project
segsp.t <- spTransform(segsp,CRSobj=lcc_proj4)

# build segment data
segdata <- data.frame(longitude      = as.numeric(seg[,16]),
                      latitude       = as.numeric(seg[,15]),
                      x              = segsp.t@coords[,1],
                      y              = segsp.t@coords[,2],
                      Effort         = 1000*as.numeric(seg[,14]),
                      Transect.Label = seg[,9],
                      Sample.Label   = seg[,13],
                      depth          = as.double(as.numeric(seg[,17]))
                     )

### distance data
# create what we want for the mrds analysis
distdata <- obsdata[obsdata$size>0, ]
distdata$Sample.Label <- NULL
distdata$detected <- rep(1,nrow(distdata))


## get other data
distdata$beaufort <- as.numeric(obs[,22][obsdata$size>0])
# include the lat/long but just for plotting
distdata$latitude <- as.numeric(obs[,15][obsdata$size>0])
distdata$longitude <- as.numeric(obs[,16][obsdata$size>0])
# make an sp object to project
distsp <- SpatialPoints(cbind(distdata$longitude,
                              distdata$latitude))

# give the sp object a projection
proj4string(distsp) <-CRS("+proj=longlat +datum=WGS84")
# re-project
distsp.t <- spTransform(distsp,CRSobj=lcc_proj4)
distdata$x <- distsp.t@coords[,1]
distdata$y <- distsp.t@coords[,2]


#### prediction data

proj4string(pred) <- CRS("+proj=longlat +datum=WGS84")

# grid spacing
grid.eps <- 0.16666667/2

# build some polygons
ID <- 0
pp <- alply(pred@coords, 1, function(x){

  this.poly <- rbind(c(x[1]+grid.eps,x[2]+grid.eps),
                     c(x[1]+grid.eps,x[2]-grid.eps),
                     c(x[1]-grid.eps,x[2]-grid.eps),
                     c(x[1]-grid.eps,x[2]+grid.eps),
                     c(x[1]+grid.eps,x[2]+grid.eps))

  ID <<- ID + 1
  Polygons(list(Polygon(this.poly)),ID=as.character(ID))
})

# make polygon object, give it a projection
pp <- SpatialPolygons(pp)
proj4string(pp) <- CRS("+proj=longlat +datum=WGS84")


# convert projection to lcc
# for the squares
pp.t <- spTransform(pp,CRSobj=lcc_proj4)
# for points
pred.t <- spTransform(pred,CRSobj=lcc_proj4)

# test plot
#par(mfrow=c(1,2))
#plot(pred, pch=19,cex=0.2)
#plot(pp,add=TRUE)
#plot(pred.t, pch=19,cex=0.2)
#plot(pp.t,add=TRUE)

pred.polys <- pp.t

preddata <- data.frame(latitude  = pred@coords[,2],
                       longitude = pred@coords[,1],
                       x         = pred.t@coords[,1],
                       y         = pred.t@coords[,2],
                       depth     = as.double(pred.d[,8]),
                       area      = ldply(pp.t@polygons[[1]]@Polygons,
                                   function(x) x@area)[[1]])


#preddata$width <- 1000*
#(lr.tmp$km.e[(length(preddata$latitude)+1):length(lr.tmp$km.e)]-
#lr.tmp$km.e[1:length(preddata$latitude)])
#preddata$height <- 1000*
#(tb.tmp$km.n[(length(preddata$latitude)+1):length(tb.tmp$km.n)]-
#tb.tmp$km.n[1:length(preddata$latitude)])

#rm(lr, tb, lr.tmp, tb.tmp, pred.tmp)

obsdata <- obsdata[obsdata$size>0,]

mexdolphins <- list(segdata     = segdata,
                    obsdata     = obsdata,
                    distdata    = distdata,
                    preddata    = preddata,
                    survey.area = survey.area,
                    pred.polys  = pred.polys)
# save everything to file
save(mexdolphins, file="../data/dolphins.RData")


# write the predictions as shapefiles
row.names(preddata) <- as.character(1:nrow(preddata))
pp.df <- SpatialPolygonsDataFrame(pp.t, preddata)
writeSpatialShape(pp.df, "prediction-grid")

proj4string(survey.area) <- CRS("+proj=longlat +datum=WGS84")
survey.area <- spTransform(survey.area,CRSobj=lcc_proj4)
writeSpatialShape(survey.area,"survey-area")

