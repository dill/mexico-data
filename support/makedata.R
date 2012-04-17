### RUN THIS FROM THE TOP LEVEL!!!

# three csv files that are copy and pasted from Distance data viewer
# mexico-observation.csv
# mexico-segment.csv
# mexico-transect.csv -- pretty sure we can ignore this one

mex.obs<-read.csv("data/csv/mexico-observation.csv",header=FALSE)
mex.seg<-read.csv("data/csv/mexico-segment.csv",header=FALSE)
#mex.tran<-read.csv("data/csv/mexico-transect.csv",header=FALSE)

# create obsdata
# * ``object`` - object id
# * ``Segment.Label`` - the segment the observation occurred in
# * ``group.size`` - group size for the observation
# * ``distance`` - perpendicular/radial distance to observation

# get rid of the column titles put there by Distance
obs.tmp <- mex.obs[6:nrow(mex.obs),]
# take only those rows with observations in them (this is checking
# the "detected" column
obs.tmp <- obs.tmp[obs.tmp[,24]!="",]
obsdata <- data.frame(
                      object=1:nrow(obs.tmp),
                      Sample.Label=as.character(obs.tmp[,13]),
                      group.size=as.numeric(as.character(obs.tmp[,21])),
                      distance=as.numeric(as.character(obs.tmp[,19]))
                     )

# create segdata
# * ``x`` - centreline of the transect (i.e. "across the transect")
# * ``y`` - centre in the direction of the transect (i.e. "up/down")
# * ``Effort`` - the effort expended
# * ``Transect.Label`` - identifier for the transect this segment is in
# * ``Segment.Label`` - identifier for the segment (unique!)
# * ``esw`` - effective strip width... ?

# top 5 rows are the column titles from Distance
seg.tmp <- mex.seg[6:nrow(mex.seg),]

# use lat and long for now
segdata <- data.frame(
                      latitude=as.numeric(as.character(seg.tmp[,15])),
                      longitude=as.numeric(as.character(seg.tmp[,16])),
                      Effort=as.numeric(as.character(seg.tmp[,14])),
                      Transect.Label=as.character(seg.tmp[,9]),
                      Sample.Label=as.character(seg.tmp[,13]),
                      depth=as.numeric(as.character(seg.tmp[,17]))
                     )


# create what we want for the mrds analysis
distdata <- obsdata
distdata$Sample.Label <- NULL
distdata$detected <- rep(1,nrow(distdata))
distdata$beaufort <- as.numeric(as.character(obs.tmp[,22]))
# the group size must be called size for mrds
distdata$size <- distdata$group.size
distdata$group.size <- NULL

# include the lat/long but just for plotting
distdata$latitude <- as.numeric(as.character(obs.tmp[,15]))
distdata$longitude <- as.numeric(as.character(obs.tmp[,16]))

# save everything to file
save(segdata,obsdata,distdata,file="data/dolphins.RData")



