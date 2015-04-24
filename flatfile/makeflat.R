# make a "flatfile" version of the data for use with Distance2



load("../data/dolphins.RData")


## take the segment data as the base data

# select the labels and effort only
segdata <- unique(segdata[,c("Effort","Transect.Label")])

# replace Transect.Label with Sample.Label
segdata$Sample.Label <- segdata$Transect.Label
segdata$Transect.Label <- NULL

# sum the Effort
segdata <- t(simplify2array(by(segdata, segdata$Sample.Label,
                               function(x) c(x[1,2], sum(x[,1])))))
segdata <- data.frame(Sample.Label = row.names(segdata),
                      Effort       = segdata[,2])


## collect the correct distance columns
distdata <- unique(distdata[,c("object", "size", "detected",
                               "beaufort", "distance")])


# convert Sample.Label to sample rather than segment
obsdata$Sample.Label <- sub("-\\d$", "", obsdata$Sample.Label)
obsdata <- unique(obsdata[,c("object", "Sample.Label")])

## merge the distance data onto the observations
obsdata <- merge(obsdata, distdata, by="object")

# merge the observation data onto the segment data
mexdolphins <- merge(segdata, obsdata, by="Sample.Label", all.x=TRUE)

mexdolphins$Area <- sum(apply(preddata[,c("width","height")], 1, prod))

save(mexdolphins, file="mexdolphins.RData")
