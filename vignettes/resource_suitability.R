## ----echo=FALSE, message=FALSE-------------------------------------------
# load packages
require(rsMove)
require(raster)
require(ggplot2)
require(knitr)
require(kableExtra)

## ----message=FALSE-------------------------------------------------------
data("shortMove") # movement data
ndvi <- stack(list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)) # environmental predictors

## ----message=FALSE-------------------------------------------------------
plot(shortMove)

## ----message=FALSE-------------------------------------------------------
obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time), format="%Y/%m/%d %H:%M:%S") # format observation time
reduced.samples <- moveReduce(shortMove, ndvi, obs.time, derive.raster=TRUE) # remove redundant data points
plot(reduced.samples$total.time) # show total time image

## ----message=FALSE-------------------------------------------------------
move.mask <- reduced.samples$total.time > 0 & reduced.samples$total.time < 60 # build sample mask
usable.pixels <- which.max(move.mask) # identify relevant pixels
presence.samples <- SpatialPoints(xyFromCell(move.mask, usable.pixels), proj4string=crs(shortMove)) # build shapefile from samples (presences)

## ----error=TRUE----------------------------------------------------------
sample.id <- labelSample(presence.samples, ndvi, agg.radius=60) # aggregate samples in space
sample.id # show indices
absence.samples <- backSample(presence.samples, ndvi, sample.id, sampling.method="pca") # identify absence samples
absence.samples # show samples

## ----message=FALSE-------------------------------------------------------
env.presences <- extract(ndvi, presence.samples) # extract environmental data for presences
env.absences <- extract(ndvi, absence.samples)  # extract environmental data for absences
resourceModel1 <- predictResources(env.presences, env.absences, sample.id, env.data=ndvi) # build model

## ----message=FALSE-------------------------------------------------------
plot(resourceModel1$probabilities >= 0.5) # probability map
points(presence.samples) # presences

## ----message=FALSE-------------------------------------------------------
kable_styling(kable(head(resourceModel1$f1, 2), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

## ----message=FALSE-------------------------------------------------------
absence.samples <- backSample(presence.samples, ndvi, sampling.method="random") # identify absence samples (random)
env.absences <- extract(ndvi, absence.samples)  # extract environmental data for absences
resourceModel2 <- predictResources(env.presences, env.absences, sample.id, env.data=ndvi) # build model
plot(resourceModel2$probabilities >= 0.5) # probability map
points(presence.samples) # presences
kable_styling(kable(head(resourceModel2$f1, 2), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

## ----message=FALSE-------------------------------------------------------
landCover <- raster(system.file('extdata', 'landCover.tif', package="rsMove"))

## ----message=TRUE--------------------------------------------------------
class.labels <- c("Arable land", "Land without use", "Open spaces", "Wetlands", "Permanent crops", "Extraction/Dump sites", "Industrial areas", "Green urban areas")
probMask <- stack(resourceModel1$probabilities> 0.5, resourceModel2$probabilities> 0.5) # stack of probabilities (pca and random)
ptest <- plausibilityTest(probMask, landCover, class.labels=class.labels)

## ----message=TRUE--------------------------------------------------------
ptest$relative.plot
kable_styling(kable(head(ptest$relative.count, 8), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

