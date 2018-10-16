## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------
# load packages
require(rsMove)
require(raster)
require(ggplot2)
require(knitr)
require(kableExtra)
require(caret)
require(lattice)

## ----message=FALSE-------------------------------------------------------
data("shortMove") # movement data
ndvi <- stack(list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)) # environmental predictors

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE----
plot(calc(ndvi, max, na.rm=TRUE))
points(shortMove)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center"----
obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time), format="%Y/%m/%d %H:%M:%S") # format observation time
reduced.samples <- moveReduce(shortMove, ndvi, obs.time, derive.raster=TRUE) # remove redundant data points
plot(reduced.samples$total.time) # show total time image

## ----message=FALSE-------------------------------------------------------
upper.limit <- quantile(reduced.samples$total.time, 0.95) # identify upper threshold using 95%-percentile
move.mask <- reduced.samples$total.time > 0 & reduced.samples$total.time < upper.limit # build sample mask
usable.pixels <- which.max(move.mask) # identify relevant pixels
presence.samples <- SpatialPoints(xyFromCell(move.mask, usable.pixels), proj4string=crs(shortMove)) # build shapefile from samples (presences)

## ----message=FALSE-------------------------------------------------------
sample.id <- labelSample(presence.samples, ndvi, agg.radius=60) # aggregate samples in space
absence.samples <- backSample(presence.samples, ndvi, sample.id, sampling.method="pca") # identify absence samples
absence.samples # show samples

## ----message=FALSE, warning=FALSE, results='hide'------------------------
env.presences <- extract(ndvi, presence.samples) # extract environmental data for presences
env.absences <- extract(ndvi, absence.samples)  # extract environmental data for absences
resourceModel1 <- predictResources(env.presences, env.absences, sample.id, env.data=ndvi) # build model

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center"----
plot(resourceModel1$probabilities >= 0.5) # probability map
points(presence.samples) # presences

