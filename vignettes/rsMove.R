## ----echo=FALSE, message=FALSE-------------------------------------------
# load packages
library(rsMove)
library(raster)
library(sp)
library(ggplot2)
library(knitr)
library(kableExtra)

## ----message=FALSE-------------------------------------------------------
data("longMove")
data("shortMove")

## ----message=FALSE-------------------------------------------------------
# read remote sensing data
ndvi <- stack(list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE))
landCover <- raster(system.file('extdata', 'landCover.tif', package="rsMove"))

# extract ndvi raster dates
file.name <- names(ndvi)
ndvi.dates <- as.Date(paste0(substr(file.name, 2, 5), '-', substr(file.name, 7, 8), '-', substr(file.name, 10, 11)))

## ----message=FALSE-------------------------------------------------------
sample.regions <- hotMove(xy=longMove, pixel.res=0.1, return.shp=TRUE)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center"----
par(mar=c(4,4,0,4), xpd = NA, font.lab=2)
plot(longMove@coords[,1], longMove@coords[,2], pch=16, cex=0.5, xlab="Lon", ylab="Lat", cex.lab=1, cex.axis=1)
plot(sample.regions$polygons, col=rgb(1,0,0,0.3), add=TRUE)

## ------------------------------------------------------------------------
region.stats <- hotMoveStats(region.id=sample.regions$indices, obs.time=as.Date(longMove@data$timestamp))

## ---- echo=FALSE, results='asis'-----------------------------------------
kable_styling(kable(head(region.stats$region.stats, 5), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE----
region.stats$plot

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE, message=FALSE----
df <- data.frame(id=1:21, time=region.stats$region.stats$`Total Time`)
pol <- fortify(SpatialPolygonsDataFrame(sample.regions$polygons, df))
pol <- merge(pol, df, by="id")
ggplot(pol, aes(x=long, y=lat, group=id, fill=time)) + theme_bw() + geom_polygon() + xlab("Long") + ylab("Lat")

## ------------------------------------------------------------------------
s.res <- sMoveRes(xy=shortMove, pixel.res=c(10, 30, 250))

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE, message=FALSE----
kable_styling(kable(head(s.res$stats, 3), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")
s.res$plot

## ------------------------------------------------------------------------
t.res <- tMoveRes(xy=longMove, obs.date=as.Date(longMove@data$timestamp), time.res=c(1,8,16), pixel.res=0.01)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE, message=FALSE----
kable_styling(kable(head(t.res$stats, 3), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")
t.res$plot

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
plot(ndvi[[2]])

## ---- warning=FALSE, message=FALSE---------------------------------------
s.var <- specVar(img=ndvi[[2]], pixel.res=250)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
plot(s.var$mape)
s.var$plot

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE----
plot(ndvi[[1]], ext=shortMove)
points(shortMove, type="l")

## ------------------------------------------------------------------------
obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time), format="%Y/%m/%d %H:%M:%S")
reduced.samples <- moveReduce(xy=shortMove, obs.time=obs.time, img=landCover)

