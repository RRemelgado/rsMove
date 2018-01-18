## ----echo=FALSE, message=FALSE-------------------------------------------
# load packages
library(rsMove)
library(raster)
library(sp)
require(ggplot2)

## ----message=FALSE-------------------------------------------------------
data("longMove")
data("shortMove")

## ----message=FALSE-------------------------------------------------------
# read remote sensing data
ndvi <- stack(list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE))
landCover <- raster(system.file('extdata', 'landCover.tif', package="rsMove"))

## ----message=FALSE-------------------------------------------------------
# run function
sample.regions <- hotMove(xy=longMove, pixel.res=0.1, return.shp=TRUE)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center"----
par(mar=c(4,4,0,4), xpd = NA, font.lab=2)
plot(longMove@data$long, longMove@data$lat, pch=16, cex=0.5, xlab="Lon", ylab="Lat", cex.lab=1, cex.axis=1)
plot(sample.regions$polygons, col=rgb(1,0,0,0.3), add=TRUE)

## ------------------------------------------------------------------------
region.stats <- hotMoveStats(region.id=sample.regions$indices, obs.time=as.Date(longMove@data$timestamp))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(region.stats$region.stats, 5))

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
knitr::kable(head(s.res$stats, 3), align="c")
s.res$plot

## ------------------------------------------------------------------------
t.res <- tMoveRes(xy=longMove, obs.date=as.Date(longMove@data$timestamp), time.res=c(1,8,16), pixel.res=0.01)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE, message=FALSE----
knitr::kable(head(t.res$stats, 3), align="c")
t.res$plot

## ----message=FALSE-------------------------------------------------------
s.var <- specVar(img=ndvi[[1]], pixel.res=250)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
plot(s.var$mape)
s.var$plot

## ------------------------------------------------------------------------
obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time), format="%Y/%m/%d %H:%M:%S")
seg <- moveSeg(xy=shortMove, obs.time=obs.time, env.data=landCover, data.type="cat")

## ---- echo=FALSE---------------------------------------------------------
seg$stats$time <- format(seg$stats$time, digits=3)
knitr::kable(head(seg$stats, 5), align="c")

