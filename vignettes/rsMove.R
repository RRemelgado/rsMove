## ----message=FALSE-------------------------------------------------------
# load package
library(rsMove)
library(raster)
library(sp)
require(ggplot2)

## ----message=FALSE-------------------------------------------------------
# read remote sensing data
files <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
r.stk <- stack(files)

# read animal movement data (migratory movements)
proj1 <- crs("+proj=longlat +ellps=WGS84 +no_defs")
move1 <- read.csv(system.file('extdata', 'latlon_example.csv', package="rsMove"))
move1 <- SpatialPointsDataFrame(move1[,2:3], move1, proj4string=proj1)

# read animal movement data (within-habitat movements)
proj2 <- crs(r.stk)
move2 <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
move2 <- SpatialPointsDataFrame(move2[,1:2], move2, proj4string=proj2)

## ----message=FALSE-------------------------------------------------------
# run function
sample.regions <- hotMove(xy=move1, pixel.res=0.1, return.shp=TRUE)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center"----
par(mar=c(4,4,0,4), xpd = NA, font.lab=2)
plot(move1@data$long, move1@data$lat, pch=16, cex=0.5, xlab="Lon", ylab="Lat", cex.lab=1, cex.axis=1)
plot(sample.regions$polygons, col=rgb(1,0,0,0.3), add=TRUE)

## ------------------------------------------------------------------------
region.stats <- hotMoveStats(region.id=sample.regions$indices, obs.time=as.Date(move1@data$timestamp))

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
s.res <- sMoveRes(xy=move2, pixel.res=c(10, 30, 250))

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE, message=FALSE----
knitr::kable(head(s.res$stats, 3))
s.res$plot

## ------------------------------------------------------------------------
t.res <- tMoveRes(xy=move1, obs.date=as.Date(move1@data$timestamp), time.res=c(1,8,16), pixel.res=0.01)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", echo=FALSE, message=FALSE----
knitr::kable(head(t.res$stats, 3))
t.res$plot

## ----message=FALSE-------------------------------------------------------
s.var <- specVar(img=r.stk[[1]], pixel.res=250)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
plot(s.var$mape)
s.var$plot

