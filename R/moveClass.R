#' @title moveClass
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Translates movement data into regions where the animal remained 1 or more days.
#' @param xy Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonsDataFrame}.
#' @param obs.time Object of class \emph{POSIXlt} or \emph{POSIXct}.
#' @param min.trend Element of class \emph{numeric}.
#' @param min.distance Element of class \emph{numeric}.
#' @param nr.days Element of class \emph{numeric}.
#' @param distance.method One of "m" or "deg".
#' @return A \emph{SpatialPolygonsDataFrame}.
#' @importFrom raster rasterToPoints res crs cellStats area crop
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPolygonsDataFrame
#' @importFrom grDevices chull
#' @importFrom spatialEco polyPerimeter
#' @details {The output of the function is a \emph{SpatialPolygonsDataFrame} reporting on:
#' \itemize{
#'  \item{\emph{nrPoints} - Unique polygon identifier.}
#'  \item{\emph{nrDays} - Polygon Area (in square meters).}
#'  \item{\emph{start} - Polygon perimeter (in square meters).}
#'  \item{\emph{end} - Percentage of non-NA pixels.}}}
#' @seealso \code{\link{extractTS}} \code{\link{analyzeTS}}
#' @examples {}
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

moveClass <- function(xy, obs.time, min.trend=0.5, min.distance=10000, nr.days=3, distance.method="deg") {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  if (distance.method == 'deg') {
    distance.fun <- function(x, y) {
      rc <- rbind(x*pi/180, y*pi/180)
      xDiff <- abs(rc[2,1]-rc[1,1])
      yDiff <- abs(rc[2,2]-rc[1,2])
      aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(rc[2,2]) * cos(rc[2,1]) * sin(xDiff/2.) * sin(xDiff/2.)
      cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
      return(6371000 * cCoef)}}

  if (distance.method == 'm') {distance.fun <- function(x, y) {return(sqrt((x[1]-y[1])^2 + (x[2]-y[2])^2))}}

  od <- as.Date(obs.time) # observation days
  ud <- unique(od) # unique observation dates
  nd <- length(od) # number of records

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 2. extract daily statistics
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  odf <- do.call(rbind, sapply(1:nd, function(x) {

    ii <- which(obs.date==x) # target day
    nt <- ((hour(obs.time[ii])*3600) + (minute(obs.time[ii])*60) + second(obs.time[ii])) / 86400 # normalized time
    nr <- length(ii) # numer of records
    ld <- sapply(2:length(ii), function(i) {distance.fun(xy@coords[ii[i],1], xy@coords[(ii[i]-1),])}) # linear distance
    td <- sum(ld) # total distance
    dc <- cor(nt[2:nr], td) # time / distance correlation
    cx <- median(xy@coords[,1]) # center x
    cy <- median(xy@coords[,2]) # center y
    md <- max(runmed(ld, 3)) # max distance

    return(data.frame(x=cx, y=cy, n=nr, c=dc, m=md, t=td, d=x, stringsAsFactors=FALSE))

  }))

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 3. classify daily movements (distringuish between short and long-distance)
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  mc <- vector('numeric', nd)
  mc[1] <- 0

  for (d in 2:nd) {

    if (odf$m > min.distance & ld > min.trend) {

      ld <- sapply(1:nr.days, function(i) {distance.fun(odf[d,1:2], odf[(d+i),1:2])}) # check for consistency
      mc[d] <- ifelse(min(ld) > min.distance, 1, 0) # determine if a real long-distance movement

    }

  }

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 4. label segments of short and long-distance movements
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  # rfind segments
  si <- rle(mc)$lengths

  # identify sample indices
  ur <- do.call(rbind, lapply(1:length(si), function(i) {

    ii <- (1+sum(si[(i-1)])):sum(si[1:i])
    iv <- replicate(i, length(ii))

    return(data.frame(ii=ii, iv=iv))

  }))

  # assign labels
  sc <- vector('numeric', nd)
  sc[ur$ii] <- ur$iv

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 5. build bounding boxes for each segment and characterize them
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  # derive polygons
  ur <- unique(sc)

  pc <- lapply(1:length(ur), function(u) {

    si <- which(sc == u)
    ch <- chull(xy@coords[si,])
    if (length(ch) > 2) {
      p <- Polygons(list(Polygon(xy[si[c(ch,ch[1])],])), ID=u)
      return(list(rp=p, df=data.frame(np=length(si), nd=unique(length(od[si])), st=min(od[si]), et=max(od[si]), id=u, stringsAsFactors=FALSE)))
    } else {return(NULL)}

  })

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 6. build shapefile
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  pc <- pc[sapply(pc, function(i) {!is.null(i)})] # remove unused entries
  pc <- SpatialPolygonsDataFrame(lapply(pc, function(i) {i$rp}), do.call(rbind, lapply(pc, function(i) {i$df})), proj4string=crs(xy))
  colnames(pc$data) <- c('nrPoints', 'nrDays', 'start', 'end', 'ID')

}
