#' @title moveReduce
#'
#' @description Pixel based summary of movement data that preserves periodic movements.
#' @param x Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param z Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with the observation time of \emph{x}.
#' @param y Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param preserve.revisits Logical. Should the function preserve revisit patterns?
#' @param derive.raster Should a \emph{RasterLayer} with the total time per pixel be provided?
#' @importFrom raster crs cellFromXY rasterize
#' @importFrom sp SpatialPointsDataFrame
#' @seealso \code{\link{sampleMove}} \code{\link{moveSeg}}
#' @return A \emph{list} object.
#' @details {Translates (\emph{x}) into pixel coordinates within a reference raster (\emph{y}). The
#' function identifies temporal segments corresponding to groups of consecutive observations within
#' the same pixel. In this process, revisits to recorded pixels are preserved. Once the segments are
#' identified, the function derives mean x and y coordinates for each of them and evaluates the time
#' spent within each pixel. The function reports on the start and end timestamps and the elapsed time.
#' If \emph{preserve.revisits} is FALSE, the function will then summarize the output on a pixel level
#' summing the time spent at each pixel. Additionally, if \emph{derive.raster} is TRUE, the function
#' will derive a \emph{RasterLayer} with the same configuration as \emph{y} depicting the the total
#' amount of time spent per pixel. The output of the function consists of:
#' \itemize{
#' \item{\emph{points} - \emph{SpatialPointsDataFrame} with the reduced sample set.}
#' \item{\emph{total.time} - \emph{RasterLayer} depicting the total time spent at each pixel.}}}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # observation time
#'  z <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(shortMove, r, z, derive.raster=TRUE)
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

moveReduce <- function(x, y, z, preserve.revisits=TRUE, derive.raster=FALSE) {

#----------------------------------------------------------------------------------------------------------#
# 1. check input variables
#----------------------------------------------------------------------------------------------------------#

  # samples
  if (!class(x)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"x" is not of a valid class')}

  # sample dates
  if (!class(z)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"z" is nof of a valid class')}
  if (length(z)!=length(x)) {stop('errorr: "x" and "z" have different lengths')}

  # environmental data
  if (!class(y)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    stop('"y" is not a valid raster object')}
  if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "edata" have different projections')}

#----------------------------------------------------------------------------------------------------------#
# 2. identify segments for each temporal window
#----------------------------------------------------------------------------------------------------------#

  # convert x to single pixels
  os <- order(z)
  x <- x[os,]
  z <- z[os]
  sp <- cellFromXY(y, x@coords)
  up <- unique(sp)

  rm(os)

  # search for segments and return sample indices
  pd <- rle(sp)$lengths
  sg <- vector('numeric', length(sp))
  for (p in 1:length(pd)) {sg[(sum(pd[0:(p-1)])+1):sum(pd[1:p])] <- p}

  rm(pd)

  # estimate
  odf <- do.call(rbind, lapply(1:max(sg), function(s) {
    ind <- which(sg==s)
    pp <- sp[ind[1]]
    mx <- median(x@coords[ind,1])
    my <- median(x@coords[ind,2])
    s.time <- z[ind[1]]
    e.time <- z[ind[length(ind)]]
    d.time <- as.numeric(difftime(e.time, s.time, units='mins'))
    return(data.frame(x=mx, y=my, pos=pp, start.time=s.time, end.time=e.time, elapsed.time=d.time, segment.id=s))}))

  # if preserve.revisits is FALSE, reduce to unique pixels
  if (!preserve.revisits) {
    odf <- do.call(rbind, lapply(unique(odf$pos), function(p) {
      i <- which(up==p)
      xy <- xyFromCell(y, p)
      return(data.frame(x=xy[1], y=xy[2], pos=p, start.time=min(odf$start.time[i]),
                        end.time=max(odf$end.time[i]), elapsed.time=sum(odf$elapsed.time[i])))}))

  }

  # build statistic shapefile
  r.shp <- SpatialPointsDataFrame(odf[,c("x","y")], odf, proj4string=crs(x))

#----------------------------------------------------------------------------------------------------------#
# 3. derive single raster
#----------------------------------------------------------------------------------------------------------#

  if (derive.raster) {

    # find cell positions of reduced sample set
    sp <- cellFromXY(y, odf[,c("x", "y")])
    up <- unique(sp)

    # estimate time sum per cell
    t.sum <- sapply(up, function(p) {sum(odf$elapsed.time[which(sp==p)], na.rm=TRUE)})

    # build raster
    t.sum.r <- rasterize(xyFromCell(y, up), y, t.sum)

  } else {t.sum.r <- NULL}

#----------------------------------------------------------------------------------------------------------#
# 4. build output
#----------------------------------------------------------------------------------------------------------#

  return(list(points=r.shp, total.time=t.sum.r))





}

