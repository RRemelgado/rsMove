#' @title moveReduce
#'
#' @description Pixel based summary of movement data that preserves periodic movements.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param derive.raster Should a raster with the total time per pixel be provided?
#' @importFrom raster crs cellFromXY rasterize
#' @importFrom sp SpatialPointsDataFrame
#' @seealso \code{\link{sampleMove}} \code{\link{moveSeg}}
#' @return A \emph{list} object.
#' @details {Reduces a set of input samples (\emph{xy}) based on their corresponding pixel coordinates
#' within a reference raster (\emph{img}). Using this data, the function identifies temporal segments
#' corresponding to groups of consecutive samples found within the same pixel. In this process, revisits
#' to recorded pixels are preserved. Once the segments are identified, the function derives mean x and y
#' coordinates for each of them and evaluates the time spent within each pixel. The function reports on
#' the start and end timestamps, the mean timestamp and the elapsed time. The output of the function
#' consists of:
#' \itemize{
#' \item{\emph{r.shp} - Shapefile with reduced sample set and its corresponding temporal information.}
#' \item{\emph{total.time} - Raster showing the total time spent at each pixel (if \emph{derive.raster} is TRUE).}
#' \item{\emph{indices} - Indices for each sample in \emph{xy} showing which samples were aggregated.}}
#'
#' }
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
#'  obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(xy=shortMove, obs.time=obs.time, img=r)
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

moveReduce <- function(xy=xy, obs.time=obs.time, img=img, derive.raster=FALSE) {

#----------------------------------------------------------------------------------------------------------#
# 1. check input variables
#----------------------------------------------------------------------------------------------------------#

  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy) # output projection

  # sample dates
  if (!is.null(obs.time)) {
    if (!class(obs.time)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"obs.time" is nof of a valid class')}
    if (length(obs.time)!=length(xy)) {stop('errorr: "xy" and "obs.time" have different lengths')}}

  # environmental data
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    stop('"img" is not a valid raster object')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "edata" have different projections')}

#----------------------------------------------------------------------------------------------------------#
# 2. identify segments
#----------------------------------------------------------------------------------------------------------#

  # convert xy to single pixels
  os <- order(obs.time)
  xy <- xy[os,]
  obs.time <- obs.time[os]
  sp <- cellFromXY(img, xy@coords)

  rm(os)

  # search for segments and return sample indices
  pd <- rle(sp)$lengths
  sg <- vector('numeric', length(sp))
  for (p in 1:length(pd)) {sg[(sum(pd[0:(p-1)])+1):sum(pd[1:p])] <- p}

  rm(pd)

  # estimate
  df <- do.call(rbind, lapply(1:max(sg), function(s) {
    ind <- which(sg==s)
    mx <- median(xy@coords[ind,1])
    my <- median(xy@coords[ind,2])
    s.time <- obs.time[ind[1]]
    e.time <- obs.time[ind[length(ind)]]
    d.time <- as.numeric(difftime(e.time, s.time, units='mins'))
    return(data.frame(x=mx, y=my, start.time=s.time, end.time=e.time, diff.time=d.time, segment.id=s))}))
  colnames(df) <- c("x", "y", "Timeststamp (start)", "Timeststamp (end)",
                    "Elapsed time (minutes)", "Segment ID")

  # build statistic shapefile
  r.shp <- SpatialPointsDataFrame(df[,1:2], df, proj4string=crs(xy))

#----------------------------------------------------------------------------------------------------------#
# 3. derive single raster
#----------------------------------------------------------------------------------------------------------#

  if (derive.raster) {

    # estimate time sum per cell
    up <- unique(sp)
    t.sum <- sapply(up, function(p) {sum(df$'Elapsed time (minutes)'[which(sp==p)], na.rm=TRUE)})

    # build raster
    t.sum.r <- rasterize(xyFromCell(img, up), crop(img, extent(xy)), t.sum)

  } else {t.sum.r <- NULL}

#----------------------------------------------------------------------------------------------------------#
# 4. build output
#----------------------------------------------------------------------------------------------------------#

  return(list(points=r.shp, total.time=t.sum.r, indices=sg))

}

