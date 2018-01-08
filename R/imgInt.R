#' @title imgInt
#'
#' @description {Temporal linear interpolation of environmental data using
#' a \emph{raster}, \emph{SpatialPointsDataFrames} or \emph{data frames}.}
#' @param env.data Object of class \emph{RasterStack}, \emph{RasterBrick} or \emph{data.frame}.
#' @param target.dates Object of class \emph{Date} with target dates.
#' @param env.dates Object of class \emph{Date} with dates of \emph{env.data}.
#' @param time.buffer A two-element vector with temporal search buffer (expressed in days).
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @importFrom raster crs nlayers brick
#' @importFrom stats lm
#' @seealso @seealso \code{\link{dataQuery}} \code{\link{timeDir}} \code{\link{spaceDir}} \code{\link{moveSeg}}
#' @return A \emph{RasterBrick} or a \emph{data frame}. If a \emph{RasterBrick}, each layer represents a date. If a \emph{data.frame}, columns represent dates and rows represent samples.
#' @details {Performs a pixel-wise linear interpolation over a raster for a given set of dates (\emph{target.dates}).
#' A temporal buffer (\emph{time.buffer}) is required to limit the search for reference data points (\emph{time.buffer}).
#' This is defined by a two element vector which limits the search in the past and future. If \emph{xy} is provided and
#' \emph{env.data} is a \emph{raster} object the function only considers the pixels that overlap with the shapefile.
#' Otherwise, all pixels are considered providing a \emph{RasterBrick}. However, if \emph{env.data} is a \emph{data.frame},
#' \emph{xy} is ignored.}
#' @examples {
#'
#'  require(raster)
#'
#'  #'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'  rsStk <- stack(rsStk, rsStk, rsStk) # dummy files for the example
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130805-20130811.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[1:10,1:2], moveData[1:10,], proj4string=crs(rsStk))
#'
#'  # raster dates
#'  env.dates = seq.Date(as.Date("2012-01-01"), as.Date("2012-12-31"), 45)
#'
#'  # target dates
#'  target.dates = as.Date("2012-04-01")
#'
#'  # interpolate raster data to target dates
#'  i.env.data <- imgInt(env.data=rsStk, env.dates=env.dates, target.dates=target.dates, time.buffer=c(60,60), xy=moveData)
#'
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------#

imgInt <- function(env.data=env.data, env.dates=env.dates, target.dates=target.dates, time.buffer=time.buffer, xy=NULL) {

#-----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------------#

  # check temporal variables
  if(class(target.dates)!='Date') {stop('"target.dates" is not a "Date" object')}
  if(class(env.dates)!='Date') {stop('"env.dates" is not a "Date" object')}
  if (is.null(time.buffer)) {stop('"time.buffer" is missing')}
  if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}

  # check environmnetal information
  if (!class(env.data)[1]%in%c('RasterStack', 'RasterBrick', 'data.frame')) {stop('"env.data" is not of a valid class')}
  if (!class(env.data)[1]%in%c('RasterStack', 'RasterBrick')) {
    if (!is.null(xy)) {
      if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"shp is nor a valid point shapefile object')}
      if (crs(xy)@projargs!=crs(env.data)@projargs) {stop('"xy" and "env.data" have different projections')}
      if (nlayers(env.data)!=length(env.dates)) {stop('length of "env.data" and "env.dates" do not match')}}
    processRaster=TRUE}
  if (is.data.frame(env.data)) {
    if (ncol(env.data)!=length(env.dates)) {stop('mismatch in "env.data" and "env.dates" dimensions')}
    processRaster=FALSE}

#-----------------------------------------------------------------------------------------------------------------------------------#
# 2. build interpolation function
#-----------------------------------------------------------------------------------------------------------------------------------#

  otd <- ''
  int <- function(x) {
    di <- which(env.dates==otd & !is.na(x))
    if (length(di)>0) {return(mean(x[di]))} else {
      bi <- rev(which(!is.na(x) & env.dates < otd & env.dates >= (otd-time.buffer[1])))
      ai <- which(!is.na(x) & env.dates > otd & env.dates <= (otd+time.buffer[2]))
      if (length(bi)>=1 & length(ai)>=1) {
        lc <- lm(c(x[bi[1]],x[ai[1]])~as.numeric(c(env.dates[bi[1]],env.dates[ai[1]])))
        return(as.numeric(otd)*lc$coefficients[2]+lc$coefficients[1])
      } else {return(NA)}}}

#-----------------------------------------------------------------------------------------------------------------------------------#
# 3. read environmental data (if required)
#-----------------------------------------------------------------------------------------------------------------------------------#

  if (!is.null(xy) & processRaster) {
    env.data <- extract(env.data, xy)
    processRaster=FALSE}

#-----------------------------------------------------------------------------------------------------------------------------------#
# 4. interpolate
#-----------------------------------------------------------------------------------------------------------------------------------#

  # apply function (if raster)
  if (processRaster) {
    out <- brick(env.data[[1]], nl=length(target.dates))
    for (r in 1:length(target.dates)) {
      otd <- target.dates[r]
      out[[r]] <- calc(env.data, int)}}

  # apply function (if data frame)
  if (!processRaster) {
    out <- matrix(0, nrow(env.data), length(target.dates))
    colnames(out) <- as.character(target.dates)
    for (r in 1:length(target.dates)) {
      otd <- target.dates[r]
      out[,r] <- apply(env.data, 1, int)}}

  # provide output
  return(out)

}
