#' @title imgInt
#'
#' @description {Temporal linear interpolation of raster data for
#' \emph{rasters}, \emph{SpatialPointsDataFrames} or \emph{data frames}.}
#' @param img Object of class \emph{RasterStack} or \emph{RasterBrick}.
#' @param obs.dates Target dates. Object of class \emph{Date}.
#' @param r.dates Raster dates. Object of class \emph{Date}.
#' @param time.buffer wo element vector with temporal search buffer (expressed in days).
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param edata Object of class \emph{data frame} or \emph{matrix} with remote sensing data.
#' @import raster sp rgdal
#' @importFrom stats lm
#' @seealso @seealso \code{\link{dataQuery}} \code{\link{timeDir}} \code{\link{spaceDir}} \code{\link{moveSeg}}
#' @return A \emph{RasterBrick} or a \emph{data frame}.
#' @details {Performs a pixel-wise linear interpolation over a raster for a given set of dates (\emph{obs.dates}).
#' A teporal buffer (\emph{time.buffer}) is required to limit the search for reference data points (\emph{r.dates}). \emph{time.buffer}
#' is defined by a two element vector which limits the search of images in the past and future. If \emph{xy} is
#' provided the function only considers the pixels that overlap with the these sample points. Otherwise, a RasterBrick
#' is provided. However, if \emph{edata} is provided, \emph{xy} and \emph{img} are ignored. and a data frame is provided.}
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
#'  r.dates = seq.Date(as.Date("2012-01-01"), as.Date("2012-12-31"), 45)
#'
#'  # target dates
#'  obs.dates = as.Date("2012-04-01")
#'
#'  # interpolate raster data to target dates
#'  i.img <- imgInt(img=rsStk, r.dates=r.dates, obs.dates=obs.dates, time.buffer=c(60,60), xy=moveData)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------#

imgInt <- function(img=NULL, r.dates=r.dates, obs.dates=obs.dates, time.buffer=NULL, xy=NULL, edata=NULL) {

#-------------------------------------------------------------------------------------------------#
# 1. check input variables
#-------------------------------------------------------------------------------------------------#

  # check variables (if edata is provided)
  if(class(obs.dates)!='Date') {stop('"obs.dates" is not a "Date" object')}
  if(class(r.dates)!='Date') {stop('"r.dates" is not a "Date" object')}
  if (!is.null(edata)) {
    if (class(edata)[1]!='data.frame') {stop('"edata" is neither a matrix or a data frame')}
    if (ncol(edata)!=length(r.dates)) {stop('mismatch in "edata" and "r.dates" dimensions')}
    xy <- NULL
    img <- NULL
    time.buffer <- NULL
  } else {
    if (is.null(img)) {stop('"img" is missing')}
    if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
    if (!is.null(xy)) {
      if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {
        stop('"shp is nor a valid point shapefile object')}
      if (crs(xy)@projargs!=crs(img)@projargs) {
        stop('"xy" and "img" have different projections')}}
    if (nlayers(img)!=length(r.dates)) {stop('length of "img" and "r.dates" do not match')}
    if (!is.null(time.buffer)) {if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}}}

#-------------------------------------------------------------------------------------------------#
# 2. interpolate values
#-------------------------------------------------------------------------------------------------#

  # interpolation function
  otd <- ''
  int <- function(x) {
    di <- which(r.dates==otd & !is.na(x))
    if (length(di)>0) {return(mean(x[di]))} else {
      bi <- rev(which(!is.na(x) & r.dates < otd & r.dates >= (otd-time.buffer[1])))
      ai <- which(!is.na(x) & r.dates > otd & r.dates <= (otd+time.buffer[2]))
      if (length(bi)>=1 & length(ai)>=1) {
        lc <- lm(c(x[bi[1]],x[ai[1]])~as.numeric(c(r.dates[bi[1]],r.dates[ai[1]])))
        return(as.numeric(otd)*lc$coefficients[2]+lc$coefficients[1])
      } else {return(NA)}}}

  # apply function
  if (is.null(xy) & !is.null(img)) {
    out<-brick(img[[1]], nl=length(obs.dates))
    for (r in 1:length(obs.dates)) {
      otd <- obs.dates[r]
      out[[r]] <- calc(img, int)}}
  if (!is.null(xy)) {
    if (is.null(edata)) {edata <- extract(img, xy)}
    out <- matrix(0, length(xy), length(obs.dates))
    colnames(out) <- as.character(obs.dates)
    for (r in 1:length(obs.dates)) {
      otd <- obs.dates[r]
      out[,r] <- apply(edata, 1, int)}}

  # provide output
  return(out)

}
