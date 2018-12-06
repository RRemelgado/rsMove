#' @title imgInt
#'
#' @description {Temporal linear interpolation of environmental data using
#' a \emph{raster}, \emph{SpatialPointsDataFrames} or \emph{matrix}/\emph{data.frame}.}
#' @param x Object of class \emph{RasterStack}, \emph{RasterBrick} or \emph{data.frame}.
#' @param y Object of class \emph{Date} with target dates. Alternatively, a \emph{RasterStack} or \emph{RasterBrick} with julian days for each pixel.
#' @param x.dates Object of class \emph{Date} with dates of \emph{x}.
#' @param time.buffer A two-element vector with temporal search buffer (expressed in days).
#' @param smooth Logical argument. Default is TRUE.
#' @param smooth.fun Smoothing function. uses \code{\link{runmean2}} by default.
#' @importFrom raster crs nlayers brick
#' @importFrom stats lm
#' @importFrom pryr mem_used
#' @importFrom utils memory.size
#' @importFrom lubridate is.Date
#' @useDynLib rsMove
#' @importFrom Rcpp sourceCpp evalCpp
#' @seealso \code{\link{dataQuery}} \code{\link{timeDir}} \code{\link{spaceDir}} \code{\link{moveSeg}}
#' @return A \emph{RasterBrick} or a \emph{data frame}. If a \emph{RasterBrick}, each layer represents a date in \emph{y}. If a \emph{data.frame}/\emph{matrix}, columns represent dates and rows represent samples.
#' @details {Wrapper for the function \code{\link{intime}} that performs a time-sensitive, linear interpolation
#' of a multi-band raster. The output dates are specified by \emph{y} and can differ from the dates of the input,
#' specified by \emph{x.dates}. \emph{time.buffer} controls the search for dates to interpolate from specifying
#' the maximum number of days that the selected data points can differ from the target date(s) in \emph{y}.
#' \emph{time.buffer} is provided as a two element vector which limits the search in the past and future. If
#' \emph{smooth} is TRUE, the function will also smooth the interpolated time series. \emph{fun} determines
#' which function to use. By default, \code{\link{runmean2}} is used whcih is an NA-sensitive, c++ implementation
#' of a simple running mean.}
#' @examples {
#'
#'  require(raster)
#'
#'  #'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)
#'  r.stk <- stack(file)
#'  r.stk <- stack(r.stk, r.stk, r.stk) # dummy files for the example
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # raster dates
#'  file.name <- names(r.stk)
#'  x.dates <- as.Date(paste0(substr(file.name, 2, 5), '-',
#'  substr(file.name, 7, 8), '-', substr(file.name, 10, 11)))
#'
#'  # interpolate raster data to target dates
#'  out <- imgInt(r.stk[1:50,1:50,drop=FALSE], x.dates, as.Date("2013-08-10"), c(60,60))
#'
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------#

imgInt <- function(x, x.dates, y, time.buffer, smooth=TRUE, smooth.fun=function(j) {runmean2(as.numeric(j), 1)}) {

#-----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------------#


  # check environmnetal data
  if (!class(x)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"x" is not of a valid class')}
  if (is.Date(x.dates)) {
    if (nlayers(x)!=length(x.dates)) {stop('"x" and "x.dates" have different lengths')}
    int.method <- 1}
  if (class(x.dates) %in% c('RasterStack', 'RasterBrick')) {
    if (nlayers(x)!=nlayers(x.dates)) {stop('"x" and "x.dates" have different lengths')}
    int.method <- 2}
  if (missing(int.method)) {stop('"x.dates" is not of a valid class')}

  # check temporal variables
  if(class(y)!='Date') {stop('"y" is not a "Date" object')}
  if (is.null(time.buffer)) {stop('"time.buffer" is missing')}
  if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}
  nl <- length(y) > 1 # will the output be a multi-band raster?

  # smoothing function
  if (!is.logical(smooth)) {stop('"smooth" is not a logical argument')}
  if (smooth) {if (!is.function(smooth.fun)) {stop('"smooth.fun" is not a function')}}

#-----------------------------------------------------------------------------------------------------------------------------------#
# 2. interpolate
#-----------------------------------------------------------------------------------------------------------------------------------#

  out <- brick(x[[1]], nl=length(y)) # output image stack

  # interpolation when dates are fixed
  if (int.method==1) {
    ov <- getValues(x) # import values into memory
    ov <- intime(ov, as.numeric(x.dates), as.numeric(y), time.buffer) # interpolate values
  }

  # interpolation when each pixel has its own date
  if (int.method==2) {
    v1 <- getValues(x) # import values into memory
    v2 <- getValues(x.dates)
    ov <- intime2(v1, v2, as.numeric(y), time.buffer) # interpolate values
    rm(v1, v2)
  }

  # smooth and convert to raster
  if (smooth & nl) {ov <- t(apply(out, 1:2, smooth.fun))}
  out <- setValues(out, ov) # assign data values
  names(out) <- as.character(y) # assign band names

  # provide output
  return(out)

}
