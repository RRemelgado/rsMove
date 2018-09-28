#' @title imgInt
#'
#' @description {Temporal linear interpolation of environmental data using
#' a \emph{raster}, \emph{SpatialPointsDataFrames} or \emph{matrix}/\emph{data.frame}.}
#' @param x Object of class \emph{RasterStack}, \emph{RasterBrick} or \emph{data.frame}.
#' @param y Object of class \emph{Date} with target dates.
#' @param x.dates Object of class \emph{Date} with dates of \emph{x}.
#' @param time.buffer A two-element vector with temporal search buffer (expressed in days).
#' @importFrom raster crs nlayers brick
#' @importFrom stats lm
#' @importFrom pryr mem_used
#' @importFrom utils memory.size
#' @seealso \code{\link{dataQuery}} \code{\link{timeDir}} \code{\link{spaceDir}} \code{\link{moveSeg}}
#' @return A \emph{RasterBrick} or a \emph{data frame}. If a \emph{RasterBrick}, each layer represents a date in \emph{y}. If a \emph{data.frame}/\emph{matrix}, columns represent dates and rows represent samples.
#' @details {Performs a pixel-wise linear interpolation over a raster for a given set of dates (\emph{y}).
#' A temporal buffer (\emph{time.buffer}) is required to limit the search for reference data points (\emph{time.buffer}).
#' This is defined by a two element vector which limits the search in the past and future.}
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
#'  out <- imgInt(extract(r.stk, shortMove), x.dates,
#'  as.Date("2013-08-10"), c(60,60))
#'
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------#

imgInt <- function(x, x.dates, y, time.buffer) {

#-----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------------#

  # check temporal variables
  if(class(y)!='Date') {stop('"y" is not a "Date" object')}
  if(class(x.dates)!='Date') {stop('"x.dates" is not a "Date" object')}
  if (is.null(time.buffer)) {stop('"time.buffer" is missing')}
  if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}

  # check environmnetal information
  if (!class(x)[1]%in%c('RasterStack', 'RasterBrick', 'data.frame', 'matrix')) {stop('"x" is not of a valid class')}
  if (class(x)[1]%in%c('RasterStack', 'RasterBrick')) {
    if (nlayers(x) != length(x.dates)) {stop('"x" and "x.dates" have a different dimensions')}
    processRaster <- TRUE}
  if (class(x)[1] %in% c("data.frame", "matrix")) {
    if (ncol(x)!=length(x.dates)) {stop('"x" and "x.dates" have different dimensions')}
    processRaster=FALSE}

#-----------------------------------------------------------------------------------------------------------------------------------#
# 2. build interpolation function
#-----------------------------------------------------------------------------------------------------------------------------------#

  intTime <- function(x) {

    tmp <- sapply(y, function(d) {

      di <- which(x.dates==d & !is.na(x))

      if (length(di) > 0) {return(mean(x[di]))} else {

        bi <- rev(which(!is.na(x) & x.dates < d & x.dates >= (d-time.buffer[1])))
        ai <- which(!is.na(x) & x.dates > d & x.dates <= (d+time.buffer[2]))

        if (length(bi)>=1 & length(ai)>=1) {
          lc <- lm(c(x[bi[1]],x[ai[1]])~as.numeric(c(x.dates[bi[1]],x.dates[ai[1]])))
          return(as.numeric(d)*lc$coefficients[2]+lc$coefficients[1])
        } else {return(NA)}

      }})

    return(tmp)

  }

#-----------------------------------------------------------------------------------------------------------------------------------#
# 3. interpolate
#-----------------------------------------------------------------------------------------------------------------------------------#

  # apply function (if raster)
  if (processRaster) {

      out <- brick(x[[1]], nl=length(y)) # output image stack
      v <- getValues(x) # import values into memory
      v <- t(apply(v, 1, intTime)) # interpolate values
      out <- setValues(out, v) # assign data values
      names(out) <- as.character(y) # assign band names

  }

  # apply function (if data frame/matrix)
  if (!processRaster) {

    out <- as.data.frame(apply(x, 1, intTime))
    if (length(y) > 1) {out <- t(out)}
    colnames(out) <- as.character(y)

  }

  # provide output
  return(out)

}
