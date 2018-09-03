#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs using the nearest non NA value in time.
#' @param x Object of class \emph{RasterStack}, \emph{RasterBrick}.
#' @param y Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param x.dates Object of class \emph{Date} with \emph{x} observation dates.
#' @param y.dates Object of class \emph{Date} with \emph{y} observation dates.
#' @param time.buffer Two element vector with a temporal search buffer (expressed in days).
#' @param spatial.buffer Spatial buffer size used to smooth the returned values. The unit depends on the spatial projection.
#' @param smooth.fun Smoothing function applied with \emph{spatial.buffer}.
#' @importFrom raster crs extract nlayers
#' @importFrom stats median
#' @seealso \code{\link{sampleMove}} \code{\link{backSample}}
#' @return A \emph{data.frame} with the selected values and their corresponding dates.
#' @details {Returns environmental variables from a multi-layer raster object \emph{x} for a given set of coordinates (\emph{y}) depending
#' on the temporal distance between the observation dates (\emph{y.dates}) and the date on which each layer in the environmental data was
#' collected (\emph{x.dates}). \emph{time.buffer} controls the search for non-NA values in time and is adjusted to the observation date of
#' each element in \empH{y}. The user may also provide \emph{spatial.buffer} to spatially smooth the selected environmental information. In
#' this case, for each sample, the function will consider the neighboring pixels within the selected acquisition and apply a smoothing function
#' defined by \emph{smooth.fun}. If \emph{smooth.fun} is not specified, a weighted mean will be returned by default.}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
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
#'  # sample dates
#'  y.dates <- as.Date(shortMove@data$date)
#'
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(r.stk, shortMove, x.dates, y.dates, c(10,10))
#'
#' }
#'
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(x, y, x.dates, y.dates, time.buffer, spatial.buffer=NULL, smooth.fun=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  ev <- try(extract(x, y), silent=TRUE) # extract raster values (stop if failed)
  if (class(ev)[1] == "try-error") {stop('"x" and/or "y" are not valid inputs')}
  ev <- as.data.frame(ev)

  # check temporal information
  if (class(x.dates)[1]!='Date') {stop('"x.dates" is not of a valid class')}
  if (class(y.dates)[1]!='Date') {stop('"y.dates" is not of a valid class')}

  # auxiliary
  if (!is.null(spatial.buffer)) {if (!is.numeric(spatial.buffer)) {stop('"spatial.buffer" assigned but not numeric')}} else {smooth.fun=NULL}
  if (!is.null(spatial.buffer) & is.null(smooth.fun)) {smooth.fun <- function(x) {sum(x*x) / sum(x)}} else {
  if (!is.null(smooth.fun)) {if (!is.function(smooth.fun)) {stop('"smooth.fun" is not a valid function')}}}
  if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}
  if (length(time.buffer)!=2) {stop('"time.buffer" should be a two element vector')}

#-------------------------------------------------------------------------------------------------------------------------------#
# extract environmental data
#-------------------------------------------------------------------------------------------------------------------------------#

  # function to select pixels
  qf <- function(i) {
    ind <- which(!is.na(ev[i,]))
    if (length(ind)!=0) {
      diff <- abs(y.dates[i]-x.dates[ind])
      loc <- x.dates[ind] >= (y.dates[i]-time.buffer[1]) & x.dates[ind] <= (y.dates[i]+time.buffer[2])
      if (sum(loc)>0) {
        loc <- which(diff==min(diff[loc]))[1]
        return(list(value=ev[i,ind[loc]], date=x.dates[ind[loc]]))
      } else {return(list(value=NA, date=NA))}}
    else {return(list(value=NA, date=NA))}}

  # retrieve values
  ev <- lapply(1:nrow(ev), qf)
  orv <- do.call('c', lapply(ev, function(x) {x$value}))
  ord <- do.call('c', lapply(ev, function(x) {x$date}))

  # derive shapefile
  return(data.frame(value=orv, date=ord))

}
