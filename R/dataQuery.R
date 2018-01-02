#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs using the nearest non NA value in time.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.date Object of class \emph{Date} with \emph{xy} observation dates.
#' @param env.data Object of class \emph{RasterStack}, \emph{RasterBrick} or \emph{data.frame}.
#' @param env.date Object of class \emph{Date} with \emph{env.data} observation dates.
#' @param time.buffer Two element vector with temporal search buffer (expressed in days).
#' @param spatial.buffer Spatial buffer size used to smooth the returned values. The unit depends on the spatial projection.
#' @param smooth.fun Smoothing function applied with \emph{spatial.buffer}.
#' @importFrom raster crs extract nlayers
#' @importFrom stats median
#' @seealso \code{\link{sampleMove}} \code{\link{backSample}}
#' @return An object of class \emph{data.frame} with the selected values and their corresponding dates.
#' @details {Returns environmental variables from a raster object for a given set of x and y coordinates depending on the
#' temporal distance between the sample observation date (\emph{obs.date}) and the date on which the environmental data was
#' collected (\emph{env.date}). Within the buffer specified by \emph{time.buffer}, the function will search for the nearest
#' non \emph{NA} value with the shortest temporal distance. The user can adjust \emph{time.buffer} to control which pixels
#' are considred in this analysis. For example, \emph{time.buffer} can be set to c(30,0) prompting the function to ignore
#' environmental information acquired after the sample observation date and limit the search to -30 days. If \emph{time.buffer}
#' is set to null all acquisitions are considered. The user may also provide \emph{spatial.buffer} to spatially smooth the selected
#' environmental information. In this case, for each sample, the function will consider the neighboring pixels within the selected
#' acquisition and aplly a smoothing function defined by \emph{smooth.fun}. If \emph{smooth.fun} is not specified, a weighted mean
#' will be returned by default. If \emph{env.data} is a \emph{data.frame} \emph{spatial.buffer} and \emph{smooth.fun} are ignored and
#' \emph{env.date} should refer to each column.}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'  rsStk <- stack(rsStk, rsStk, rsStk) # dummy files for the example
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130805-20130811.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[1:10,1:2], moveData[1:10,], proj4string=crs(rsStk))
#'
#'  # raster dates
#'  env.date <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'
#'  # sample dates
#'  obs.date <- as.Date(moveData@data$date)
#'
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=moveData, obs.date=obs.date, env.data=rsStk, env.date=env.date, time.buffer=c(30,30))
#'
#' }
#'
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(xy=xy, obs.date=obs.date, env.data=env.data, env.date=env.date, time.buffer=NULL, spatial.buffer=NULL, fun=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}

  # raster
  if (!exists('env.data')) {stop('"env.data" is missing')}
  if (!class(env.data)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"env.data" is not of a valid class')}
  if (crs(xy)@projargs!=crs(env.data)@projargs) {stop('"xy" and "env.data" have different projections')}

  # check temporal information
  if (!exists('env.date')) {stop('"env.date" is missing')}
  if (class(env.date)[1]!='Date') {stop('"env.date" is not of a valid class')}
  if (length(env.date)!=nlayers(env.data)) {stop('lengths of "env.date" and "env.data" differ')}
  if (is.null(obs.date)) {stop('"obs.date" is missing')}
  if (class(obs.date)[1]!='Date') {stop('"obs.date" is not of a valid class')}
  if (length(obs.date)!=length(xy)) {stop('lengths of "obs.date" and "xy" differ')}

  # auxiliary
  if (!is.null(spatial.buffer)) {if (!is.numeric(spatial.buffer)) {stop('"spatial.buffer" assigned but not numeric')}} else {fun=NULL}
  if (!is.null(spatial.buffer) & is.null(smooth.fun)) {smooth.fun <- function(x) {sum(x*x) / sum(x)}} else {
  if (!is.null(smooth.fun)) {if (!is.function(smooth.fun)) {stop('"smooth.fun" is not a valid function')}}}
  if (!is.null(time.buffer)) {
    if (!is.numeric(time.buffer)) {stop('"time.buffer" is not numeric')}
    if (length(time.buffer)!=2) {stop('"time.buffer" should be a two element vector')}}

#-------------------------------------------------------------------------------------------------------------------------------#
# extract environmental data
#-------------------------------------------------------------------------------------------------------------------------------#

  # read data
  if (!is.data.frame(env.data) {env.data <- as.data.frame(extract(env.data, xy@coords, buffer=spatial.buffer, fun=smooth.fun, na.rm=TRUE))}

  # function to select pixels
  if (is.null(time.buffer)) {
    qf <- function(i) {
      ind <- which(!is.na(env.data[i,]))
      if (length(ind)!=0) {
        diff <- abs(obs.date[i]-env.date[ind])
        loc <- which(diff==min(diff))[1]
        return(list(value=env.data[i,ind[loc]], date=env.date[ind[loc]]))
      } else{return(list(value=NA, date=NA))}}
  } else {
    qf <- function(i) {
      ind <- which(!is.na(env.data[i,]))
      if (length(ind)!=0) {
        diff <- abs(obs.date[i]-env.date[ind])
        loc <- env.date[ind] >= (obs.date[i]-time.buffer[1]) & env.date[ind] <= (obs.date[i]+time.buffer[2])
        if (sum(loc)>0) {
          loc <- which(diff==min(diff[loc]))[1]
          return(list(value=env.data[i,ind[loc]], date=env.date[ind[loc]]))
        } else {return(list(value=NA, date=NA))}}
      else {return(list(value=NA, date=NA))}}}

  # retrieve values
  env.data <- lapply(1:length(xy), qf)
  orv <- do.call('c', lapply(env.data, function(x) {x$value}))
  ord <- do.call('c', lapply(env.data, function(x) {x$date}))

  # derive shapefile
  return(data.frame(value=orv, date=ord))

}
