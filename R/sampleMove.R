#' @title sampleMove
#'
#' @description Remote sensing oriented sampling of stops along a movement track.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with the same length as \emph{xy}.
#' @param error Distance (in meters).
#' @param method How should the disntance be estimated? One of 'm' or 'deg'. Default is 'm'.
#' @param tUnit Time unit to estimate elapsed time. See \code{\link[base]{difftime}} for keywords. Default is \emph{mins}.
#' @import raster rgdal
#' @importFrom stats median
#' @return A \emph{SpatialPointsDataFrame} with a reduced sample set.
#' @seealso \code{\link{labelSample}} \code{\link{backSample}} \code{\link{dataQuery}}
#' @details {This function finds location where an animal showed little or no movement based on GPS tracking data.
#' It looks at the distance among consecutive samples and places pointer when the distance is bellow the defined
#' threshold (\emph{error}). When a pointer is found, the function looks at the distance between the pointer and
#' the following samples. While this is below the distance threshold, the samples are assigned to the same segment.
#' Then, for each segment, the function summarizes the corresponding samples deriving mean coordinates, the start,
#' end and total time spent and the total number of samples per segment ('count'). The user should selected \emph{method}
#' in accordance with the projection system associated to the data. If 'm' it bases this analysis on the the ecludian distance.
#' However, if 'deg' it set, the function uses the haversine formula.}
#' @examples {
#'
#'  require(raster)
#'
#' # reference data
#' sprj <- crs("+proj=longlat +ellps=WGS84 +no_defs")
#' moveData <- read.csv(system.file('extdata', 'latlon_example.csv', package="rsMove"))
#' moveData <- SpatialPointsDataFrame(moveData[,2:3], moveData, proj4string=sprj)
#'
#'  # sampling without reference grid
#'  obs.time = strptime(moveData$timestamp, "%Y-%m-%d %H:%M:%S")
#'  output <- sampleMove(xy=moveData, obs.time=obs.time, error=7, method='deg')
#'
#'  # compare original vs new samples
#'  plot(moveData, col="black", pch=16)
#'  points(output$x, output$y, col="red", pch=15)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

sampleMove <- function(xy=xy, obs.time=obs.time, error=error, method='m', tUnit=NULL) {

#-----------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------#

  # check input variables
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (!exists('obs.time')) {stop('"obs.time" is missing')}
  if (!class(obs.time)[1]%in%c('Date', 'POSIXlt', 'POSIXct')) {stop('"obs.time" is nof of a valid class')}
  if (length(obs.time)!=length(xy)) {stop('"xy" and "obs.time" have different lengths')}
  io <-order(obs.time) # index used to check data order
  xy <- xy[io,]
  obs.time <- obs.time[io]
  rm(io)
  if (method!='m' & method!='deg') {stop(paste0('method ', method, ' not recognized'))}
  if (!is.null(tUnit)) {tUnit<-'days'}

#-----------------------------------------------------------------------------------------------------------------------------#
# 2. extract samples
#-----------------------------------------------------------------------------------------------------------------------------#

  # Identify time segments
  sc <- list()
  sp0 <- 0
  for (r in 2:length(xy)) {

    # Estimate distance (harvesine method)
    if (method=='deg') {
      if (sp0==0) {rc<-xy@coords[(r-1):r,]*pi/180} else {rc<-rbind(xy@coords[sp0,],xy@coords[r,])*pi/180}
      xDiff <- abs(rc[2,1]-rc[1,1])
      yDiff <- abs(rc[2,2]-rc[1,2])
      aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(rc[2,2]) * cos(rc[2,1]) * sin(xDiff/2.) * sin(xDiff/2.)
      cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
      lDist <- 6371000 * cCoef}

    # estimate distance (ecludian method)
    if (method=='m') {
      if (sp0==0) {rc <- xy@coords[(r-1):r,]} else {rc<-rbind(xy@coords[sp0,], xy@coords[r,])}
      lDist <- sqrt((rc[2,1]-rc[1,1])^2 + (rc[2,2]-rc[1,2])^2)}

    # determine if the sample belongs to a new segment
    if (lDist < error & sp0==0) {sp0 <- r-1}
    if (lDist > error & sp0>0) {
      sc[[length(sc)+1]] <- c(sp0,(r-1))
      sp0 <- 0}

  }

#-----------------------------------------------------------------------------------------------------------------------------#
# 3. summarize samples and derive statistics
#-----------------------------------------------------------------------------------------------------------------------------#

  # continue if segments were detected
  ns <- length(sc)
  if (ns > 0) {

    # summarize time segments
    xs <- 1:ns
    ys <- 1:ns
    st <- vector('list', ns)
    et <- vector('list', ns)
    td <- 1:ns
    ss <- 1:ns
    for (r in 1:ns) {
      loc <- sc[[r]]
      xs[r] <- median(xy@coords[loc[1]:loc[2],1])
      ys[r] <- median(xy@coords[loc[1]:loc[2],2])
      st[[r]] <- obs.time[loc[1]]
      et[[r]] <- obs.time[loc[2]]
      td[r] <- as.numeric(difftime(obs.time[loc[2]], obs.time[loc[1]], units=tUnit))
      ss[r] <- length(loc[1]:loc[2])}

    st <- do.call('c', st)
    et <- do.call('c', et)

#-----------------------------------------------------------------------------------------------------------------------------#
# 4. build output
#-----------------------------------------------------------------------------------------------------------------------------#

    # if no layer is provided return the original sample set
    os <- data.frame(x=xs, y=ys, start.time=st, end.time=et, timeSum=td, count=ss, stringsAsFactors=F)
    os <- SpatialPointsDataFrame(cbind(xs,ys), os, proj4string=crs(xy))

    return(os)

  } else {stop('no samples found')}
}
