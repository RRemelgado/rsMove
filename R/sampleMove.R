#' @title sampleMove
#'
#' @description Sampling of possible stops along a movement track.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with the same length as \emph{xy}.
#' @param error Distance (in meters).
#' @param method How should the disntance be estimated? One of 'm' or 'deg'. Default is 'm'.
#' @param tUnit Time unit to estimate elapsed time. See \code{\link[base]{difftime}} for keywords. Default is \emph{mins}.
#' @import raster grDevices rgdal
#' @importFrom stats median
#' @return A \emph{SpatialPointsDataFrame}.
#' @details {This function offers a simple approach to sample from locati where an animal showed little or no movement 
#' based on GPS tracking data. It looks at the distance among consecutive samples (\emph{error}) and estimates mean coordinates 
#' for the temporal segments where the animal moved less than the defined distance from the first location of the segment. 
#' The user should selected \emph{method} in accordance with the projection system associated to the data. If 'm' it estimates 
#' the ecludian distance. If 'deg' it uses the haversine formula. The output reports on the mean sample coordinates for 
#' the sample locations ('x' and 'y'), the total time spent per sample ('time' expressed in minutes) and the total number 
#' of observations per sample ('count').}
#' @examples {
#'  
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'  
#'  # reference data
#'  file <- system.file('extdata', 'latLon_example.shp', package="rsMove")
#'  moveData <- shapefile(file)
#' 
#'  # sampling without reference grid
#'  ot = strptime(moveData$timestamp, "%Y-%m-%d %H:%M:%S")
#'  output <- sampleMove(xy=moveData, ot=ot, error=10, method='deg')
#'  
#'  # compare original vs new samples
#'  plot(moveData, col="black", pch=16)
#'  points(output$x, output$y, col="red", pch=15)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

sampleMove <- function(xy=xy, ot=ot, error=error, method='m', tUnit=NULL) {

#-----------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------------------------------------------------#

  # check input variables
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (!exists('ot')) {stop('"ot" is missing')}
  if (!class(ot)[1]%in%c('Date', 'POSIXlt', 'POSIXct')) {stop('"ot" is nof of a valid class')}
  if (length(ot)!=length(xy)) {stop('"xy" and "ot" have different lengths')}
  io <-order(ot) # index used to check data order
  xy <- xy[io,]
  ot <- ot[io]
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
    tt <- vector('character', ns)
    td <- 1:ns
    ss <- 1:ns
    for (r in 1:ns) {
      loc <- sc[[r]]
      xs[r] <- median(xy@coords[loc[1]:loc[2],1])
      ys[r] <- median(xy@coords[loc[1]:loc[2],2])
      tt[[r]] <- as.character(ot[loc[1]])
      td[r] <- as.numeric(difftime(ot[loc[2]], ot[loc[1]], units=tUnit))
      ss[r] <- length(loc[1]:loc[2])}
    
#-----------------------------------------------------------------------------------------------------------------------------#
# 4. build output
#-----------------------------------------------------------------------------------------------------------------------------#
    
    # if no layer is provided return the original sample set
    os <- data.frame(x=xs, y=ys, timestamp=tt, timeSum=td, count=ss, stringsAsFactors=F)
    os <- SpatialPointsDataFrame(cbind(xs,ys), os, proj4string=crs(xy))
    
    return(os)
    
  } else {stop('no samples found')}
}