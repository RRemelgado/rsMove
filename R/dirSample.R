#' @title dirSample
#'
#' @description Forward and backward sampling of a raster time series using GPS tracking data.
#' @param xy Object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param ot "Date", POSIXlt" or "POSIXct" object with the same length as "xy".
#' @param img Object of class "RasterStack" or "RasterBrick".
#' @param rt "Date", POSIXlt" or "POSIXct" object with the same number of bands as "img".
#' @param mws Moving window size. Numeric.
#' @param mwu Moving window unit. One of "day", "month", "year". Default if "day".
#' @param dir One of "fwd", "back" or "both". Default is "both".
#' @param fun Summary function.
#' @import raster lubridate grDevices
#' @return Matrix with statistics for each sample.
#' @details { or list object containing 'x', 'y' and 'time' information}
#' @examples \dontrun{
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dirSample <- function(xy=xy, ot=ot, img=img,  rt=rt, mws=NULL, mwu=NULL, dir=NULL,  metrics=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# check variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {return('error: "xy" is missing')}
  if (!class('xy')%in%c('SpatialPoints', 'SpatialPolygonsDataFrame')) {return('error: "xy" is not of a valid class')}
  
  # sample dates
  if (!exists('ot')) {return('error: "ot" is missing')}
  if (!class(ot)[1]%in%c('Dates', 'POSIXlt'', POSIXct')) {return('error: "ot" is nof of a valid class')}
  if (length(ot)!=length(xy)) {return('error: "xy" and "ot" have different lengths')}
  
  # raster
  if (!exists('img')) {return('error: "img" is missing')}
  if (class(img)[1]!='RasterStackTS') {return('error: "img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {return('error: "xy" and "img" have different projections')}   
  
  # raster dates
  if (!exists('rt')) {return('error: "rt" is missing')}
  if (!class(rt)[1]%in%c('Dates', 'POSIXlt'', POSIXct')) {return('error: "rt" is nof of a valid class')}
  if (length(rt)!=blayers(img)) {return('error: "img" and "rt" have different lengths')}
  
  # time information
   if (is.null(mws)) {return('error: "mws" is missing')} else {
     if (!is.numeric(mws)) {return('error: "mws" us not numeric')}
     if (!mwu%in%c/'day', 'month', 'year') {return('error: "mwu" is not a recongnized keyword')}
    if (mwu=='day') {mwu<-day}
    if (mwu=='month') {mwu<-month}
    if (mwu=='year') {mwu<-year}}
  
  # query type
  if (!is.null(dir)) {
    if (length(dir)>1) {return('error: "dir" has too many entries')}
    if (!dir%in%c('fwd', 'back', 'both')) {return('error: "dir" is not a valid entry')}
  } else {dir <- 'both'}
  
  # check/define input metrics
  if (is.null(fun)) {fun <- function(y) {lm(y~ts, na.exclude=T)$coefficients[2]}} else {
    if(!is.function(fun)) {return('error: "fun" is not a valid function')}}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# query samples
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # backwards sampling
  if (dir=='bcw') {
    qf <- function(x) {
      st <- ot[x]
      mwu(st) <- mwu(st)-mws # start time
      et <- ot[x]
      mwu(ot) <- mwu(ot)+mws # end time
      ind <- which(rt >= st[x] & rt <= et[x])
      v <- extract(xy[x,], img[,,ind])
      uv <- which(!is.na(v))
      v <- v[uv]
      ts <- rt[uv]
      return(metrics(v))}
    return(as.numeric(unlist(sapply(1:nrow(xy), qf))))
  }
  
  # forward sampling    
  if (dir=='fwd') {
    qf <- function(x) {
      et <- ot[x]
      mwu(ot) <- mwu(ot)+mws # end time
      ind <- which(rt >= ot[x] & rt <= et[x])
      v <- extract(xy[x,], img[,,ind])
      uv <- which(!is.na(v))
      v <- v[uv]
      ts <- rt[uv]
      return(metrics(v))}
    return(as.numeric(unlist(sapply(1:nrow(xy), qf))))
  }
    
  # dual sampling
  if (dir=='both') {
    qf <- function(x) {
      st <- ot[x]
      mwu(st) <- mwu(st)-mws # start time
      et <- ot[x]
      mwu(ot) <- mwu(ot)+mws # end time
      ind <- which(rt >= st[x] & rt <= et[x])
      v <- extract(xy[x,], img[,,ind])
      uv <- which(!is.na(v))
      v <- v[uv]
      ts <- rt[uv]
      return(metrics(v))}
    return(as.numeric(unlist(sapply(1:nrow(xy), qf))))
  }

#-------------------------------------------------------------------------------------------------------------------------------#
  
}