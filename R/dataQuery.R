#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param st Object of class \emph{Date} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param rt Object of class \emph{Date} with \emph{img} observation dates.
#' @param tb Two element vector with temporal search buffer, expressed in days.
#' @param bs Buffer size (unit depends on the raster projection).
#'  @param fun Passes an external function.
#' @import raster rgdal
#' @importFrom stats median
#' @seealso \code{\link{sampleMove}} \code{\link{backSample}}
#' @return A n object of class \emph{vector} or \emph{data.frame}.
#' @details {Returns environmental variables from a raster object for a given set of x and y coordinates.
#'          A buffer size (\emph{bs}) and a user defined function (\emph{fun}) can be specified to sample 
#'          within an area. The defaut is to estimate a weighted mean. If raster acquisition times are provided 
#'          (\emph{rt}) and the date of sampling (\emph{st}). In this case, the function will treat the raster 
#'          data as a time series and search for clear pixel in time within the contraints of a temporal buffer 
#'          defined by \emph{tb}. \emph{tb} passes two values which represent the size of the buffer in two 
#'          directions: before and after the target date. This allows for bacward and forward sampling.}
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
#'  r.date <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'  
#'  # sample dates
#'  o.date <- as.Date(moveData@data$date)
#'  
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=moveData, st=o.date, img=rsStk, rt=r.date, tb=c(30,30))
#' 
#' }
#' 
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(xy=xy, st=NULL, img=img, rt=NULL, tb=NULL, bs=NULL, fun=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# check variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  
  # raster
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]%in%c('RasterLayer','RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}   
  
  # check if raster is a ts
  if (!is.null(rt)) {
    if (class(rt)[1]!='Date') {stop('"rt" is not of a valid class')}
    if (length(rt)!=nlayers(img)) {stop('lengths of "rt" and "img" differ')}
    if (is.null(st)) {stop('"st" is missing')}
    if (class(st)[1]!='Date') {stop('"st" is not of a valid class')}
    if (length(st)!=length(xy)) {stop('lengths of "st" and "xy" differ')}
    processTime <- TRUE
  } else {processTime <- FALSE}
  
  # auxiliary
  if (!is.null(bs)) {if (!is.numeric(bs)) {stop('"bs" assigned but not numeric')}} else {fun=NULL}
  if (!is.null(bs) & is.null(fun)) {fun <- function(x) {sum(x*x) / sum(x)}} else {
  if (!is.null(fun)) {if (!is.function(fun)) {stop('"fun" is not a valid function')}}}
  if (!is.null(tb)) {
    if (!is.numeric(tb)) {stop('"tb" is not numeric')}
    if (length(tb)!=2) {stop('"tb" should be a two element vector')}}
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # output projection
  op <- crs(xy)
  
  # read data
  edata <- as.data.frame(extract(img, xy@coords, buffer=bs, fun=fun, na.rm=T))
  
  # extract environmental data
  if (processTime) {
    
    # number of samples
    ns <- nrow(xy@coords)
    
    # function to select pixels
    if (is.null(tb)) {
      qf <- function(i) {
        ind <- which(!is.na(edata[i,]))
        if (length(ind)!=0) {
          diff <- abs(st[i]-rt[ind])
          loc <- which(diff==min(diff))[1]
          return(list(value=edata[i,ind[loc]], date=rt[ind[loc]]))
        } else{return(list(value=NA, date=NA))}}
    } else {
      qf <- function(i) {
        ind <- which(!is.na(edata[i,]))
        if (length(ind)!=0) {
          diff <- abs(st[i]-rt[ind])
          loc <- rt[ind] >= (st[i]-tb[1]) & rt[ind] <= (st[i]+tb[2])
          if (sum(loc)>0) {
            loc <- which(diff==min(diff[loc]))[1]
            return(list(value=edata[i,ind[loc]], date=rt[ind[loc]]))
          } else {return(list(value=NA, date=NA))}}
        else {return(list(value=NA, date=NA))}}}
    
    # retrieve values
    edata <- lapply(1:ns, qf)
    orv <- do.call('c', lapply(edata, function(x) {x$value}))
    ord <- do.call('c', lapply(edata, function(x) {x$date}))
    
    # derive shapefile
    return(data.frame(value=orv, date=ord))
      
  } else {
    
    # simple query
    return(edata)
    
  }
}