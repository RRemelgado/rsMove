#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param st Object of class \emph{Date} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param rt Object of class \emph{Date} with \emph{img} observation dates.
#' @param type One of \emph{exact} or \emph{nearest}.
#' @param bs Buffer size (unit depends on the raster projection).
#' @param rd Logical. Should the function ignore duplicated pixels? Default if FALSE.
#' @param fun Passes an external function.
#' @import raster rgdal
#' @importFrom stats median
#' @seealso \code{\link{sampleMove}} \code{\link{backSample}}
#' @return A SpatialPointsDataDataFrame.
#' @details {Returns environmental variables from a raster object for a given set of x and y coordinates.
#'          A buffer size (\emph{bs}) and a user defined function (\emph{fun}) can be specified to sample within an 
#'          area. The defaut is to estimate a weighted mean. If acquisition times are provided (\emph{rt}) the 
#'          raster data is treated as a time series. In this case, the function applies 
#'          one of two sampling approaches: \emph{exact} or \emph{nearest}. If \emph{exact}, the function attempts to map the 
#'          dates of the raster time series with the observation dates of the samples (\emph{ot}). If nearest, 
#'          it searches for the nearest time step. If \emph{rd} is set, the function will account for duplicated 
#'          pixels. The samples will be transposed to pixel coordinates and, for each unique pixel, median 
#'          coordinates will be estimated for each pixel and used to build the output shapefile.}
#' @examples {
#'  
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'  
#'  # read movement data
#'  file <- system.file('extdata', 'konstanz_20130805-20130811.shp', package="rsMove")
#'  moveData <- shapefile(file)
#'  
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'  rsStk <- stack(rsStk, rsStk, rsStk) # dummy files for the example
#'  
#'  # raster dates
#'  r.date <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'  
#'  # sample dates
#'  o.date <- as.Date(moveData@data$date)
#'  
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=moveData, st=o.date, img=rsStk, rt=r.date, type='nearest')
#' 
#' }
#' 
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(xy=xy, st=NULL, img=img, rt=NULL, type=NULL, bs=NULL, rd=FALSE, fun=NULL) {
  
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
  
  # duplicate removal
  if (!is.logical(rd)) {stop('"rd" is not a valid logical argument')}
  if (rd) {
    ext <- extent(img)
    nr <- dim(img)[1]
    pxr <- res(img)[1]
    sp <- (round((ext[4]-xy@coords[,2])/pxr)+1) + nr * round((xy@coords[,1]-ext[1])/pxr)
    dr <- !duplicated(sp)
    op <- crs(xy)
    up <- sp[dr]
    xr <- vector('numeric', length(up))
    yr <- vector('numeric', length(up))
    for (r in 1:length(up)) {
      ind <- which(sp==up[r])
      xr[r] <- median(xy@coords[ind,1])
      yr[r] <- median(xy@coords[ind,2])
    }
    xy <- cbind(xr, yr)
  } else {
    ns <- length(xy)
    op <- crs(xy)
    xy <- xy@coords}
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # extract environmental data
  if (processTime) {
    
    # function to determine target indices
    ifun <- function(x) {
      if (type=='exact') {
        diff <- abs(x-rt)
        return(which(diff==min(diff)[1]))}
      if (type=='nearest') {
        loc <- which(rt==x)
        if (length(loc)>0) {return(loc[1])} else {return(NA)}}}
    
    # retrieve indices per sample
    ind <- sapply(st, ifun)
    ui <- unique(ind)
    
    # output variables
    orv <- vector('numeric', ns)
    odv <- vector('numeric', ns)
    class(odv) <- 'Date'
    
    # query data
    for (i in 1:length(ui)) {
      loc <- which(ind==ui[i])
      if(!is.na(ui[i])) {
        orv[loc] <- extract(img[[ui[i]]], xy[loc,], buffer=bs, fun=fun, na.rm=T)
        odv[loc] <- rt[ui[i]]
      } else {
        orv[loc] <- NA
        odv[loc] <- NA}}
    
    # derive output
    return(SpatialPointsDataFrame(xy, as.data.frame(date=odv, value=orv), proj4string=op))
    
  } else {
    
    # simple query
    orv <- extract(img, xy, buffer=bs, fun=fun, na.rm=T)
    return(SpatialPointsDataFrame(xy, as.data.frame(value=orv), proj4string=op))
    
  }
}