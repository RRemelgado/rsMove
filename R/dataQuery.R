#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param st Object of class \emph{Date} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param rt Object of class \emph{Date} with \emph{img} observation dates.
#' @param tb Two element vector with temporal search buffer, expressed in days.
#' @param type One of \emph{exact} or \emph{nearest}.
#' @param bs Buffer size (unit depends on the raster projection).
#' @param remove.dup Logical. Should the function ignore duplicated pixels? Default if FALSE.
#' @param fun Passes an external function.
#' @import raster rgdal
#' @importFrom stats median
#' @seealso \code{\link{sampleMove}} \code{\link{backSample}}
#' @return A SpatialPointsDataDataFrame.
#' @details {Returns environmental variables from a raster object for a given set of x and y coordinates.
#'          A buffer size (\emph{bs}) and a user defined function (\emph{fun}) can be specified to sample 
#'          within an area. The defaut is to estimate a weighted mean. If acquisition times are provided 
#'          (\emph{rt}) the raster data is treated as a time series. In this case, the function applies 
#'          one of two sampling approaches: \emph{exact} or \emph{nearest}. If \emph{exact}, the function 
#'          attempts to map the dates of the raster time series with the observation dates of the samples 
#'          (\emph{ot}). If nearest, it searches for the nearest time step. In this case, a temporal buffer, 
#'          defined by \emph{tb}, can be defined to restric the search. \emph{tb} passes two values which 
#'          represent the size of the buffer in two directions: before and after the target date. This allows 
#'          for bacward of forward sampling. If \emph{remove.dup} is set, the function will account for duplicated 
#'          pixels. The samples will be transposed to pixel coordinates and, for each unique pixel, median 
#'          coordinates will be estimated for each pixel and used to build the output shapefile.}
#' @examples {
#'  
#'  require(raster)
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
#'  rsQuery <- dataQuery(xy=moveData, st=o.date, img=rsStk, rt=r.date, tb=c(30,30), type='nearest')
#' 
#' }
#' 
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(xy=xy, st=NULL, img=img, rt=NULL, tb=NULL, type=NULL, bs=NULL, remove.dup=FALSE, fun=NULL) {
  
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
  if (!type%in%c('nearest', 'exact')) {stop('"type" is not a recognized keyword')}
  if (!is.null(tb)) {
    if (!is.numeric(tb)) {stop('"tb" is not numeric')}
    if (length(tb)!=2) {stop('"tb" should be a two element vector')}}
  
  # duplicate removal
  if (!is.logical(remove.dup)) {stop('"remove.dup" is not a valid logical argument')}
  if (remove.dup) {
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
    op <- crs(xy)
    xy <- xy@coords}
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # read data
  edata <- as.data.frame(extract(img, xy, buffer=NULL, fun=fun, na.rm=T))
  
  # extract environmental data
  if (processTime) {
    
    # number of samples
    ns <- nrow(xy)
    
    if (type=='nearest') {
      
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
      return(SpatialPointsDataFrame(xy, data.frame(value=orv, date=ord), proj4string=op))
      
    }
    
    if (type=='exact') {
      
      # function to extract values
      qf <- function(i) {
        ind <- which(!is.na(edata[i,]))
        if (length(ind)!=0) {
          diff <- abs(st[i]-rt[ind])
          loc <- which(diff==0)[1]
          if (length(loc) > 0) {return(edata[i,ind[loc]])
            } else {return(NA)}} else {return(NA)}}
      orv <- as.numeric(sapply(1:ns, qf))
      
      # derive shapefile
      return(SpatialPointsDataFrame(xy, data.frame(value=orv), proj4string=op))}
      
  } else {
    
    # simple query
    return(SpatialPointsDataFrame(xy, as.data.frame(edata), proj4string=op))
    
  }
}