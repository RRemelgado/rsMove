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
#' @import raster grDevices
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
#' @examples \dontrun{
#' 
#'  # read movement data
#'  moveData <- shapefile(system.file('extdata', 'konstanz_20130805-20130811.shp', package="rsMove"))
#'  
#'  # read remote sensing data
#'  rsStack <- stack(list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=T))
#'  
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=moveData,img=rsStack)
#' 
#' }
#' 
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(xy=xy, st=NULL, img=img, rt=NULL, type=NULL, bs=NULL, rd=F, fun=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# check variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('error: "xy" is missing')}
  if (!class(xy)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('error: "xy" is not of a valid class')}
  
  # raster
  if (!exists('img')) {stop('error: "img" is missing')}
  if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('error: "img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('error: "xy" and "img" have different projections')}   
  
  # check if raster is a ts
  if (!is.null(rt)) {
    if (class(rt)[1]!='Date') {stop('error: "rt" is not of a valid class')}
    if (length(rt)!=length(xy)) {stop('error: lengths of "rt" and "xy" differ')}
    if (is.null(st)) {stop('error: "st" is missing')}
    if (class(st)[1]!='Date') {stop('error: "st" is not of a valid class')}
    if (length(st)!=length(xy)) {stop('error: lengths of "st" and "xy" differ')}
    processTime <- TRUE
  } else {processTime <- FALSE}
  
  # auxiliary
  if (!is.null(bs)) {if (!is.numeric(bs)) {stop('error: "bs" assigned but not numeric')}} else {fun=NULL}
  if (!is.null(bs) & is.null(fun)) {fun <- function(x) {sum(x*x) / sum(x)}} else {
  if (!is.null(fun)) {if (!is.function(fun)) {stop('error: "fun" is not a valid function')}}}
  
  # duplicate removal
  if (!is.logical(rd)) {stop('error: "rd" is not a valid logical argument')}
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
  } else {xy <- xy@coords}
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # extract environmental data
  if (processTime) {
    
    # function for exact date query
    qf <- if (type=='exact') {
      function(x) {
        ind <- which(rt==st[x])
        if (length(ind)>1) {return(extract(xy[x,], img[[ind]], buffer=bs, fun=fun, na.rm=T))}
      }
    }
    
    # function for nearest date query
    qf <- if (type=="nearest") {
      function(x) {
        td <- abs(st[x]-rt[x])
        ind <- which(td==min(td))
        return(extract(img[[ind]], xy[x,], buffer=bs, fun=fun, na.rm=T))
      }
    }
    
    # apply query function
    ov <- unlist(sapply(1:nrow(xy), qf))
    
  } else {
    
    # simple query
    ov <- extract(img, xy, buffer=bs, fun=fun, na.rm=T)
  
  }
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # return output
  return(SpatialPointsDataFrame(xy, as.data.frame(ov), proj4string=op))
  
}