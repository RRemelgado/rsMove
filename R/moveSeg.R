#' @title moveSeg
#'
#' @description Remote sensing based point segmentation
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param type Raster data type. One of \emph{cont} (continues) or \emph{cat} (for categorical).
#' @param threshold Percent change threshold.
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates. 
#' @param fun Summary function.
#' @import raster rgdal
#' @seealso \code{\link{timeDirSample}} \code{\link{dataQuery}}
#' @return A \emph{list}.
#' @details {Segmentation of a point shapefile based on the spatial variability 
#' of a raster dataset. When the \emph{method} is set to \emph{'cont'}, the raster 
#' data is assumed to be continuous. Then, the function determines the percentual 
#' change between each pair of two consecutive coordinate pairs. If this change is 
#' above a predifined \emph{threshold}, a new pointer is added and the previous 
#' sequence of samples is labeled as a unique segment. If \emph{method} is set as 
#' \emph{'cat'}, the function assumes the raster data is categorical ignoring the 
#' \emph{theshold} keyword. In this case, a new segment is identified if the any 
#' change is observed between two consecutife points. The output consists of a list 
#' containing a \emph{SpatialPointsDataFrame} (\emph{$points}) reporting on the segment 
#' ID (\emph{sid}) associated to each sample and a data frame (\emph{$report}) with the 
#' amount of points in each region and the value returned by \emph{fun}. If \emph{ot} is 
#' provided, the function also provides the elapsed time within each segment. If \emph{fun} 
#' is set by the user, the provided function will be used to summarize the raster values 
#' at each segment. Also, if \emph{img} is a \emph{RasterStack} or a \emph{RasterBrick}, 
#' the \emph{fun} is used to reduce the multi-layered object to a single layer. By default, 
#' the maximum value is used. Like for \emph{threshold}, \emph{fun} is ignored if 
#' \emph{method} is \emph{cat}.}
#' @examples {
#'  
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'  
#'  # read movement data
#'  moveData <- shapefile(system.file('extdata', 'konstanz_20130804.shp', package="rsMove"))
#'  
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'  
#'  # perform directional sampling
#'  seg <- moveSeg(xy=moveData, img=r, type="cont", threshold=0.1)
#'  
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------#

moveSeg <- function(xy=xy, img=img, type='cont', ot=NULL, threshold=0.1, fun=NULL) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input variables  
#---------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  
  # sample dates
  if (!is.null(ot)) {
    if (!class(ot)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"ot" is nof of a valid class')}
    if (length(ot)!=length(xy)) {stop('errorr: "xy" and "ot" have different lengths')}}
  
  # raster
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
  if (type=='cat' & nlayers(img)>1) {stop(paste0('type was set to ', type, '. raster should be a single layer'))}
  
  # compare projections
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}
  
  # change threshold
  if (!is.numeric(threshold)) {stop('"threshold" is not numeric')}

  # summariy function
  if (type=='cont') {
  if (!type%in%c('cont', 'cat')) {stop('type is not  avalid keyword')}
  if (!is.null(fun)) {if (!is.function(fun)) {
    stop('"fun" is not a function')}} else {fun = function(x) {return(max(x, na.rm=T))}}}
  if (type=='cat') {
    threshold <- 1
    fun <- function(x){return(x[1])}}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. identify segments
#---------------------------------------------------------------------------------------------------------------------#
  
  edata <- as.data.frame(extract(img, xy)) # read rs data
  if (nlayers(img)) {edata <- as.data.frame(apply(edata, 1, fun))}
  pxr <- res(img)[1] # pixel resolution
  rProj <- crs(img)
  
  # search for segments
  r0 <- 1
  li <- 1
  id <- list() # segment id
  rv <- list() # segment value
  
  for (r in 2:length(xy)) {
    if (type=='cont') {diff <- abs(1-(edata[r,1]/edata[(r-1),1]))}
    if (type=='cat') {diff <- abs(edata[r,1]-edata[(r-1),1])}
    if (!is.na(diff)) {
      if (diff >= threshold) {
        ep <- r-1
        rv[[li]] <- fun(edata[c(r0:ep),1])
        id[[li]] <- replicate(length(c(r0:ep)), li)
        r0 <- r
        li <- li + 1
        if (r==length(xy)) {
          id[[li]] <- li
          rv[[li]] <- edata[r,1]}
      } else {if (r==length(xy)) {
        ep <- r
        rv[[li]] <- fun(edata[c(r0:ep),1])
        id[[li]] <- replicate(length(c(r0:ep)), li)}}}}
  rv <- unlist(rv)
  id <- unlist(id)
 
#---------------------------------------------------------------------------------------------------------------------#
# 3. build/return output
#---------------------------------------------------------------------------------------------------------------------#

  # update original shapefile
  p.shp <- SpatialPointsDataFrame(xy@coords, data.frame(sid=id), proj4string=rProj)
  
  # build region report
  uid <- sort(unique(id))
  if (!is.null(ot)) {
    f <- function(x) {
      ind <- which(id==x)
      et <- difftime(ot[ind[length(ind)]], ot[ind[1]], units="mins")
      np <- length(ind)
      return(list(time=as.numeric(et), count=np))}
    sstat <- lapply(uid, f)
    df <- data.frame(sid=uid, count=sapply(sstat, function(x) {x$count}), 
          time=sapply(sstat, function(x) {x$time}), value=rv)
  } else {df <- data.frame(sid=uid, count=sapply(uid, function(x){sum(id==x)}))}
  
  # derive output
  return(list(points=p.shp, report=df))
  
}