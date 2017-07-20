#' @title moveSeg
#'
#' @description Remote sensing based point segmentation
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param edata Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param type Raster data type. One of \emph{cont} (continues) or \emph{cat} (for categorical).
#' @param threshold Change threshold.
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates. 
#' @param r.fun Raster summary function. Default is pca.
#' @param s.fun Output summary function. Default is mean.
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
#' \emph{theshold} and \emph{s.fun} keywords. In this case, a new segment is identified 
#' if the any change is observed between two consecutife points. The output consists of 
#' a list containing a \emph{SpatialPointsDataFrame} (\emph{$points}) reporting on the 
#' segment ID (\emph{sid}) associated to each sample and a data frame (\emph{$report}) 
#' with the amount of points in each region and the value returned by \emph{s.fun}. If 
#' \emph{ot} is provided, the function also provides the elapsed time within each segment. 
#' If \emph{fun} is set by the user, the provided function will be used to summarize the 
#' raster values at each segment. Also, if \emph{edata} is a \emph{RasterStack} or a 
#' \emph{RasterBrick}, \emph{r.fun} is used to reduce the multi-layered object to a single 
#' layer. he user can either use a Principal Component Analysis (PCA) analysis by setting 
#' \emph{r.fun} to \emph{pca} or provide a function. If \emph{pca} is selected, the first 
#' Principal Conponent (PC) \emph{threshold} will be set automatically to the standard 
#' deviation of the first PC. This is the default. If \emph{method} is \emph{cat} the 
#' function will assume the data is categorical and wil define segments whenever a change occurres.}
#' @examples {
#'  
#'  require(raster)
#'  
#'  # read movement data
#'  moveData <- shapefile(system.file('extdata', 'konstanz_20130804.shp', package="rsMove"))
#'  
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'  
#'  # perform directional sampling
#'  seg <- moveSeg(xy=moveData, edata=r, type="cont", threshold=0.1)
#'  
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------#

moveSeg <- function(xy=xy, edata=edata, type='cont', ot=NULL, threshold=0.1, r.fun=NULL, s.fun=NULL) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input variables  
#---------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy) # output projection
  
  # sample dates
  if (!is.null(ot)) {
    if (!class(ot)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"ot" is nof of a valid class')}
    if (length(ot)!=length(xy)) {stop('errorr: "xy" and "ot" have different lengths')}}
  
  # environmental data
  if (class(edata)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    if (crs(xy)@projargs!=crs(edata)@projargs) {stop('"xy" and "edata" have different projections')}
    nvar <- nlayers(edata)
    } else {
    if (!class(edata)[1]%in%c('data.frame')) {stop('"edata" is neither a raster or a data frame')}
    if (nrow(edata)!=length(xy)) {stop('number of elements in "xy" and "edata" do not match')}
      nvar <- ncol(edata)}
  if (type=='cat' & nvar>1) {stop(paste0('type was set to ', type, '. raster should be a single layer'))}
  
  # check threshold
  if (!is.numeric(threshold)) {stop('"threshold" is not numeric')}
  
  # check query type
  if (!type%in%c('cont', 'cat')) {stop('type is not  avalid keyword')}
  if (type=='cont') {
    if(nvar>1) {
      if (is.null(r.fun)) {r.fun <- 'pca'}
      if (!is.null(r.fun) & r.fun!='pca') {if (!is.function(r.fun)){stop('"r.fun" is not a valid keyword or function')}}}
    if (!is.null(s.fun)) {if (!is.function(s.fun)){stop('"s.fun" is not a valid keyword or function')}}
    if (is.null(s.fun)) {s.fun <- function(x) {return(mean(x, na.rm=T))}}}
  if (type=='cat') {
    threshold <- 1
    s.fun <- function(x){return(x[1])}}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. identify segments
#---------------------------------------------------------------------------------------------------------------------#
  
  # read rs data
  if (class(edata)[1]!='data.frame') {edata <- as.data.frame(extract(edata, xy))}
  
  # summarize data (if needed)
  if (nvar>1) {
    if (is.function(r.fun)) {edata <- as.data.frame(apply(edata, 1, r.fun))}
    if (r.fun=='pca') {
      edata <- prcomp(edata, scale=T, center=T)
      threshold <- edata$sdev[1]
      edata <- as.data.frame(edata$x[,1])}}
  
  # search for segments
  r0 <- 1
  li <- 1
  id <- list() # segment id
  rv <- list() # segment value
  
  for (r in 2:length(xy)) {
    diff <- abs(edata[r,1]-edata[(r-1),1])
    if (!is.na(diff)) {
      if (diff >= threshold) {
        ep <- r-1
        rv[[li]] <- s.fun(edata[c(r0:ep),1])
        id[[li]] <- replicate(length(c(r0:ep)), li)
        r0 <- r
        li <- li + 1
        if (r==length(xy)) {
          id[[li]] <- li
          rv[[li]] <- edata[r,1]}
      } else {if (r==length(xy)) {
        ep <- r
        rv[[li]] <- s.fun(edata[c(r0:ep),1])
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