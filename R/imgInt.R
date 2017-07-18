#' @title imgInt
#'
#' @description Temporal linear interpolation of raster data.
#' @param img Object of class \emph{RasterStack} or \emph{RasterBrick}.
#' @param td Target dates. Object of class \emph{Date}.
#' @param rd Raster dates. Object of class \emph{Date}.
#' @param bs Temporal buffer size (in days).
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @import raster sp rgdal
#' @importFrom stats lm
#' @seealso @seealso \code{\link{dataQuery}} \code{\link{timeDirSample}} \code{\link{spaceDirSample}}
#' @return A \emph{RasterBrick} or a \emph{data frame}.
#' @details {Performs a pixel-wise linear interpolation over a raster 
#' for a given set of dates (\emph{td}). A teporal buffer (\emph{bs}) is required 
#' to limit the search for reference data points (\emph{rd}). If \emph{xy} is 
#' provided the function only considers the pixels that overlap with 
#' the these sample points. Otherwise, a RasterBrick is provided.}
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
#'  rd = seq.Date(as.Date("2012-01-01"), as.Date("2012-12-31"), 45)
#'  
#'  # target dates
#'  td = seq.Date(as.Date("2012-04-01"), as.Date("2012-08-31"), 45)
#'  
#'  # interpolate raster data to target dates
#'  i.img <- imgInt(img=rsStk, rd=rd, td=td, bs=60, xy=moveData)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------#

imgInt <- function(img=img, rd=rd, td=td, bs=NULL, xy=NULL) {
  
#-------------------------------------------------------------------------------------------------#
# 1. check input variables  
#-------------------------------------------------------------------------------------------------#
  
  # check variables
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
  if (!is.null(xy)) {
    if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {
      stop('"shp is nor a valid point shapefile object')}
      if (crs(xy)@projargs!=crs(img)@projargs) {
        stop('"xy" and "img" have different projections')}}
  if (!exists('td')) {stop('"td" is missing')}
  if(class(td)!='Date') {stop('"td" is not a "Date" object')}
  if (!exists('rd')) {stop('"rd" is missing')}
  if(class(rd)!='Date') {stop('"rd" is not a "Date" object')}
  if (is.null(bs)) {stop('"bs" is missing')}
  if (!is.numeric(bs)) {stop('"bs" is not numeric')}
  
#-------------------------------------------------------------------------------------------------#
# 2. interpolate values
#-------------------------------------------------------------------------------------------------#
  
  # interpolation function
  int <- function(x) {
    di <- which(rd==otd & !is.na(x))
    if (length(di)>0) {return(mean(x[di]))} else {
      bi <- rev(which(!is.na(x) & rd < otd & rd >= (otd-bs)))
      ai <- which(!is.na(x) & rd > otd & rd <= (otd+bs))
      if (length(bi)>=1 & length(ai)>=1) {
        lc <- lm(c(x[bi[1]],x[ai[1]])~as.numeric(c(rd[bi[1]],rd[ai[1]])))
        return(as.numeric(otd)*lc$coefficients[2]+lc$coefficients[1])
      } else {return(NA)}}}
    
  # apply function
  np <- nrow(iData)
  if (is.null(xy)) {
    out<-brick(img[[1]], nl=length(td))
    for (r in 1:length(td)) {
      otd <- td[r]
      out[[r]] <- calc(img, int)}}
  if (!is.null(xy)) {
    idata <- extract(img, xy)
    out<-data.frame(date=td, matrix(0,length(td),np))
    for (r in 1:nrow(idata)) {
      out[r,2:ncol(out)] <- apply(iData, 1, int)}}

  # provide output
  return(out)
  
}