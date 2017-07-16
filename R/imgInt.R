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
  if (!is.null(xy) {if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}}   
  if (!exists('td')) {stop('"td" is missing')}
  if(class(td)!='Date') {stop('"td" is not a "Date" object')}
  if (!exists('rd')) {stop('"rd" is missing')}
  if(class(rd)!='Date') {stop('"rd" is not a "Date" object')}
  if (is.null(bs)) {stop('"bs" is missing')}
  if (!is.numeric(bs)) {stop('"bs" is not numeric')}
  if (!is.null(xy)) {if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {
    stop('"shp is nor a valid point shapefile object')}}
  
#-------------------------------------------------------------------------------------------------#
# 2. extract data
#-------------------------------------------------------------------------------------------------#
  
  # raster dimensions
  ext <- extent(img) # raster extent
  gd <- dim(img) # raster dimensions  
  op <- crs(img) # output projection
  
  if (!is.null(xy)) {iData <- extract(img, xy)} else {iData <- getValues(img)}
  
#-------------------------------------------------------------------------------------------------#
# 3. interpolate values
#-------------------------------------------------------------------------------------------------#
  
  # convert dates to julian
  td0 <- td # used for the output
  td <- julian(td)
  rd <- julian(rd)
  
  # interpolation function
  int <- function(x) {
    di <- which(rd==td[r])
    if (length(di)>0) {return(mean(x[di]))} else {
      bi <- rev(which(!is.na(x) & rd < td[r] & rd >= (td[r]-bs)))
      ai <- which(!is.na(x) & rd > td[r] & rd <= (td[r]+bs))
      if (length(bi)>0 & length(ai)>0) {
        lc <- lm(c(x[bi[1]],x[ai[1]])~c(rd[bi[1]],rd[ai[1]]))
        return(as.numeric(td[r]*lc$coefficients[2]+lc$coefficients[1]))
      } else {return(NA)}}}
    
  # apply function
  np <- nrow(iData)
  if (is.null(xy)) {
    out<-brick(extent(img), nrow=gd[1], ncol=gd[2], crs=op, nl=length(td))
    for (r in 1:length(td)) {out[[r]] <- setValues(img[[1]], apply(iData, 1, int))}
  } else {
    out<-data.frame(date=td0, matrix(0,length(td0),np))
    for (r in 1:length(td)) {out[r,2:ncol(out)] <- apply(iData, 1, int)}}

  # provide output
  return(out)
  
}