#' @title imgInt
#'
#' @description Temporal linear interpolation of raster data.
#' @param img Object of class \emph{RasterStack} or \emph{RasterBrick}.
#' @param td Target dates. Object of class \emph{Date}.
#' @param rd Raster dates. Object of class \emph{Date}.
#' @param bs Temporal buffer size (in days).
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param edata Object of class \emph{data frame} or \emph{matrix} with remote sensing data.
#' @import raster sp rgdal
#' @importFrom stats lm
#' @seealso @seealso \code{\link{dataQuery}} \code{\link{timeDir}} \code{\link{spaceDir}} \code{\link{moveSeg}}
#' @return A \emph{RasterBrick} or a \emph{data frame}.
#' @details {Performs a pixel-wise linear interpolation over a raster 
#' for a given set of dates (\emph{td}). A teporal buffer (\emph{bs}) is required 
#' to limit the search for reference data points (\emph{rd}). If \emph{xy} is 
#' provided the function only considers the pixels that overlap with 
#' the these sample points. Otherwise, a RasterBrick is provided. However, 
#' if \emph{edata} is provided, \emph{xy} and \emph{img} are ignored. and a 
#' data frame is provided.}
#' @examples {
#'  
#'  require(raster)
#'  
#'  # read movement data
#'  file <- system.file('extdata', 'konstanz_20130805-20130811.shp', package="rsMove")
#'  moveData <- shapefile(file)[1:10,]
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
#'  td = as.Date("2012-04-01")
#'  
#'  # interpolate raster data to target dates
#'  i.img <- imgInt(img=rsStk, rd=rd, td=td, bs=60, xy=moveData)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------#

imgInt <- function(img=NULL, rd=rd, td=td, bs=NULL, xy=NULL, edata=NULL) {
  
#-------------------------------------------------------------------------------------------------#
# 1. check input variables  
#-------------------------------------------------------------------------------------------------#
  
  # check variables (if edata is provided)
  if(class(td)!='Date') {stop('"td" is not a "Date" object')}
  if(class(rd)!='Date') {stop('"rd" is not a "Date" object')}
  if (!is.null(edata)) {
    if (class(edata)[1]!='data.frame') {stop('"edata" is neither a matrix or a data frame')}
    if (ncol(edata)!=length(rd)) {stop('mismatch in "edata" and "rd" dimensions')}
    xy <- NULL
    img <- NULL
    bs <- NULL
  } else {
    if (is.null(img)) {stop('"img" is missing')}
    if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
    if (!is.null(xy)) {
      if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {
        stop('"shp is nor a valid point shapefile object')}
      if (crs(xy)@projargs!=crs(img)@projargs) {
        stop('"xy" and "img" have different projections')}}
    if (nlayers(img)!=length(rd)) {stop('length of "img" and "rd" do not match')}
    if (!is.null(bs)) {if (!is.numeric(bs)) {stop('"bs" is not numeric')}}}
  
#-------------------------------------------------------------------------------------------------#
# 2. interpolate values
#-------------------------------------------------------------------------------------------------#
  
  # interpolation function
  otd <- ''
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
  if (is.null(xy) & !is.null(img)) {
    out<-brick(img[[1]], nl=length(td))
    for (r in 1:length(td)) {
      otd <- td[r]
      out[[r]] <- calc(img, int)}}
  if (!is.null(xy)) {
    if (is.null(edata)) {edata <- extract(img, xy)}
    out <- matrix(0, length(xy), length(td))
    colnames(out) <- as.character(td)
    for (r in 1:length(td)) {
      otd <- td[r]
      out[,r] <- apply(edata, 1, int)}}

  # provide output
  return(out)
  
}