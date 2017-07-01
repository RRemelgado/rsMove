#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param xy Object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param ot OValid "POSIXlt" or "POSIXct" object with the same length as "xy".
#' @param img Object of class "RasterLayer", RasterStack", "RasterBrick" or "RasterStackTS".
#' @param type One of "exact" or "nearest".
#' @param bs Buffer size (unit depends on the raster projection).
#' @param remove.duplicates Logical. Default if FALSE.
#' @param fun Passes an external function.
#' @import raster rts grDevices
#' @return Matrix or vector of environmental variables for each sample.
#' @details {Returns environmental variables from a raster object for a given set of x and y coordinates.
#'          A buffer size ("bs") and a user defined function ("fun") can be specified to sample within an area.
#'          The defaut is to estimate a weighted mean. If a raster time series is provided the function applies 
#'          one of two sampling approaches: "exact" or "nearest". If "exact", the function attempts to map the 
#'          dates of the raster time series with the observation dates of the samples ("ot"). If nearest, 
#'          it searches for the nearest time step.}
#' @examples \dontrun{
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(xy, img=img, ot=NULL, type=NULL, bs=NULL, remove.duplicates=F, fun=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# check variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('error: "xy" is missing')}
  if (!class('xy')%in%c('SpatialPoints', 'SpatialPolygonsDataFrame')) {stop('error: "xy" is not of a valid class')}
  
  # raster
  if (!exists('img')) {stop('error: "img" is missing')}
  if (!class(img)[1]%in%c('RasterStack', 'RasterBrick', 'RasterStackTS')) {stop('error: "img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('error: "xy" and "img" have different projections')}   
  
  # check if raster is a ts
  if (class(img)[1]=='RasterStackTS') {
    if (is.null(ot)) {stop('error: "ts" object provided. "ot" is missing')}
    if (!class(ot)[1]%in%c('Date', 'POSIXlt'', POSIXct')) {stop('error: "ot" is not of a valid class')}
    if (length(ot)!=length(xy)) {stop('error: lengths of "ot" and "xy" differ')}
    processTime <- TRUE
  } else {processTime <- FALSE}
  
  # auxiliary
  if (!is.null(bs)) {if (!is.numeric(bs)) {stop('error: "bs" assigned but not numeric')}} else {fun=NULL}
  if (!is.null(bs) & is.null(fun)) {fun <- function(x) {sum(x*x) / sum(x)}} else {
    if (!is.function(fun)) {stop('error: "fun" is not a valid function')}
    tf <- fun(c(1,2,3)) # test function
    if (length(tf)!=1) {stop('error: provided function returns more than 1 output')}
  }
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # extract environmental data
  if (processTime) {
    
    r.dates <- index(img) # extract rts object dates
    s.dates <- as.Date(ot) # convert sample timestamp to "Date"
    
    # function for exact date query
    qf <- if (type=='exact') {
      function(x) {
        ind <- which(r.dates==s.dates[x])
        if (length(ind)>1) {return(extract(img[[ind]], xy@coords[x], buffer=bs, fun=fun, na.rm=T))}
      }
    }
    
    # function for nearest date query
    qf <- if (type=="nearest") {
      function(x) {
        td <- abs(s.dates[x]-r.dates[x])
        ind <- which(td==min(td))
        return(extract(img[[ind]], xy@coords[x], buffer=bs, fun=fun, na.rm=T))
      }
    }
    
    # apply query function
    ov <- unlist(sapply(1:nrow(xy), qf))
    
  } else {
    
    # simple query
    ov <- extract(img, xy@coords, buffer=bs, fun=fun, na.rm=T)
  
  }
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # return output
  return(ov)
  
}