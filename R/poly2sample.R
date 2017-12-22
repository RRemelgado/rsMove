#' @title poly2sample
#'
#' @description {Selection of samples from a reference grid based on a per 
#' pixel estimate of the percentage of the pixel covered by polygons.}
#' @param pol Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonDataFrame}.
#' @param ref.ext Object of class \emph{Extent} or a raster object from which an extent can be derived.
#' @param min.cover Minimum pixel cover (0-100). Default is 100.
#' @param pixel.res Pixel resolution.
#' @import sp raster rgdal
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{SpatialPointsDataFrame}.
#' @details {For each pixel overlaped by polygons (\emph{"pol"}) the functions returns the percent of the 
#' pixel covered and returns the corresponding samples if a minimum percent is reached (\emph{min.cover}). 
#' the function requires a reference extent (\emph{ref.ext}) and pixel resolution (\emph{pixel.res}) over 
#' which this analysis will be performed. This function serves as a wrapper for the function 
#' \code{\link[raster]{rasterize}} from the raster package.}
#' @examples {
#'
#'  require(raster)
#'
#'  # load example probability image
#'  file <- system.file('extdata', 'konstanz_probabilities.tif', package="rsMove")
#'  img <- raster(file)
#'
#'  # load area of interest
#'  file <- system.file('extdata', 'konstanz_roi.shp', package="rsMove")
#'  roi <- shapefile(file)
#'
#'  # segment probabilities
#'  samples <- poly2sample(pol=roi, ref.ext=img)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------#

poly2sample <- function(pol=pol, ref.ext=NULL, min.cover=NULL, pixel.res=NULL) {

#-------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-------------------------------------------------------------------------------------------------------------------------#
  
  # check shapefile
  if(is.null('pol')) {stop('"pol" is missing')}
  if(!class(pol)[1]%in%c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
    stop('"pol" is not a valid input.')}
  
  # check pixel resolution
  if (is.null(pixel.res)) {stop('"pixel.res"is missing')}
  if (!is.numeric(pixel.res)) {stop('"pixel.res" is not numeric')}
  
  # check/derive reference extent
  if (!is.null(ref.ext)) {
    if (class(ref.ext)!='Extent') {stop('"ref.ext" is not a valid input')}
  } else {ref.ext <- extent(pol)}
  
  # check cover value
  if (is.null(min.cover)) {min.cover <- 100}
  if (min.cover < 0 | min.cover > 100) {stop('"min.cover" should be between 0 and 100')}

  # build reference raster
  rr <- raster(ref.ext, res=pixel.res, crs=crs(pol))
  
#-------------------------------------------------------------------------------------------------------------------------#
# 2. evaluate polygons
#-------------------------------------------------------------------------------------------------------------------------#
  
  lf <- function(i) {
    r <- crop(rr, extent(pol[i,]))
    r <- rasterToPoints(rasterize(pol[i,], r, getCover=T))
    ind <- which(r[,3] > 0) # usable pixels
    return(data.frame(x=r[ind,1], y=r[ind,2], c=r[ind,3]))}
  df0 <- do.call(rbind, lapply(1:length(pol), lf))

#-------------------------------------------------------------------------------------------------------------------------#
# 3. remove duplicated pixels and update pecent count
#-------------------------------------------------------------------------------------------------------------------------#
  
  pp <- cellFromXY(rr, df0[,1:2]) # pixel positions
  up <- unique(pp) # unique positions
  pc <- sapply(up, function(x) {sum(df0[which(pp==x),3])}) # update percentages
  pc[which(pc>100)] <- 100 # account for miss-calculations (related to e.g. overlapping polygons)
  xy <- xyFromCell(rr, up) # convert unique positions to coordinates
  df0 <- data.frame(x=xy[,1], y=xy[,2], cover=pc) # build final data frame

  # apply percent cover filter
  ind <- which(pc >= min.cover)
  df0 <- df0[ind,]

  rm(pp, up, pc, xy, ind)
  
  # return output
  return(SpatialPointsDataFrame(df0[,1:2], df0, proj4string=crs(pol)))

#------------------------------------------------------------------------------------------------------------------------#

}
