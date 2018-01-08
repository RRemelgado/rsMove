#' @title poly2sample
#'
#' @description {Converts a raster grid to points depending on how much each pixel is covered by a polygon.}
#' @param pol.shp Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonDataFrame}.
#' @param ref.ext Object of class \emph{Extent} or a \emph{raster} object from which an extent can be derived.
#' @param min.cover Minimum percent a pixel should be covered by a polygon for sampling (0-100). Default is 100.
#' @param pixel.res Pixel resolution. Required if \emph{ref.ext} is an \emph{Extent} object. Unit depends on spatial projection.
#' @import sp raster rgdal
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{SpatialPointsDataFrame} with sampled pixels reporting on polygon percent coverage.
#' @details {\emph{poly2Sample} extends on the \code{\link[raster]{rasterize}} function from the raster package making it more efficient
#' over large areas and converting its output into point samples rather than a raster object. For each polygon in (\emph{"pol.shp"}),
#' \emph{poly2sample} extracts the overlapping pixels derived from \emph{ref.ext}. Then, for each pixel, the function estimates the
#' percentage of it that is covered by the reference polygon. Finnally, the function extracts coordinate pairs for pixels that has a
#' percent coverage equal to or greater than \emph{min.cover}.}
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
#'  samples <- poly2sample(pol.shp=roi, ref.ext=img)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------#

poly2sample <- function(pol.shp=pol.shp, ref.ext=NULL, min.cover=NULL, pixel.res=NULL) {

#-------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-------------------------------------------------------------------------------------------------------------------------#

  # check shapefile
  if(is.null('pol.shp')) {stop('"pol.shp" is missing')}
  if(!class(pol.shp)[1]%in%c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
    stop('"pol.shp" is not a valid input.')}

  # check pixel resolution
  if (is.null(pixel.res)) {stop('"pixel.res"is missing')}
  if (!is.numeric(pixel.res)) {stop('"pixel.res" is not numeric')}

  # check/derive reference extent
  if (!is.null(ref.ext)) {
    if (class(ref.ext)!='Extent') {stop('"ref.ext" is not a valid input')}
  } else {ref.ext <- extent(pol.shp)}

  # check cover value
  if (is.null(min.cover)) {min.cover <- 100}
  if (min.cover < 0 | min.cover > 100) {stop('"min.cover" should be between 0 and 100')}

  # build reference raster
  rr <- raster(ref.ext, res=pixel.res, crs=crs(pol.shp))

#-------------------------------------------------------------------------------------------------------------------------#
# 2. evaluate polygons
#-------------------------------------------------------------------------------------------------------------------------#

  lf <- function(i) {
    r <- crop(rr, extent(pol.shp[i,]))
    r <- rasterToPoints(rasterize(pol.shp[i,], r, getCover=T))
    ind <- which(r[,3] > 0) # usable pixels
    return(data.frame(x=r[ind,1], y=r[ind,2], c=r[ind,3]))}
  df0 <- do.call(rbind, lapply(1:length(pol.shp), lf))

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
  return(SpatialPointsDataFrame(df0[,1:2], df0, proj4string=crs(pol.shp)))

#------------------------------------------------------------------------------------------------------------------------#

}
