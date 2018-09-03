#' @title checkOverlap
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Reports on how much two spatial objects overlap.
#' @param x A spatial object.
#' @param y A spatial object.
#' @return A two element numeric \emph{vector}.
#' @importFrom raster extent intersect ncell raster
#' @details {Uses \link[raster]{intersect} to report on the percentage of the area
#' of \emph{x} and \emph{y} that coincides with their common spatial overlap.}
#' @examples {
#'
#'  require(raster)
#'
#'  # load example probability image
#'  file <- system.file('extdata', 'probabilities.tif', package="rsMove")
#'  img <- raster(file)
#'
#'  # load area of interest
#'  file <- system.file('extdata', 'roi.shp', package="rsMove")
#'  roi <- shapefile(file)
#'
#'  # extract samples
#'  checkOverlap(img, roi)
#'
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

checkOverlap <- function(x, y) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  r1 <- tryCatch(extent(x), error=function(e) return(FALSE)) # is x a spatial object?
  r2 <- tryCatch(extent(y), error=function(e) return(FALSE)) # is y a spatial object?
  if (isFALSE(r1) & isTRUE(r2)) {return(warning('"x" is not a spatial object'))}
  if (isTRUE(r1) & isFALSE(r2)) {return(warning('"y" is not a spatial object'))}
  if (isFALSE(r1) & isFALSE(r2)) {return(warning('neither "x" and "y" are spatial object'))}

  oc <- tryCatch(intersect(r1,r2), error=function(e) return(FALSE)) # do they overlap?
  if (isFALSE(oc)) {return(c(0,0))}

  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. Compare input variables
  #-----------------------------------------------------------------------------------------------------------------------------------------------#

  ac <- ncell(raster(extent(oc), res=1)) # how much is the overlap area?
  a1 <- round((ac / ncell(raster(r1, res=1))) * 100) # what proportion of x is in the overlap?
  a2 <- round((ac / ncell(raster(r2, res=1))) * 100) # what proportion of y is in the overlap?

  return(c(a1,a2))

}
