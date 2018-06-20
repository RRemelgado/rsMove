#' @title checkOverlap
#-----------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Reports on how much two spatial objects overlap.
#' @param x A spatial object.
#' @param y A spatial object.
#' @return A two element numeric \emph{vector}.
#' @importFrom raster extent intersect ncell raster
#' @details {Uses \link[raster]{intersect} to estimate the amount of overlap between two spatial objects (i.e. \emph{x}
#' and \emph{y}). The function reports on the percentage of the area in each spatial object represented by their overlap.}
#' @seealso \code{\link{derivePlots}} \code{\link{rankPlots}}
#' @examples {}
#' @export

#-----------------------------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------------------------#

checkOverlap <- function(x, y) {

#-----------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-----------------------------------------------------------------------------------------------------------------------------------------------#

  r1 <- tryCatch(!is.null(extent(x)), error=function(e) return(FALSE)) # is x a spatial object?
  r2 <- tryCatch(!is.null(extent(y)), error=function(e) return(FALSE)) # is y a spatial object?
  if (!r1 & r2) {return(warning('"x" is not a spatial object'))}
  if (!r2 & r1) {return(warning('"y" is not a spatial object'))}
  if (!r1 & !r2) {return(warning('neither "x" and "y" are spatial object'))}

  oc <- tryCatch(!is.null(intersect(x,y)), error=function(e) return(FALSE)) # do they overlap?
  if (!oc) {return(c(0,0))}

  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. Compare input variables
  #-----------------------------------------------------------------------------------------------------------------------------------------------#

  ac <- ncell(raster(extent(intersect(x,y)), res=1)) # how much is the overlap area?
  a1 <- round((ac / ncell(raster(extent(x), res=1))) * 100) # what proportion of x is in the overlap?
  a2 <- round((ac / ncell(raster(extent(y), res=1))) * 100) # what proportion of y is in the overlap?

  return(c(a1,a2))

}
