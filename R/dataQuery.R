#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param var Object of class RasterStack or RasterBrick containing environmental layers.
#' @import raster, grDevices
#' @return Matrix of environmental variables for each sample. If a sample does not overlap NA is returned.
#' @examples \dontrun{
#' }

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(x=x, y=y, var=var) {
  
  # check variables
  if (!exists(xy)) {stop('error: "xy" is missing')}
  if (!exists(var)) {stop('error: "var" is missing')}
  if (class(var)!='raster') {stop('error: "var" is not a valid raster object')}
  if (nlength(x)!=length(y)) {stop('error: x and y have different lengths')}
  
  # derive pixel coordinates
  ext <- raster::extent(var) # raster extent
  pxr <- raster::res(var) # raster resolution
  nc <- round((ext[2]-ext[1]) / pxr) + 1 # number of columns
  nr <- round((ext[4]-ext[3]) / pxr) + 1 # number of rows
  sp <- (round((ext[4]-y)/pxr)+1) + nr * round((x-ext[1])/pxr) # convert coordinates to pixel positions
  
  # identify pixels with valid/invalid samples
  np <- nr*nc # number of pixels
  if (max(sp) < 0 | min(sp) > np) {stop('error: samples do not overlap with the environmental data (is the projection correct?)')}
  ind0 <- which(sp < 0 | sp > np) # pixels outside the matrix
  ind1 <- which(sp > 0 | sp < np) # pixels inside the valid range
  
  # retrieve/return environmental data
  ov <- matrix(0, length(x), nlayers(var))
  ov[ind0,] <- NA
  ov[ind1,] <- var[sp[ind1]]
  return(ov)
  
}