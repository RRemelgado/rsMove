#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param var Object of class RasterStack or RasterBrick containing environmental layers.
#' @import raster grDevices
#' @return Matrix of environmental variables for each sample. If a sample does not overlap NA is returned.
#' @examples \dontrun{
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(x=x, y=y, var=var) {
  
  # check variables
  if (!exists(x)) {stop('error: "x" is missing')}
  if (!exists(y)) {stop('error: "y" is missing')}
  if (!exists(var)) {stop('error: "var" is missing')}
  if (class(var)!='RasterStack' & class(var)!='RasterBrick') {stop('error: "var" is not a valid raster object')}
  if (length(x)!=length(y)) {stop('error: x and y have different lengths')}
  
  # derive pixel coordinates
  ext <- raster::extent(var) # raster extent
  dims <- dim(var) # image dimensions
  pxr <- raster::res(var) # raster resolution
  sp <- (round((ext[4]-y)/pxr)+1) + dims[1] * round((x-ext[1])/pxr) # convert coordinates to pixel positions
  
  # identify pixels with valid/invalid samples
  np <- dims[1] * dims[2] # number of pixels
  if (max(sp) < 0 | min(sp) > np) {stop('error: samples do not overlap with the environmental data (is the projection correct?)')}
  ind0 <- which(sp < 0 | sp > np) # pixels outside the matrix
  ind1 <- which(sp > 0 | sp < np) # pixels inside the valid range
  
  # retrieve/return environmental data
  ov <- matrix(0, length(x), raster::nlayers(var))
  ov[ind0,] <- NA
  ov[ind1,] <- var[sp[ind1]]
  return(ov)
  
}