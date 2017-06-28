#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param var Object of class RasterStack or RasterBrick containing environmental layers.
#' @param remove.duplicates Optional. If TRUE, checks for and removes duplicated samples that fall in the same pixels.
#' @import raster grDevices
#' @return Matrix of environmental variables for each sample. If a sample does not overlap NA is returned.
#' @examples \dontrun{
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(x=x, y=y, var=var, remove.duplicates=F) {
  
  # check variables
  if (!exists('x')) {stop('error: "x" is missing')}
  if (!exists('y')) {stop('error: "y" is missing')}
  if (!exists('var')) {stop('error: "var" is missing')}
  if (class(var)!='RasterStack' & class(var)!='RasterBrick') {stop('error: "var" is not a valid raster object')}
  if (length(x)!=length(y)) {stop('error: x and y have different lengths')}
  
  # derive pixel coordinates
  ext <- raster::extent(var) # raster extent
  dims <- dim(var) # image dimensions
  pxr <- raster::res(var) # raster resolution
  sp <- (round((ext[4]-y)/pxr)+1) + dims[1] * round((x-ext[1])/pxr) # convert coordinates to pixel positions
  if (remove.duplicates) {
    sp <- unique(sp) # return unique records
    id <- which(!duplicate(sp)) # return indices of unique records
  } else {
    id <- 1:length(sp)  
  }
  
  # identify pixels with valid/invalid samples
  np <- dims[1] * dims[2] # number of pixels
  if (max(sp) < 0 | min(sp) > np) {stop('error: samples do not overlap with the environmental data (is the projection correct?)')}
  ind0 <- which(sp < 0 | sp > np) # pixels outside the matrix
  ind1 <- which(sp > 0 | sp < np) # pixels inside the valid range
  
  # retrieve/return environmental data
  ov <- matrix(0, length(sp), (raster::nlayers(var)+1))
  ov[ind0,] <- NA
  ov[,1] <- id
  ov[ind1,2:ncol(ov)] <- var[sp[ind1]]
  return(ov)
  
}