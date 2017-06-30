#' @title dataQuery
#'
#' @description Query environmental data for coordinate pairs.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param var Object of class "RasterStack" or "RasterBrick".
#' @param bs Buffer size (unit depends on the raster projection).
#' @param remove.duplicates Logical. Default if FALSE.
#' @param fun Passes an external function.
#' @import raster grDevices
#' @return Matrix of environmental variables for each sample.
#' @details {Returns environmental variables from a raster object for a given set of x and y coordinates.
#'          If "remove.duplicates" is TRUE, the function checks for and removes duplicated samples that 
#'          fall in the same pixels. The pixel ID, contained in the output, informs on which samples were kept.
#'          If a sample does not overlap with the reference raster NA is returned.
#'          If "bs" is set, the function uses a buffer arround each sample and extracts a summary metric.
#'          If "fun" is NULL, the function uses a weighted mean by default.}
#' @examples \dontrun{
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

dataQuery <- function(x=x, y=y, var=var, bs=NULL, remove.duplicates=F, fun=NULL) {
  
  # check variables
  if (!exists('x')) {stop('error: "x" is missing')}
  if (!exists('y')) {stop('error: "y" is missing')}
  if (!exists('var')) {stop('error: "var" is missing')}
  if (class(var)!='RasterStack' & class(var)!='RasterBrick') {stop('error: "var" is not a valid raster object')}
  if (length(x)!=length(y)) {stop('error: x and y have different lengths')}
  if (!is.null(bs)) {if (!is.numeric(bs)) {stop('error: "bs" assigned but not numeric')}}
  if (!is.null(fun)) {fun <- function(x) {sum(x*x) / sum(x)}} else {
    if (!is.function(fun)) {stop('error: "fun" is not a valid function')}
    tf <- fun(c(1,2,3)) # test function
    if (length(tf)!=1) {stop('error: provided function returns more than 1 output.')}
  }
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # derive pixel coordinates
  ext <- raster::extent(var) # raster extent
  dims <- dim(var) # image dimensions
  pxr <- raster::res(var)[1] # raster resolution
  sp <- (round((ext[4]-y)/pxr)+1) + dims[1] * round((x-ext[1])/pxr) # convert coordinates to pixel positions
  if (remove.duplicates) {
    sp <- unique(sp) # return unique records
    id <- which(!duplicate(sp)) # return indices of unique records
  } else {
    id <- 1:length(sp)  
  }

#-------------------------------------------------------------------------------------------------------------------------------#
  
  # identify pixels with valid/invalid samples
  np <- dims[1] * dims[2] # number of pixels
  if (max(sp) < 0 | min(sp) > np) {stop('error: samples do not overlap with the environmental data (is the projection correct?)')}
  ind0 <- which(sp < 0 | sp > np) # pixels outside the matrix
  ind1 <- which(sp > 0 | sp < np) # pixels inside the valid range

#-------------------------------------------------------------------------------------------------------------------------------#
  
  # extract environmental data
  ov <- matrix(0, length(sp), (raster::nlayers(var)+1))
  ov[ind0,] <- NA
  ov[,1] <- id
  if (!is.null(bb)) {
    bs <- round(bs / pxr)
    bf <- function(x) {
      xp <- x / dims[1]
      yp <- x %% dims[1]
      sr <- yp-bs
      if (sr < 1) {sr<-1}
      er <- yp+bs
      if (er > dims[1]) {er=dims[1]}
      sc <- xp-bs
      if (sc < 1) {sc<-1}
      ec <- xp+bs
      if (ec > dims[2]) {ec<-dims[2]}
      return(fun(var[sc:ec,sr:er]))}
    ov[ind1,2:ncol(ov)] <- unlist(sapply(ind1, bf))
  } else {ov[ind1,2:ncol(ov)] <- var[sp[ind1]]}
  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # return output
  return(ov)
  
}