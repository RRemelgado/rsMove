#' @title labelSample
#'
#' @description {Pixel-based labeling of spatially connected groups of points in a \emph{SpatialPoints} object.}
#' @param x Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param y Pixel resolution or a valid raster layer.
#' @param agg.radius Minimum radius for pixel aggregation. Unit depends on the projection of the data.
#' @param nr.points Minimum number of observations per pixel.
#' @param nr.pixels Minimum number of pixels per region.
#' @return A numeric \emph{vector} with region identifiers for each observation in \emph{x} to their correspondent pixel region. Filtered observations are returned as \emph{NA}.
#' @details {First, the observations are converted to pixel coordinates and pixels with a corresponding number of observations greater than \emph{nr.points}
#' are filtered out. Then, if \emph{nr.pixels} is set, the function evaluates the spatial connectivity of the pixels and regions with a pixel count smaller
#' than \emph{nr.pixels} are filtered out. Then, the algorithm aggregates nearby regions within the distance specified by \emph{agg.radius}. The final region
#' identifiers are then assigned back to the original observations in \emph{x} based on their corresponding pixel coordinates. This analysis is based on the spatial
#' extent of \emph{x} and a given pixel resolution (\emph{y}). Alternatively, the user may assign a raster object as \emph{y} assuring that the
#' final output is aligned with it.}
#' @importFrom raster crs cellFromXY extent raster res freq clump rowFromCell colFromCell focal
#' @seealso \code{\link{sampleMove}} \code{\link{hotMove}}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # derive region labels
#'  labels <- labelSample(shortMove, r, agg.radius=60)
#'
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------------------------------#

labelSample <- function(x, y, nr.points=1, nr.pixels=NULL, agg.radius=NULL) {

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # check input variables
  if (!class(x)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"x" is not of a valid class')}
  if (is.null(crs(x)@projargs)) {stop('"x" is missing a valid projection')}
  if (!is.null(nr.pixels)) {if (!is.numeric(nr.pixels) | length(nr.pixels)!=1) {stop('"nr.pixels" is not a valid input')}}
  if (!is.null(nr.points)) {if (!is.numeric(nr.points) | length(nr.points)!=1) {stop('"nr.points" is not a valid input')}}
  if (!is.null(agg.radius)) {if (!is.numeric(agg.radius) | length(agg.radius)!=1) {stop('"agg.radius" is not a valid input')}}

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 2. convert samples ot pixel coordinates
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # build sample mask (from extent)
  if (is.numeric(y)) {
    ext <- extent(x)
    y <- raster(ext, res=y, crs=crs(x), vals=NA)
    sp <- cellFromXY(y, x)
    up <- unique(sp)
    up <- up[!is.na(up)]
    y[up] <- 1}

  # build sample mask (from raster)
  if (class(y)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "y" have different projections')}
    sp <- cellFromXY(y, x@coords)
    up <- unique(sp)
    up <- up[!is.na(up)]
    y <- raster(extent(y), res=res(y), crs=crs(y), vals=NA)
    y[up] <- 1}

  if (length(up)==1) {stop('warning: only one pixel in the data. Processing aborted (is y correct?)')}

  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # 3. region label (phase I)
  #--------------------------------------------------------------------------------------------------------------------------------------------#

  # filter based on the number of pixels
  if (!is.null(nr.points)) {
    count <- sapply(up, function(i) {sum(sp==i)})
    up <- up[which(count >= nr.points)]
    rm(count)}

  # filter samples based on the size of pixel groups
  if (!is.null(nr.pixels)) {
    regions <- clump(y)
    px.freq <- freq(regions,useNA="no")
    ind <- which(px.freq[,2] < nr.pixels)
    ind <- as.vector(px.freq[ind,1])
    if (!is.na(ind[1])) {for (i in 1:length(ind)) {y[regions==ind[i]] <- NA}}
    rm(regions, px.freq, ind)
    up <- up[y[up] > 0]}

  # control sample amountth
  if (length(up)==0) {stop(paste0('there are no regions with >= ', as.character(nr.pixels), ' pixels. Consider reducing "nr.pixels"'))}

  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # 4. region label (phase II)
  #--------------------------------------------------------------------------------------------------------------------------------------------#

  #determine radius (in pixels)
  agg.radius <- (round((agg.radius/res(y)[1])+0.1)*2) + 1

  # dilate samples
  if (agg.radius > 0) {
    y <- focal(y, matrix(0,agg.radius, agg.radius), function(j) {sum(!is.na(j))}, pad=TRUE, padValue=NA) > 0
    y[y == 0] <- NA}

  # evaluate sample connectivity
  y <- clump(y)

  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # 5. return ID's
  #--------------------------------------------------------------------------------------------------------------------------------------------#

  return(y[sp])

}
