#' @title labelSample
#'
#' @description {Pixel-based labeling of spatially connected groups of samples for splitting them between training and validation.}
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param agg.radius Minimum radius for pixel aggregation. Unit depends on the projection of the data.
#' @param nr.points Minimum number of samples per pixel.
#' @param nr.pixels Minimum number of pixels per region.
#' @param pixel.res Pixel resolution or a valid raster layer.
#' @references \href{10.1002/rse2.70}{Remelgado, R., Leutner, B., Safi, K., Sonnenschein, R., Kuebert, C. and Wegmann, M. (2017), Linking animal movement and remote sensing - mapping resource suitability from a remote sensing perspective. Remote Sens Ecol Conserv.}
#' @return A \emph{vector} of unique identifiers assigning each point in \emph{xy} to their correspondent pixel region. Filtered observations are returned as \emph{NA}.
#' @details {First, the samples are converted to pixel coordinates and removes pixels with a corresponding number of points greater
#' than \emph{nr.points}. Then, if \emph{nr.pixels} is set, the connectivity between neighboring samples is evaluated. Internally, the
#' function will label groups of pixels based on their connectivity and regions with a pixel count smaller than the one specified by
#' \emph{nr.pixels} are excluded. Then, the algorithm aggregates nearby regions using a dilation algorithm within the radius specified
#' by \emph{agg.radius} and proceeds to reliable the pixels covered by samples. Finally, this information is used to label the original
#' samples provided by \emph{xy} based on their corresponding pixel coordinates. This analysis is based on the spatial extent of \emph{xy}
#' and a given pixel resolution (\emph{pixel.res}). Alternatively, the user may assign a raster object to \emph{pixel.res}.}
#' @importFrom raster crs cellFromXY extent raster res freq clump rowFromCell colFromCell
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

labelSample <- function(xy, pixel.res, agg.radius=NULL, nr.points=NULL, nr.pixels=NULL) {

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # check input variables
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!is.null(nr.pixels)) {if (!is.numeric(nr.pixels) | length(nr.pixels)!=1) {stop('"nr.pixels" is not a valid input')}}
  if (!is.null(nr.points)) {if (!is.numeric(nr.points) | length(nr.points)!=1) {stop('"nr.points" is not a valid input')}}
  if (!is.null(agg.radius)) {if (!is.numeric(agg.radius) | length(agg.radius)!=1) {stop('"agg.radius" is not a valid input')}}

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 2. convert samples ot pixel coordinates
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # build sample mask (from extent)
  if (is.numeric(pixel.res)) {
    ext <- extent(xy)
    pixel.res <- raster(ext, res=pixel.res, crs=crs(xy), vals=NA)
    sp <- cellFromXY(pixel.res, xy)
    up <- unique(sp)
    up <- up[!is.na(up)]
    pixel.res[up] <- 1}

  # build sample mask (from raster)
  if (class(pixel.res)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    if (crs(xy)@projargs!=crs(pixel.res)@projargs) {stop('"xy" and "pixel.res" have different projections')}
    sp <- cellFromXY(pixel.res, xy@coords)
    up <- unique(sp)
    up <- up[!is.na(up)]
    pixel.res <- raster(extent(pixel.res), res=res(pixel.res), crs=crs(pixel.res), vals=NA)
    pixel.res[up] <- 1}

  if (length(up)==1) {stop('warning: only one pixel in the data. Processing aborted (is pixel.res correct?)')}

  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # 3. region label (phase I)
  #--------------------------------------------------------------------------------------------------------------------------------------------#

  # filter based on the number of pixels
  if (!is.null(nr.points)) {
    count <- sapply(up, function(x) {sum(sp==x)})
    up <- up[which(count >= nr.points)]
    rm(count)}

  # filter samples based on the size of pixel groups
  if (!is.null(nr.pixels)) {
    regions <- clump(pixel.res)
    px.freq <- freq(regions,useNA="no")
    ind <- which(px.freq[,2] < nr.pixels)
    ind <- as.vector(px.freq[ind,1])
    if (!is.na(ind[1])) {for (i in 1:length(ind)) {pixel.res[regions==ind[i]] <- NA}}
    rm(regions, px.freq, ind)
    up <- up[pixel.res[up] > 0]}

  # control sample amountth
  if (length(up)==0) {stop(paste0('there are no regions with >= ', as.character(nr.pixels), ' pixels. Consider reducing "nr.pixels"'))}

  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # 4. region label (phase II)
  #--------------------------------------------------------------------------------------------------------------------------------------------#

  #determine radius (in pixels)
  agg.radius <- round((agg.radius/res(pixel.res)[1])+0.1)

  # dilate samples
  if (agg.radius > 0) {
    pixel.res <- focal(pixel.res, matrix(0,agg.radius, agg.radius), function(j) {sum(!is.na(j))}, pad=TRUE, padValue=NA) > 0
    pixel.res[pixel.res == 0] <- NA}

  # evaluate sample connectivity
  pixel.res <- clump(pixel.res)

  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # 5. return ID's
  #--------------------------------------------------------------------------------------------------------------------------------------------#

  return(pixel.res[sp])

}
