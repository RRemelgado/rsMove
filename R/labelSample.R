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
#' @import raster rgdal
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
#'  labels <- labelSample(xy=shortMove, agg.radius=90, nr.pixels=2, pixel.res=30)
#'
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------------------------------#

labelSample <- function(xy=xy, agg.radius=agg.radius, nr.points=NULL, nr.pixels=NULL, pixel.res=pixel.res) {

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # check input variables
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!is.null(nr.pixels) & !is.null(nr.points)) {stop('"nr.pixels" and "nr.points" are both assigned. Choose one')}
  if (!is.null(nr.pixels)) {if (!is.numeric(nr.pixels) | length(nr.pixels)!=1) {stop('"nr.pixels" is not a valid input')}}
  if (!is.null(nr.points)) {if (!is.numeric(nr.points) | length(nr.points)!=1) {stop('"nr.points" is not a valid input')}}
  if (!exists('agg.radius')) {stop('"agg.radius" is missing')}
  if (is.null(pixel.res)) {stop('provide a resolution or a raster')}

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 2. convert samples ot pixel coordinates
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # extract extent of study area
  if (is.numeric(pixel.res)) {
    ext <- extent(xy) # reference extent
    nr <- round((ext[4]-ext[3]) / pixel.res)+1 # number of rows
    nc <- round((ext[2]-ext[1]) / pixel.res)+1 # number of columns
    sp <- (round((ext[4]-xy@coords[,2])/pixel.res)+1) + nr * round((xy@coords[,1]-ext[1])/pixel.res) # convert coordinates to pixel positions
    up <- unique(sp)} # unique pixel positions

  if (class(pixel.res)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    if (crs(xy)@projargs!=crs(pixel.res)@projargs) {stop('"xy" and "pixel.res" have different projections')}
    nr <- dim(pixel.res)[1] # numer of rows
    nc <- dim(pixel.res)[2] # number of columns
    sp <- cellFromXY(pixel.res, xy@coords) # pixel positions
    up <- unique(sp) # unique pixel positions
    pixel.res <- res(pixel.res)[1]} # pixel resolution

  # derive pixel coordinates

  if (length(up)==1) {stop('warning: only one pixel with data found. Processing aborted (is pixel.res correct?)')}

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 3. region label (phase I)
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # filter based on the number of pixels
  if (!is.null(nr.points)) {
    count <- sapply(up, function(x) {sum(sp==x)})
    up <- up[which(count >= nr.points)]}

  # filter samples based on the size of pixel groups
  if (!is.null(nr.pixels)) {

    # evaluate pixel connectivity
    regions <- matrix(0, nr, nc)
    for (r in 1:length(up)) {
      rp <- ((up[r]-1) %% nr)+1
      cp <- ((up[r]-1) %/% nr)+1
      if (cp > 1) {sc<-cp-1} else {sc<-cp}
      if (cp < nc) {ec<-cp+1} else {ec<-cp}
      if (rp > 1) {sr<-rp-1} else {sr<-rp}
      if (rp < nr) {er<-rp+1} else {er<-rp}
      if (max(regions[sr:er,sc:ec])>0) {
        uv <- unique(regions[sr:er,sc:ec])
        uv <- uv[which(uv > 0)]
        mv <- min(uv)
        regions[rp,cp] <- mv
        for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
      } else {regions[rp,cp] <- max(regions)+1}
    }

    # estimate per region pixel count
    uv <- unique(regions[which(regions>0)])
    count <- sapply(uv, function(x) {sum(regions==x)})

    # remove samples related to regions with a pixel count bellow nr.pixels
    uv <- uv[which(count < nr.pixels)]
    for (r in 1:length(uv)) {regions[which(regions==uv[r])]=0}
    up <- up[which(regions[up]>0)]

    rm(regions, count, uv)

  }

  # control sample amountth
  if (length(up)==0) {stop(paste0('there are no regions with >= ', as.character(nr.pixels), ' pixels. Consider reducing "nr.pixels"'))}

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 4. region label (phase II)
#--------------------------------------------------------------------------------------------------------------------------------------------#

  #determine radius (in pixels)
  agg.radius <- round((agg.radius/pixel.res)+0.1)

  # dilate samples
  if (agg.radius > 0) {
    regions <- matrix(0, nr, nc)
    for (p in 1:length(up)) {
      rp <- ((up[p]-1)%%nr) + 1
      cp <- ((up[p]-1)%/%nr) + 1
      if (cp > agg.radius) {sc<-cp-agg.radius} else {sc<-cp}
      if (cp < (nc-agg.radius)) {ec<-cp+agg.radius} else {ec<-cp}
      if (rp > agg.radius) {sr<-rp-agg.radius} else {sr<-rp}
      if (rp < (nr-agg.radius)) {er<-rp+agg.radius} else {er<-rp}
      regions[sr:er,sc:ec]<- 1}

    # dilated sample position
    upd <- which(regions==1)

  } else {upd <- up}

  # evaluate sample connectivity
  regions <- regions * 0
  for (r in 1:length(upd)) {
    rp <- ((upd[r]-1) %% nr)+1
    cp <- ((upd[r]-1) %/% nr)+1
    if (cp > 1) {sc<-cp-1} else {sc<-cp}
    if (cp < nc) {ec<-cp+1} else {ec<-cp}
    if (rp > 1) {sr<-rp-1} else {sr<-rp}
    if (rp < nr) {er<-rp+1} else {er<-rp}
    if (max(regions[sr:er,sc:ec])>0) {
      uv <- unique(regions[sr:er,sc:ec])
      uv <- uv[which(uv > 0)]
      mv <- min(uv)
      regions[rp,cp] <- mv
      for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
    } else {regions[rp,cp] <- max(regions)+1}
  }

#--------------------------------------------------------------------------------------------------------------------------------------------#
# 5. assign region codes and return ID's
#--------------------------------------------------------------------------------------------------------------------------------------------#

  # summarize input variables
  uregions <- regions[up]
  uv <- unique(uregions)
  uv <- uv[which(uv>0)]
  for (r in 1:length(uv)) {uregions[which(regions==uv[r])]<-r}

  rid <- matrix(0,length(sp))
  for (r in 1:length(up)) {rid[which(sp==up[r])] <- uregions[r]}
  rid[which(rid==0)] <- NA
  return(as.vector(rid))

}
