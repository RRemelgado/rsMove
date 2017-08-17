#' @title labelSample
#'
#' @description Region labeling of samples based on their spatial connectivity.
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param rad Minimum radius. Unit depends on the projection of the data.
#' @param npt Minimum pixel count per pixel.
#' @param npx Minimum number of pixels.
#' @param pxr Pixel resolution os a valid raster layer.
#' @return A \emph{vector}.
#' @details {First, the samples are converted to pixel coordinates and then one  
#' of two occur: 1) if \emph{npt} is set, the function removes pixels with a pixel 
#' count smaller than the one specified; 2) If \emph{npx} is set, the connectivity 
#' between neighboring samples is evaluated and regions with a count small than the 
#' specified value are filtered. Only one option may be set at a time. Then, the 
#' remaining pixels are dilated and the samples are again labeled accounting for 
#' regions that are not connected but are near to each other. Regions within a given 
#' distance of each other (defined by \emph{rad}) are aggregated. The grid used for 
#' this analysis is built from the spatial extent of \emph{xy} and a given pixel 
#' resolution (\emph{pxr}). If \emph{pxr} is a raster, this will be used to define 
#' the dimensions of the grid. Doing so can be of use when the user has pre-select 
#' environmental predictors that will be used for modeling. Note that the finer the 
#' resolution the more independent regions are likely to be returned. The output is 
#' a vector with ID's assigning each sample to its region. Samples filtered by 
#' \emph{npt} or \emph{npx} will be returned as zeros.}
#' @import raster rgdal
#' @seealso \code{\link{sampleMove}} \code{\link{hotMove}}
#' @examples {
#'  
#'  require(raster)
#'  
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'  
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(r))
#'  
#'  # derive region labels
#'  labels <- labelSample(xy=moveData, rad=90, npx=2, pxr=30)
#'  
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------------------------------#

labelSample <- function(xy=xy, rad=rad, npt=NULL, npx=NULL, pxr=r) {
  
#--------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------------------------------------------------------------------#
  
  # check input variables
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!is.null(npx) & !is.null(npt)) {stop('"npx" and "npt" are both assigned. Choose one')}
  if (!is.null(npx)) {if (!is.numeric(npx) | length(npx)!=1) {stop('"npx" is not a valid input')}}
  if (!is.null(npt)) {if (!is.numeric(npt) | length(npt)!=1) {stop('"npt" is not a valid input')}}
  if (!exists('rad')) {stop('"rad" is missing')}
  if (is.null(pxr)) {stop('provide a resolution or a raster')}
  
#--------------------------------------------------------------------------------------------------------------------------------------------#
# 2. convert samples ot pixel coordinates  
#--------------------------------------------------------------------------------------------------------------------------------------------#
  
  # extract extent of study area
  if (is.numeric(pxr)) {
    ext <- extent(xy) # reference extent
    nr <- round((ext[4]-ext[3]) / pxr)+1 # number of rows
    nc <- round((ext[2]-ext[1]) / pxr)+1 # number of columns
    sp <- (round((ext[4]-xy@coords[,2])/pxr)+1) + nr * round((xy@coords[,1]-ext[1])/pxr) # convert coordinates to pixel positions
    up <- unique(sp)} # unique pixel positions
  
  if (class(pxr)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    if (crs(xy)@projargs!=crs(pxr)@projargs) {stop('"xy" and "pxr" have different projections')}
    nr <- dim(pxr)[1] # numer of rows
    nc <- dim(pxr)[2] # number of columns
    sp <- cellFromXY(pxr, xy@coords) # pixel positions
    up <- unique(sp) # unique pixel positions
    pxr <- res(pxr)[1]} # pixel resolution

  # derive pixel coordinates
  
  if (length(up)==1) {stop('warning: only one pixel with data found. Processing aborted (is pxr correct?)')}
  
#--------------------------------------------------------------------------------------------------------------------------------------------#
# 3. region label (phase I)
#--------------------------------------------------------------------------------------------------------------------------------------------#
  
  # filter based on the number of pixels
  if (!is.null(npt)) {
    count <- sapply(up, function(x) {sum(sp==x)})
    up <- up[which(count >= npt)]}
  
  # filter samples based on the size of pixel groups
  if (!is.null(npx)) {
    
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
    
    # remove samples related to regions with a pixel count bellow npx
    uv <- uv[which(count < npx)]
    for (r in 1:length(uv)) {regions[which(regions==uv[r])]=0}
    up <- up[which(regions[up]>0)]
    
    rm(regions, count, uv)
    
  }
  
  # control sample amountth
  if (length(up)==0) {stop(paste0('there are no regions with >= ', as.character(npx), ' pixels. Consider reducing "npx"'))}
  
#--------------------------------------------------------------------------------------------------------------------------------------------#
# 4. region label (phase II)
#--------------------------------------------------------------------------------------------------------------------------------------------#
  
  #determine radius (in pixels)
  rad <- round((rad/pxr)+0.1)
  
  # dilate samples
  regions <- matrix(0, nr, nc)
  for (p in 1:length(up)) {
    rp <- ((up[p]-1)%%nr) + 1
    cp <- ((up[p]-1)%/%nr) + 1
    if (cp > rad) {sc<-cp-rad} else {sc<-cp}
    if (cp < (nc-rad)) {ec<-cp+rad} else {ec<-cp}
    if (rp > rad) {sr<-rp-rad} else {sr<-rp}
    if (rp < (nr-rad)) {er<-rp+rad} else {er<-rp}
    regions[sr:er,sc:ec]<- 1}
  
  # dilated sample position
  upd <- which(regions==1)
  
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
  return(rid)
  
}