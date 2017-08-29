#' @title backSample
#'
#' @description Background sample selection.
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param region.id Vector of region identifyers for each sample.
#' @param method One of \emph{random} or \emph{pca}. Default is \emph{random}.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param n.samples Number of random background samples.
#' @import raster sp rgdal
#' @importFrom stats complete.cases prcomp median
#' @return A \emph{SpatialPoints} or a \emph{SpatialPointsDataFrame}.
#' @details {The function selects a set of n background samples where n is determined by
#' \emph{n.samples}. If \emph{n.samples} is not provided all background pixels are considered. If the
#' method is \emph{random} these samples are returned in the form of a \emph{SpatialPoints}
#' object. If \emph{pca} is used, a PCA analysis is performed over all available samples
#' (provided and random). For each unique region ID (\emph{region.id}) the median and the Median
#' Absolute Deviation (MAD) is estimated for a given Principal Component (PC). Background
#' samples where the difference between their variance and the variance of the region samples
#' exceeds the MAD are selected. Only the samples that are selected by all unique regions are
#' kept. This process is repeated for all PC's with eigenvalues greater than 1 (Kaiser rule).
#' The final set of samples corresponds to the unique records selected by the different PC's.
#' These samples are returned as a SpatialPointsDataFrame containing the variables extracted
#' from \emph{img}.}
#' @seealso \code{\link{labelSample}} \code{\link{hotMove}} \code{\link{dataQuery}}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130805-20130811.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(rsStk))
#'
#'  # find sample regions
#'  label <- labelSample(xy=moveData, rad=500, npx=2, pxr=30)
#'
#'  # select background samples
#'  ind <- which(label>0) # selected samples
#'  bSamples <- backSample(xy=moveData[ind,], region.id=label[ind], img=rsStk, method='random')
#'
#' }
#' @export
#'
#-------------------------------------------------------------------------------------------------------------------------------#

backSample <- function(xy=xy, region.id=region.id, method=method, img=img, n.samples=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-------------------------------------------------------------------------------------------------------------------------------#

  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!exists('region.id')) {stop('"region.id" is missing')}
  if (length(region.id)!=length(xy)) {stop('"xy" and "region.id" have different lengths')}
  if (!exists('method')) {method<-'random'}
  if (!method%in%c('random', 'pca')) {stop('"method" is not a valid keyword')}
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}
  np <- ncell(img[[1]]) # number of pixels
  op <- crs(img) # output projection
  if (!is.null(n.samples)) {if (!is.numeric(n.samples) | length(n.samples)!=1) {stop('"n.samples" is not a valid input')}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. extract random background samples
#-------------------------------------------------------------------------------------------------------------------------------#

  # convert presences to pixel positions
  sp <- cellFromXY(img[[1]], xy@coords)

  # remove duplicated records
  dr <- !duplicated(sp)
  sp <- sp[dr]
  region.id <- region.id[dr]

  # derice background samples
  ind <- which(!(1:np)%in%sp)
  if (!is.null(n.samples)) {ind <- ind[sample(1:length(ind), n.samples, replace=TRUE)]}
  xy <- rbind(xyFromCell(img[[1]], sp), xyFromCell(img[[1]], ind))
  region.id <- c(region.id, replicate(length(ind), 0))

  rm(sp, ind)

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. select background samples
#-------------------------------------------------------------------------------------------------------------------------------#

  if (method=='pca') {

    # extract environmental information
    edata <- as.data.frame(extract(img, xy))
    cc <- complete.cases(edata) # index to remove NA's
    edata <- edata[cc,]
    region.id <- region.id[cc]
    xy <- xy[cc,]

    # kaiser rule
    pcf = function(x) {which((x$sdev^2) > 1)}

    # estimate pca and apply kaiser rule
    pca <- prcomp(edata, scale=TRUE, center=TRUE)
    npc <- pcf(pca)
    pca <- data.frame(pca$x[,npc])

    # select absences
    uv <- unique(region.id[which(region.id>0)])
    i0 = which(region.id==0)
    ai = vector('list', ncol(pca))
    for (p in 1:length(npc)) {
      usr = vector('list', length(uv))
      for (z in 1:length(uv)) {
        ri = which(region.id==uv[z])
        s1 = median(pca[ri,p])
        s2 = median(abs(pca[ri,p]-s1))
        usr[[z]] = i0[which(abs(pca[i0,p]-s1) > s2)]
      }
      usr = unlist(usr)
      ui = unique(usr)
      count = vector('numeric', length(ui))
      for (z in 1:length(ui)) {count[z] = length(which(usr==ui[z]))}
      ai[[p]] = ui[which(count==length(uv))]
    }

    # return samples
    ai <- unique(unlist(ai))
    return(SpatialPointsDataFrame(xy[ai,], as.data.frame(edata[ai,]), proj4string=op))

  }

  if (method=='random') {

    # return samples
    ind <- which(region.id==0)
    return(SpatialPoints(xy[ind,], proj4string=op))

  }

}
