#' @title backSample
#'
#' @description Background sample selection.
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param region.id Vector of region identifyers for each sample.
#' @param method One of \emph{random} or \emph{pca}. Default is \emph{random}.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param nr.samples Number of random background samples.
#' @importFrom raster cellFromXY xyFromCell crs ncell
#' @importFrom sp SpatialPoints SpatialPointsDataFrame
#' @importFrom stats complete.cases prcomp median
#' @return A \emph{SpatialPoints} or a \emph{SpatialPointsDataFrame} of background samples for unique pixels in \emph{img}.
#' @details {First, the function determines the unique pixel coordinates for \emph{xy} based on the dimensions
#' of \emph{img} and retrieves n background samples where n is determined by \emph{nr.samples}. Then, the
#' selection of samples is dependent on the method chosen by the user. If \emph{method} is set to \emph{random},
#' the function will select samples randomly. However, if \emph{pca} is used, the function will use a Principal
#' Components Analysis (PCA) over \emph{img} to evaluate the similarity between the samples associated to \emph{xy}
#' and the intial set of random samples First, based on this PCA, the function selects the most important Principal Components
#' (PC's) using the kaiser rule (i.e. PC's with eigenvalues greater than 1). Then, for each PC, the function estimates
#' the median and the Median Absolute Deviation (MAD) for each unique identifier in \emph{region.id}) and selects background
#' samples where the difference between their variance and the variance of the region samples exceeds the MAD. Then, the
#' algorithm removes the background samples that were not selected by all sample regions.}
#'
#' If
#' \emph{nr.samples} is not provided all background pixels are returned
#'
#' object.
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
#'  label <- labelSample(xy=moveData, agg.radius=500, nr.pixels=2, pixel.res=30)
#'
#'  # select background samples
#'  ind <- which(label>0) # selected samples
#'  bSamples <- backSample(xy=moveData[ind,], region.id=label[ind], img=rsStk, method='random')
#'
#' }
#' @export
#'
#-------------------------------------------------------------------------------------------------------------------------------#

backSample <- function(xy=xy, region.id=region.id, method=method, img=img, nr.samples=NULL) {

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
  if (!is.null(nr.samples)) {if (!is.numeric(nr.samples) | length(nr.samples)!=1) {stop('"nr.samples" is not a valid input')}}

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
  if (!is.null(nr.samples)) {ind <- ind[sample(1:length(ind), nr.samples, replace=TRUE)]}
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
