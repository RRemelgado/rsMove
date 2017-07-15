#' @title backSample
#'
#' @description Background sample selection.
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param rid Vector of region identifyers for each sample.
#' @param method One of \emph{random} or \emph{pca}. Default is \emph{random}.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param nb Number of random background samples.
#' @import raster sp grDevices rgdal
#' @importFrom stats complete.cases prcomp median
#' @return A \emph{SpatialPoints} or a \emph{SpatialPointsDataFrame}.
#' @details {The function selects a set of n background samples where n is determined by 
#' \emph{nb}. If \emph{nb} is not provided all background pixels are considered. If the 
#' method is \emph{random} these samples are returned in the form of a \emph{SpatialPoints} 
#' object. If \emph{pca} is used, a PCA analysis is performed over all available samples 
#' (provided and random). For each unique region ID (\emph{rid}) the median and the Median 
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
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'  
#'  # read movement data
#'  files <- system.file('extdata', 'konstanz_20130805-20130811.shp', package="rsMove")
#'  moveData <- shapefile(files)
#'  
#'  # find sample regions
#'  label <- labelSample(xy=moveData, rad=500, npx=2, pxr=30)
#'  
#'  # read remote sensing data
#'  files <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(files)
#'  
#'  # select background samples
#'  ind <- which(label>0) # selected samples
#'  bSamples <- backSample(xy=moveData[ind,], rid=label[ind], img=rsStk, method='random')
#'  
#' }
#' @export
#' 
#-------------------------------------------------------------------------------------------------------------------------------#

backSample <- function(xy=xy, rid=rid, method=method, img=img, nb=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!exists('rid')) {stop('"rid" is missing')}
  if (length(rid)!=length(xy)) {stop('"xy" and "rid" have different lengths')}
  if (!exists('method')) {method<-'random'}
  if (!method%in%c('random', 'pca')) {stop('"method" is not a valid keyword')}
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}   
  pxr <- res(img)[1] # pixel resolution
  ar <- pxr / 2 # half the resolution
  nr <- dim(img)[1] # numer of rows
  ext0 <- extent(img) # sampling extent
  ext1 <- ext0 # extent used to find positions
  ext1[c(1,3)] <- ext1[c(1,3)] + ar
  ext1[c(2,4)] <- ext1[c(2,4)] - ar
  np <- ncell(img[[1]]) # number of pixels
  op <- crs(img) # output projection
  if (!is.null(nb)) {if (!is.numeric(nb) | length(nb)!=1) {stop('"nb" is not a valid input')}}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 2. extract random background samples
#-------------------------------------------------------------------------------------------------------------------------------#  
  
  # derive coordinates for background samples
  sp <- (round((ext1[4]-xy@coords[,2])/pxr)+1) + nr * round((xy@coords[,1]-ext1[1])/pxr)
  dr <- !duplicated(sp) # find duplicates
  sp <- sp[dr]
  rid <- rid[dr]
  ind <- which(!(1:np)%in%sp) # pixels with no samples
  if (!is.null(nb)) {ind <- ind[sample(1:length(ind), nb, replace=T)]}
  sp <- cbind((ext1[1]+(((sp-1) %/% nr)+1)*pxr), (ext1[4]-(((sp-1) %% nr)+1)*pxr))
  bc <- cbind((ext1[1]+(((ind-1) %/% nr)+1)*pxr), (ext1[4]-(((ind-1) %% nr)+1)*pxr))
  xy <- rbind(sp, bc)
  rid <- c(rid, matrix(0, nrow(bc)))
  
  rm(sp, bc, ind, dr)
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 3. select background samples
#-------------------------------------------------------------------------------------------------------------------------------# 
  
  if (method=='pca') {
    
    # extract environmental information
    nl <- nlayers(img) # number of layers in "img"
    if (nl>1) {edata <- extract(img, xy)}
    if (nl==1) {edata <- extract(img, xy)}
    cc <- complete.cases(edata) # index to remove NA's
    if (nl>1) {edata <- edata[cc,]} else {edata <- edata[cc]}
    rid <- rid[cc]
    
    # kaiser rule
    pcf = function(x) {which((x$sdev^2) > 1)}
    
    # estimate pca and apply kaiser rule
    pca <- prcomp(edata, scale=T, center=T)
    npc <- pcf(pca)
    pca <- data.frame(pca$x[,npc])
      
    # select absences
    uv <- unique(rid[which(rid>0)])
    i0 = which(rid==0)
    ai = vector('list', ncol(pca))
    for (p in 1:length(npc)) {
      usr = vector('list', length(uv))
      for (z in 1:length(uv)) {
        ri = which(rid==uv[z])
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
    ind <- which(rid==0)
    return(SpatialPoints(xy[ind,], proj4string=op))
    
  }

}