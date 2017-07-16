#' @title hotMove
#'
#' @description Identifies hotspots of samples using a pixel based approach.
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param pxr Grid resolution. Unit depends on \emph{xy} projection.
#' @param shp Logical. Should the function provide polygons? Default is FALSE.
#' @return A list object.
#' @import sp raster rgdal
#' @details {The function builds a matrix for a given resolution and a 
#' spatial extent derived from a set of samples (\emph{xy}). First, these 
#' samples are converted into unique pixel coordinates. Then, the 
#' spatial connectivity of these pixels is evaluated using a 8-neighboor 
#' connected component labelling algorithm. The function returns a vector 
#' with the region ID per sample (\emph{$indices}). If \emph{shp} is TRUE the 
#' function also returns a shapefile for each region defined by the 
#' convex hull of the samples within them (\emph{$polygons}).}
#' @seealso \code{\link{sampleMove}} \code{\link{hotMoveStats}}
#' @examples {
#' 
#' require(rgdal)
#' require(raster)
#' require(sp)
#' 
#' # reference data
#' file <- system.file('extdata', 'latLon_example.shp', package="rsMove")
#' moveData <- shapefile(file)
#' 
#' # extract regions
#' hm <- hotMove(xy=moveData, pxr=0.1, shp=TRUE)
#' 
#' # plot shapefile (color by region)
#' plot(hm$polygons, col=hm$indices)
#' 
#' # add new information to original shapefile
#' moveData@data <- cbind(moveData@data, hm$indices)
#' 
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

hotMove <- function(xy=xy, pxr=pxr, shp=FALSE) {
  
#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!exists('pxr')) {stop('"pxr" is missing')}
  if (!is.logical(shp)) {stop('"shp" is not a valid logical argument')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#
  
  # derive pixel coordinates
  ext <- extent(xy)
  nc <- round((ext[2]-ext[1]) / pxr) + 1 # number of columns
  if (nc < 0) {stop('number of columns negative (is x min/max correct in ext?)')}
  nr <- round((ext[4]-ext[3]) / pxr) + 1 # number of rows
  if (nr < 0) {stop('number of rows negative (is y min/max correct in ext?)')}
  sp <- (round((ext[4]-xy@coords[,2])/pxr)+1) + nr * round((xy@coords[,1]-ext[1])/pxr) # convert coordinates to pixel positions
  up <- unique(sp) # unique pixel positions
  if (length(up)==1) {stop('only one pixel with data. Processing aborted (is "pxr" correct?)')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 3. find unique sample regions
#---------------------------------------------------------------------------------------------------------------------#
  
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
  
  # update region id
  uv = unique(regions[which(regions>0)])
  uregions = regions
  for (u in 1:length(uv)) {uregions[which(regions==uv[u])] <- u}
  
#   # convert region layer to raster
#   ar <- (pxr/2)
#   or = raster(uregions)
#   extent(or) <- c(ext[1]-ar, ext[2]+ar, ext[3]-ar, ext[4]+ar) # map extent
#   res(or) <- pxr
#   crs(or) <- crs(xy)
  
#---------------------------------------------------------------------------------------------------------------------#
  # 4. produce classified samples
  #---------------------------------------------------------------------------------------------------------------------#

  # label smaples based on the value associated to each unique position
  rid <- vector('numeric', length(xy))
  for (r in 1:length(up)) {rid[which(sp==up[r])]<-uregions[up[r]]}
  
  # dummy buffer (used if a region has only 1 point)
  d.buff <- pxr*0.001
  
  # convert samples to polygons
  if (shp==T) {
    
    # uniqu values
    uv <- sort(unique(rid))
    
    # function to build polygons
    pc <- function(x) {
      ind <- which(rid==x) # target coordinates
      if (length(ind==1)) {
        ic <- xy@coords[ind,]
        ic <- rbind((ic+c(d.buff, 0)), (ic+c(d.buff, -d.buff)), (ic+c(0, d.buff)))
      } else {ic <- xy@coords[ind,]}
      ch <- chull(ic) # determine polygon edge coordinates
      return(Polygons(list(Polygon(ic[c(ch,ch[1]),])), ID=x)) # build closed polygon
    }
    
    # build/merge polygons
    uv <- sort(unique(rid))
    shp <- SpatialPolygons(lapply(uv, pc), proj4string=crs(xy))
    
    # return output
    return(list(indices=rid, polygons=shp))

  } else {return(list(indices=rid))} #  if shp=F return only the sample indices

}