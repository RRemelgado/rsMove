#' @title hotMove
#'
#' @description Detection of geographic regions of samples using a pixel based approach.
#' @param xy Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param pixel.res Grid resolution. Unit depends on \emph{xy} projection.
#' @param return.shp Logical. Should the function return polygons of the regions? Default is FALSE.
#' @return List containing a vector of region ID's per data point (\emph{$indices}) and region polygons (\emph{$polygons}).
#' @importFrom sp Polygon Polygons SpatialPolygons
#' @importFrom raster extent crs
#' @importFrom grDevices chull
#' @details {First, the function builds a raster with a resolution equal to \emph{pixel.res} and the
#' spatial extent of \emph{xy}. Then, each point in \emph{xy} is converted into pixel coordinates. Based
#' on the unique pixel coordinates, the function then evaluates the spatial connectivity of these pixels
#' using a 8-neighbor connected component labeling algorithm to detect regions. Finally, the ID's are related
#' back to each individual data point in \emph{xy} based on their pixel coordinates and - if \emph{return.shp} is TRUE -
#'  a polygon is derived from the convex hull of the points within each region.}
#' @seealso \code{\link{sampleMove}} \code{\link{hotMoveStats}}
#' @examples {
#'
#' require(raster)
#'
#' # reference data
#' data(longMove)
#'
#' # extract regions
#' hm <- hotMove(longMove, 0.1, return.shp=TRUE)
#'
#' # plot shapefile (color by region)
#' plot(hm$polygons, col=hm$indices)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

hotMove <- function(xy, pixel.res, return.shp=FALSE) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.null(crs(xy)@projargs)) {stop('"xy" is missing a valid projection')}
  if (!exists('pixel.res')) {stop('"pixel.res" is missing')}
  if (!is.logical(return.shp)) {stop('"return.shp" is not a valid logical argument')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#

  regions <- raster(extent(xy), res=pixel.res, crs=crs(xy), vals=0) # build raster from extent
  sp <- cellFromXY(regions, xy) # convert coordinates to pixel positions
  up <- unique(sp) # unique pixel positions
  if (length(up)==1) {stop('only one pixel with data. Processing aborted (is "pixel.res" correct?)')}

#---------------------------------------------------------------------------------------------------------------------#
# 3. find unique sample regions
#---------------------------------------------------------------------------------------------------------------------#

  regions[up] <- 1 # assign value to overlapping pixels
  regions <- clump(regions) # determine region connectivity

#---------------------------------------------------------------------------------------------------------------------#
  # 4. produce classified samples
  #---------------------------------------------------------------------------------------------------------------------#

  # label smaples based on the value associated to each unique position
  rid <- vector('numeric', length(xy))
  for (r in 1:length(up)) {rid[which(sp==up[r])]<-uregions[up[r]]}

  # dummy buffer (used if a region has only 1 point)
  d.buff <- pixel.res*0.001

  # convert samples to polygons
  if (return.shp==TRUE) {

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
