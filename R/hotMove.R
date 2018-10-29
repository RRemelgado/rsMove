#' @title hotMove
#'
#' @description Detection of geographic regions of samples using a pixel based approach.
#' @param x Object of class \emph{SpatialPoints} of \emph{SpatialPointsDataFrame}.
#' @param y Grid resolution. Unit depends on \emph{x} projection.
#' @param return.shp Logical. Should the function return polygons of the regions? Default is FALSE.
#' @return A List containing a vector of region ID's per data point (\emph{$indices}) and region polygons (\emph{$polygons}).
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPolygonsDataFrame
#' @importFrom ggplot2 fortify ggplot aes_string geom_polygon labs
#' @importFrom raster extent crs
#' @importFrom grDevices chull
#' @importFrom plyr join
#' @details {First, the function builds a raster with a resolution equal to \emph{y} and the
#' spatial extent of \emph{x}. Then, each point in \emph{x} is converted into pixel coordinates.
#' Based on the unique pixel coordinates, the function then evaluates the spatial connectivity of
#' these pixels using a 8-neighbor connected component labeling algorithm to detect regions. Finally,
#' the ID's are related back to each individual data point in \emph{x} based on their pixel coordinates
#' and - if \emph{return.shp} is TRUE - a polygon is derived from the convex hull of the points within
#' each region. The output of the function consists of:
#'  \itemize{
#'  \item{\emph{region.id} - Vector reporting on the region each element in \emph{x} belongs to.}
#'  \item{\emph{polygons} - Polygons for each unique region in \emph{region.id} specified by the data column \emph{region}.}
#'  \item{\emph{plot} - A plot with the output of \emph{polygons} accessible if \emph{return.shp} is TRUE.}
#'  }
#'  }
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
#' plot(hm$polygons, col=hm$region.id)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

hotMove <- function(x, y, return.shp=FALSE) {

  #---------------------------------------------------------------------------------------------------------------------#
  #  1. check inpur variables
  #---------------------------------------------------------------------------------------------------------------------#

  if (!exists('x')) {stop('"x" is missing')}
  if (!class(x)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"x" is not of a valid class')}
  if (is.null(crs(x)@projargs)) {stop('"x" is missing a valid projection')}
  if (!exists('y')) {stop('"y" is missing')}
  if (!is.logical(return.shp)) {stop('"return.shp" is not a valid logical argument')}

  #---------------------------------------------------------------------------------------------------------------------#
  # 2. determine grid coordinates for given pixels
  #---------------------------------------------------------------------------------------------------------------------#

  regions <- raster(extend(extent(x), c(y,y)), res=y, crs=crs(x), vals=0) # build raster from extent
  sp <- cellFromXY(regions, x) # convert coordinates to pixel positions
  up <- unique(sp) # unique pixel positions
  if (length(up)==1) {stop('only one pixel with data. Processing aborted (is "y" correct?)')}

  #---------------------------------------------------------------------------------------------------------------------#
  # 3. find unique sample regions
  #---------------------------------------------------------------------------------------------------------------------#

  regions[up] <- 1 # assign value to overlapping pixels
  regions <- clump(regions) # determine region connectivity

  #---------------------------------------------------------------------------------------------------------------------#
  # 4. produce classified samples
  #---------------------------------------------------------------------------------------------------------------------#

  # label smaples based on the value associated to each unique position
  rid <- vector('numeric', length(x))
  for (r in 1:length(up)) {rid[which(sp==up[r])] <- regions[up[r]]}

  # dummy buffer (used if a region has only 1 point)
  d.buff <- y*0.001

  # convert samples to polygons
  if (return.shp==TRUE) {

    # uniqu values
    uv <- sort(unique(rid))

    # function to build polygons
    pc <- function(i) {
      ind <- which(rid==i) # target coordinates
      if (length(ind==1)) {
        ic <- x@coords[ind,]
        ic <- rbind((ic+c(d.buff, 0)), (ic+c(d.buff, -d.buff)), (ic+c(0, d.buff)))
      } else {ic <- x@coords[ind,]}
      ch <- chull(ic) # determine polygon edge coordinates
      return(Polygons(list(Polygon(ic[c(ch,ch[1]),])), ID=i)) # build closed polygon
    }

    # build/merge polygons
    uv <- sort(unique(rid))
    tmp <- lapply(uv, pc)
    shp <- SpatialPolygonsDataFrame(SpatialPolygons(tmp, proj4string=crs(x)), data.frame(id=uv))
    shp.df <- fortify(shp, region="id")
    shp.df <- join(shp.df, shp@data, by="id")
    shp.df$id <- factor(shp.df$id, levels=as.character(sort(as.numeric(unique(shp.df$id)))))
    shp.p <- ggplot(shp.df, aes_string(x="long", y="lat", group="group", fill="id")) +
      geom_polygon() + labs(x="Lon", y="Lat", fill="Region ID")

    # return output
    return(list(region.id=rid, polygons=shp, plot=shp.p))

  } else {return(list(region.id=rid))} #  if shp=F return only the sample indices

}
