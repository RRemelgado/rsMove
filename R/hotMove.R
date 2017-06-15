#' @title hotMove
#'
#' @description Identifies hotspots of samples and labels these accordingly. It build a matrix for a given spatial extent and resolution converting these samples into unique pixel coordinates. Then, the spatial connectivity of these pixels is evaluated using a 8-neighboor connected component labelling algorithm.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param res Resolution of the grid used for region labelling. Unir depends on the projection of the input.
#' @param ext Optional. Four element vector with a spatial extent expressed as c(xmin, xmax, ymin, ymax). If none is provided the extent is derived from the x and y coordinates.
#' @param shp Logical. If TRUE the function returns a shapefile defined by the convex hull of the samples within each region. Default is FALSE.
#' @param proj If shp=TRUE then a valid CRS indicating the projection of the input coordinates is required.
#' @return Returns a vector containing numeric labels that represent unique sample regions. If 'shp' is in use polygons are provided in addition.
#' @seealso \code{\link{extent}}, \code{\link{sampleMove}}, \code{\link{crs}}
#' @examples \dontrun{
#' hotMove(x=whiteStork$lon, y=whiteStork=lon, res=0.01, ext)
#' }

hotMove <- function(x=x, y=y, res=res, ext=NULL, shp=F, cs=NULL) {

  #---------------------------------------------------------------------------------------------------------------------#
  #  1. check inpur variables
  #---------------------------------------------------------------------------------------------------------------------#

  if (!exists(x)) {stop('error: missing x coordinates')}
  if (!exists(y)) {stop('error: missing y coordinates')}
  if (length(x)!=length(y)) {stop('error: x and y dimensions do not match')}
  if (!exists(res)) {stop('error: a value is required for res')}
  if (is.null(ext)) {ext <- c(min(x), max(x), min(y), max(y))}
  if (shp==T) {
    if (is.null(cs)) {stop('error: cs is missing (provide a valid crs object')}
    if (class(cs)!='CRS') {stop('error: cs not a valid crs object')}
  }

  #---------------------------------------------------------------------------------------------------------------------#
  # 2. determine grid coordinates for given pixels
  #---------------------------------------------------------------------------------------------------------------------#

  # derive pixel coordinates
  nc <- round((ext[2]-ext[1]) / res) + 1 # number of columns
  if (nc < 0) {stop('error: number of columns negative (is x min/max correct in ext?)')}
  nr <- round((ext[4]-ext[3]) / res) + 1 # number of rows
  if (nr < 0) {stop('error: number of rows negative (is y min/max correct in ext?)')}
  sp <- (round((ext[4]-y)/res)+1) + nr * round((x-ext[1])/res) # convert coordinates to pixel positions
  up <- unique(sp) # unique pixel positions
  if (length(up)==1) {stop('warning: only one pixel with data found. Processing aborted (is res correct?)')}

  #---------------------------------------------------------------------------------------------------------------------#
  # 3. find unique sample regions
  #---------------------------------------------------------------------------------------------------------------------#

  # evaluate pixel connectivity
  regions <- matrix(0, nc, nr)
  for (r in 1:length(up)) {
    xp <- up %% nr
    yp <- up / nr
    if (xp > 1) {sc<-xp-1} else {sc<-xp}
    if (xp < nc) {ec<-xp+1} else {ec<-xp}
    if (yp > 1) {sr<-yp-1} else {sr<-yp}
    if (yp < nr) {er<-yr+1} else {er<-yp}
    if (max(regions[sc:ec,sr:er])>0) {
      mv <- min(regions[sc:ec,sr:er])
      uv = unique((regions[sc:ec,sr:er])[which(regions[sc:ec,sr:er]>0)])
      for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
    } else {regions[xpos,ypos] <- max(regions)+1}
  }

  # update region id
  uv = unique(regions[which(regions>0)])
  uregions = regions
  for (u in 1:length(uv)) {uregions[which(regions==uv[u])] <- u}

  rm(regions)

  #---------------------------------------------------------------------------------------------------------------------#
  # 4. produce classified samples
  #---------------------------------------------------------------------------------------------------------------------#

  # label smaples based on the value associated to each unique position
  r <- matrix(0, length(x))
  for (u in 1:length(up)) {
    loc = where(sp==up[u])
    r[loc] = uregions[up[u]]
  }
  rm(uregions, loc)

  if (shp==T) {

    # function to convert samples to polygons
    pc <- function(x) {
      ind <- which(r==x) # target coordinates
      xy <- cbind(x[ind],y[ind])
      ch <- grDevices::chull(xy) # locate polygon vertices
      return(sp::Polygon(xy[c(ch,ch[1]),])) # build closed polygon
    }

    # merge polygons
    uv = unique(r)
    shp <- sp::SpatialPolygons(lapply(uv, pc))
    sp::crs(shp) <- cs

    # return output
    return(list(indices=r, polygons=shp))

  } else {list(indices=r)} #  if shp=F return only the sample indices

}
