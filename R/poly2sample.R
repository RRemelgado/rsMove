#' @title poly2sample
#'
#' @description Convert spatial polygons into point samples.
#' @param pol Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonDataFrame}.
#' @param re Object of class \emph{Extent} or a raster object from which an extent can be derived.
#' @param mpc Minimum pixel cover (0-100). Default is 100.
#' @param pr Pixel resolution.
#' @import sp raster rgdal
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{SpatialPointsDataFrame}.
#' @details {Determines coordinates of pixels within a given extent. If \emph{re} is missing the
#'            target extent corresponds to the extent of the polygon layer. In this case, \emph{pr}
#'            is required. If \emph{re} is an \emph{Extent} object or \emph{re} is missing \emph{rp} is required. 
#'            \emph{mpc} can be used to filter pixels with low purity, i.e. pixels where the 
#'            percentage of area cover by a polygon is below the defined threshold. 
#'            The output provides the coordinates (\emph{x} and \emph{y}) the 
#'            pixel percent cover (\emph{cover}) for each selected pixel.}
#' @examples {
#'  
#'  require(raster)
#'  
#'  # load example probability image
#'  file <- system.file('extdata', 'konstanz_probabilities.tif', package="rsMove")
#'  img <- raster(file)
#'  
#'  # load area of interest
#'  file <- system.file('extdata', 'konstanz_roi.shp', package="rsMove")
#'  roi <- shapefile(file)
#'  
#'  # segment probabilities
#'  samples <- poly2sample(pol=roi, re=img)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------#

poly2sample <- function(pol=pol, re=NULL, mpc=NULL, pr=NULL) {

#-------------------------------------------------------------------------------------------------------------------------#
  
  # check input
  if(is.null('pol')) {stop('"pol" is missing')}
  if(!class(pol)[1]%in%c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
    stop('"pol" is not a valid input.')
  }
  if (!is.null(re)) {
    if(!class(re)[1]%in%c("RasterLayer", "RasterStack", "RasterBrick")) {
      if (class(re)!='Extent') {stop('"ras" is not a valid input.')}
      if (is.null(pr)) {stop('"re" is of class "Extent" object. "pr" is required.')}
      if (is.null(rp)) {stop('"re" is of class "Extent" object. "rp" is required.')}
      rr <- raster(re, crs=rp, res=pr) # build reference raster
      nr <- round((re[4]-re[3]) / pr) + 1 # number of rows in raster
    } else {
      pr <- res(re)[1] # check pixel resolution
      rp <- crs(re) # check raster projection
      nr <- dim(re)[1] # number of rows in raster
      rr <- re # set variable with reference raster
      re <- extent(rr) # extraxt extent
    }
    if (rp@projargs!=crs(pol)@projargs) {stop('using different projections')}
  } else {
    if (is.null(pr)) {stop('provide "pr" since "re" is missing')}
    re <- extent(pol) # extent derived from 
    rr <- raster(re, crs=crs(pol), res=pr)
    nr <- round((re[4]-re[3]) / pr) + 1 # number of rows in raster
  }
  if (is.null(mpc)) {mpc <- 100}
  if (mpc < 0 | mpc > 100) {stop('"mpc" should be between 0 and 100')}

  # update extent object
  # (correspond x/y to pixel center)
  ar <- pr / 2 # half resolution
  re[1] <- re[1] + ar
  re[2] <- re[2] - ar
  re[3] <- re[3] + ar
  re[4] <- re[4] - ar
  
#-------------------------------------------------------------------------------------------------------------------------#
    
  # function to apply
  lf <- function(i) {
    r <- crop(rr, extent(pol[i,]))
    r <- rasterToPoints(rasterize(pol[i,], r, getCover=T))
    ind <- which(r[,3] > 0) # usable pixels
    pp <- (round((re[4]-r[ind,2])/pr)) + nr * (round((r[ind,1]-re[1])/pr))
    pc <- r[ind,3] # percent cover
    return(list(pp=pp, pc=pc))
  }
  np <- length(pol) # number of polygons
  x <- lapply(1:np, lf) # evaluate polygons
  pp <- unlist(sapply(x, function(x) {x$pp})) # positions
  pc <- unlist(sapply(x, function(x) {x$pc})) # percent cover
  
  rm(x)
  
#-------------------------------------------------------------------------------------------------------------------------#
  
  # remove duplicated pixels and update pecent count
  
  up <- unique(pp) # unique ppositions
  pc <- sapply(up, function(x) {sum(pc[which(pp==x)])}) # update percentages
  pc[which(pc>100)] <- 100 # account for miss-calculations (related to e.g. overlapping polygons)
  
  rm(pp)
  
#-------------------------------------------------------------------------------------------------------------------------#
  
  # apply percent cover filter  
  ind <- which(pc >= mpc)
  up <- up[ind]
  pc <- pc[ind]
  
#------------------------------------------------------------------------------------------------------------------------#
  
  # build/return output
  xc <- re[1] + (round(up/nr)*pr) # convert positions to x coordinates
  yc <- re[4] - (round(up%%nr)*pr) # convert positions to y coordinates
  
  os <- data.frame(x=xc, y=yc, index=up, cover=pc)
  return(SpatialPointsDataFrame(cbind(xc,yc), os, proj4string=crs(pol)))
  
#------------------------------------------------------------------------------------------------------------------------#
  
}