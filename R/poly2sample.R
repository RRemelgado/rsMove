#' @title poly2sample
#'
#' @description Convert spatial polygons into point samples.
#' @param pol Object of class "SpatialPolygons" or "SpatialPolygonDataFrame".
#' @param re Object of class "Extent" or a raster object from which an extent can be derived.
#' @param mpc Minimum pixel cover (0-100).
#' @param pr Pixel resolution.
#' @param rp A valid "CRS" object.
#' @import raster grDevices
#' @return Matrix containing unique samples.
#' @details {Determines coordinates of pixels within a given extent. If "re" is missing the
#'            target extent corresponds to the extent of the polygon layer. In this case, "pr"
#'            is required. If "re" is an "Extent" object or "re" is missing rp" is required. 
#'            "mpc" can be used to filter pixels with low purity, i.e. pixels where the 
#'            percentage of area cover by a polygon is below the defined threshold. 
#'            The output provides the coordinates ("x" and "y") the 
#'            pixel percent cover ("cover") for each selected pixel.}
#' @examples \dontrun{
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------#

poly2sample(pol=pol, re=NULL, mpc=0, pr=NULL) {
  
  # check input
  if(is.null('pol')) {stop('error: "pol" is missing')}
  if(!class(pol)[1]%in%c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
    stop('error: "pol" is not a valid input.')
  }
  if (!is.null(re)) {
    if(!class(re)[1]%in%c("RasterLayer", "RasterStack", "RasterBrick")) {
      if (class(re)!='Extent') {stop('error: "ras" is not a valid input.')}
      if (is.null(pr)) {stop('error: "re" is of class "Extent" object. "pr" is required.')}
      if (is.null(rp)) {stop('error: "re" is of class "Extent" object. "rp" is required.')}
      rr <- round((re[4]-re[3]) / pr) + 1 # number of rows in raster
    } else {
      pr <- raster::res(re)[1] # check pixel resolution
      rp <- raster::crs(re) # check raster projection
      rr <- dim(re)[1] # number of rows in raster
      re <- raster::extent(re) # extraxt extent
    }
    if (rp@projargs!=raster::crs(pol)@projargs) {stop('error: using different projections')}
  } else {
    if (is.null(pr)) {stop('error: provide "pr" since "re" is missing')}
    re <- raster::extent(pol) # extent derived from 
    rr <- round((re[4]-re[3]) / pr) + 1 # number of rows in raster
  }
  if (!exists('mpc')) {mpc <- 0}
  if (mpc>100 | mpc<0) {stop('error: "mpc" should be between 0 and 100.')}

#-------------------------------------------------------------------------------------------------------------------------#
    
  # function to apply
  lf <- function(i) {
    te <- raster::extent(pol[i,]) # target extent
    nc <- (round((te[2]-re[1])/pr)+1)-(round((te[1]-re[1])/pr)+1) # number of columns
    nr <- (round((re[4]-te[3])/pr)+1)-(round((re[4]-te[4])/pr)+1) # number of rows
    r <-raster(nrows=nr, ncols=nc, xmn=te[1], xmx=te[2], ymn=te[3], ymx=te[4], crs=rp, resolution=pr) # target raster
    r = as.matrix(raster::rasterize(pol[i,], r, getCover=T)) # rasterize polygon
    nr <- dim(r)[1] # number of columns
    ind <- which (r > 0) # identify usable pixels
    xp <- round(((te[1]+((ind/nr)*pr))-ras[1])/pr)+1 # x image coordinates
    yp <- round((ras[4]-(te[4]-((ind%%nr)*pr)))/pr)+1 # y image coordinates
    return(list(pp=(yp+rr*xp), pc=r[ind])) # return metrics
  }
  np <- length(pol) # number of polygons
  x <- lapply(1:np, lf) # evaluate polygons
  sp <- unlist(sapply(x, function(x) {x$pp})) # pixel positions in raster
  pc <- unlist(sapply(x, function(x) {x$pc})) # pixel percent cover
  
  rm(x)
  
#-------------------------------------------------------------------------------------------------------------------------#
  
  # remove duplicated pixels and update pecent count
  up <- unique(sp) # unique pixels
  pc <- sapply(up, function(x) {sum(pc[which(sp==x)])}) # update percentages
  pc[which(pc>100)] <- 100 # account for miss-calculations (related to e.g. overlapping polygons)
  
  rm(sp)
  
#-------------------------------------------------------------------------------------------------------------------------#
  
  # apply percent cover filter  
  ind <- which(pc >= mpc)
  up <- up[ind]
  pc <- pc[ind]
  
#------------------------------------------------------------------------------------------------------------------------#
  
  # build/return output
  xc <- re[1] + ((up/rr)*pr) # convert positions to x coordinates
  yc <- re[4] - ((up%%rr)*pr) # convert positions to y coordinates
  
  return(data.frame(x=xc, y=yc, index=up, cover=pc))
  
#------------------------------------------------------------------------------------------------------------------------#
  
}