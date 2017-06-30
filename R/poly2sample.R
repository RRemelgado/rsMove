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

poly2sample <- function(pol=pol, re=NULL, mpc=0, pr=NULL) {
  
  #check if package is isntalled
  pkgTest <- function(x){
    if (!require(x,character.only = TRUE)){
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
  pkgTest("raster")
  
#-------------------------------------------------------------------------------------------------------------------------#
  
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
      rr <- raster(re, crs=rp, res=pr) # build reference raster
      nr <- round((re[4]-re[3]) / pr) + 1 # number of rows in raster
    } else {
      pr <- res(re)[1] # check pixel resolution
      rp <- crs(re) # check raster projection
      nr <- dim(re)[1] # number of rows in raster
      rr <- re # set variable with reference raster
      re <- extent(rr) # extraxt extent
    }
    if (rp@projargs!=crs(pol)@projargs) {stop('error: using different projections')}
  } else {
    if (is.null(pr)) {stop('error: provide "pr" since "re" is missing')}
    re <- extent(pol) # extent derived from 
    rr <- raster(re, crs=crs(pol), res=pr)
    nr <- round((re[4]-re[3]) / pr) + 1 # number of rows in raster
    
  }
  if (!exists('mpc')) {mpc <- 0}
  if (mpc>100 | mpc<0) {stop('error: "mpc" should be between 0 and 100.')}

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
  
  return(data.frame(x=xc, y=yc, index=up, cover=pc))
  
#------------------------------------------------------------------------------------------------------------------------#
  
}