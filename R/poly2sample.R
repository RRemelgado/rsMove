#' @title poly2sample
#'
#' @description Convert spatial polygons into x and y coordinates.
#' @param pol Object of class "SpatialPolygons" or "SpatialPolygonDataFrame".
#' @param ras Object of class "Raster", "RasterBrick" or "RasterStack".
#' @param mpc Minimum pixel cover (0-100).
#' @param pr Numeric value giving the pixel resolution.
#' @param rp A valid "CRS" object.
#' @param ccheck Logical. Shound the function check for duplicated pixels? Default is FALSE.
#' @import raster grDevices
#' @return Matrix containing unique samples.
#' @details {Determines coordinates of pixels within a given extent or reference raster layer.
#'            "mpc" can be used to filter pixels with low purity, i.e. pixels where the 
#'            percentage of area cover by a polygon is below the defined threshold.
#'            If an extent object is provided, a reference projection system ("rp") 
#'            and a target resoluton ("pr") are required. The output reports on sample 
#'            coordinates ("x" and "y"), polygon id ("id") and pixel percent cover ("cover").
#'            If "ccheck" is TRUE the function checks if pixels were overlapped by multiple polygons.
#'            If so, the code updates the percent cover estimations. This might be slow.}
#' @examples \dontrun{
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

poly2sample(pol=pol, ras=ras, mpc=0, pr=NULL, rp=NULL, ccheck=F) {
  
  # check input
  if(!exists('pol')) {stop('error: "pol" is missing')}
  if(!exists('ras')) {stop('error: "ras" is missing')}
  if(!class(pol)[1]%in%c('SpatialPolygons', 'SpatialPolygonDataFrame')) {
    stop('error: "pol" is not a valid input.')
  }
  if(!class(ras)[1]%in%c("Raster", "RasterStack", "RasterBrick")) {
    if (class(ras)!='Extent') {stop('error: "ras" is not a valid input.')}
    if (is.null(pr)) {stop('error: "ras" is of class "Extent" object. "pr" is required.')}
    if (is.null(rp)) {stop('error: "ras" is of class "Extent" object. "rp" is required.')}
    rr <- raster::raster(ras, res=pr, crs=rp) # create reference raster
  } else {
    rr <- raster::raster(ras) # read reference raster
    pr <- raster::res(ras) # check pixel resolution
    rp <- raster::crs(ras) # check raster projection
    ras <- raster::extent(rr) # extraxt extent
  }
  if (rp@projargs!=raster::crs(pol)@projargs) {stop('error: using different projections')}
  if (!exists('mpc')) {mpc <- 0}
  if (mpc>100 | mpc<0) {stop('error: "mpc" should be between 0 and 100.')}

#----------------------------------------------------------------------------------------------------------#
    
  # check each polygon
  xc <- list()
  yc <- list()
  id <- list()
  pc <- list()
  lf <- function(i) {
    te <- raster::extent(pol[i,]) # target extent
    tr <- raster::crop(rr, te) # target raster
    r = as.matrix(raster::rasterize(pol[i,], tr, getCover=T)) # rasterize polygon
    nr <- dim(r)[1] # number of columns
    ind <- which (r > 0) # identify usable pixels
    xc[[i]] <<- te[1]+((ind / nr)*pr) # x coordinates
    yc[[i]] <<- te[4]-((ind %% nr)*pr) # y coordinates
    id[[i]] <<- replicate(length(ind), i) # polygon ids
    pc[[i]] <<- r[ind] # pixel percent cover
    return()
  }
  x <- sapply(1:np, lf)
  xc <- unlist(xc)
  yc <- unlist(yc)
  id <- unlist(id)
  pc <- unlist(pc)
  
  # convert coordinates to pixel indices
  sp <- (round((ras[2]-yc)/pr)+1) + dim(rr)[1] * round((xc-ras[1])/pr)
  
  rm(x, rr)
  
#----------------------------------------------------------------------------------------------------------#
  
  if (ccheck) {
    
    # remove duplicates (use unique pixel indices as reference)
    # if duplicates exist, sum percent cover over repeated records
    ui <- duplicated(sp) # which are the dupllicates?
    ud <- unique(sp[ui]) # unique duplicates
    ui <- !ui# records to keep
    lf <- function (i) {
      pc[[which(ui==T & sp==i)]] <<- sum(pc[which(sp==i)])
      return()
    }
    x <- sapply(ud, lf)
    xc <- xc[ui]
    yc <- yc[ui]
    sp <- sp[ui]
    id <- id[ui]
    pc <- pc[ui]
    
    rm(x, ud, ui)
    
  }
  
#----------------------------------------------------------------------------------------------------------#
  
  # apply percent cover filter  
  ind <- which(pc > mpc)
  xc <- xc[ind]
  yc <- yc[ind]
  sp <- sp[ind]
  id <- id[ind]
  pc <- pc[ind]
  
  # build/return output
  return(data.frame(x=xc, y=yc, index=sp, id=id, cover=pc))
  
}