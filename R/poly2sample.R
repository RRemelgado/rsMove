#' @title poly2sample
#'
#' @description Convert spatial polygons into x and y coordinates.
#' @param pol Object of class "SpatialPolygons" or "SpatialPolygonDataFrame".
#' @param ras Object of class "Raster", "RasterBrick" or "RasterStack".
#' @param mpc Minimum pixel cover (0-100).
#' @param pr Numeric value giving the pixel resolution.
#' @param rp A valid "CRS" object.
#' @import raster grDevices
#' @return Matrix containing unique samples.
#' @details {Determines coordinates of pixels within a given extent or reference raster layer.
#'            "mpc" can be used to filter pixels with low purity, i.e. pixels where the 
#'            percentage of area cover by a polygon is below the defined threshold.
#'            If an extent object is provided, a reference projection system ("rp") 
#'            and a target resoluton ("pr") are required. The output reports on sample 
#'            coordinates ("x" and "y"), polygon id ("id") and pixel percent cover ("cover").}
#' @examples \dontrun{
#' }
#' @export

#--------------------------------------------------------------------------------------------------------------------------------------------------#

poly2sample(pol=pol, ras=ras, mpc=0, pr=NULL, rp=NULL) {
  
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
  }
  if (rp@projargs!=raster::crs(pol)@projargs) {stop('error: using different projections')}
  if (!exists('mpc')) {mpc <- 0}
  if (mpc>100 | mpc<0) {stop('error: "mpc" should be between 0 and 100.')}
  
  # output variables
  xc <- list() # x coordinates
  yc <- list() # y coordinates
  id <- list() # polygon id coordinates
  pc <- list() # pixel cover
  
  # check each polygon
  lp <- 1
  for (p in 1:length(pol)) {
    te <- raster::extent(pol[p,]) # target extent
    tr <- crop(rr, te) # target raster
    r = as.matrix(raster::rasterize(pol[p,], tr, getCover=T)) # rasterize polygon
    nr <- dim(r)[1] # number of columns
    ind <- which (r > mpc) # identify usable pixels
    if (length(ind)>0) {
      xc[[lp]] <- te[1]+((ind / nr)*pr) # convert pixel coordinates to x
      yc[[lp]] <- te[4]-((ind %% nr)*pr) # convert pixel coordinates to y
      id[[lp]] <- replicate(length(ind), p) # assign polygon ID ot coordinates
      pc[[lp]] <- r[ind] # assign pixel cover information
      lp <- lp + 1
    }
    rm(te, tr, r, ind)
  }
  
  # build /return data frame
  xc <- unlist(xc)
  yc <- unlist(yc)
  id <- unlist(id)
  pc <- unlist(pc)
  df <- data.frame(x=xc, y=yc, id=id, cover=pc, stringsAsFactors=F)
  return(df)
  
}