#' @title poly2sample
#'
#' @description Convert spatial polygons into x and y coordinates.
#' @param pol Object of class "SpatialPolygons" or "SpatialPolygonDataFrame".
#' @param ras Object of class "Raster".
#' @param mpc Minimum pixel cover. Percent proportion a pixel should be covered by a polygon so that that sample is kept. Default is 1%.
#' @import raster grDevices
#' @return Matrix containing unique samples. The output reports on sample coordinates ("x" and "y"), polygon id ("id") and pixel % cover ("cover").
#' @examples \dontrun{
#' }
#' @export

#--------------------------------------------------------------------------------------------------------------------------------------------------#

poly2sample(pol, ras, mpc=1) {
  
  if(!class(pol)[1]%in%c('SpatialPolygons', 'SpatialPolygonDataFrame')) {}
  
  
  # reference raster specifications
  rr <- raster::raster(ras) # raster
  pr <- raster::res(ras) # pixel resolution
  re <- raster:extent(ras)[c(1,3)] # extent
  nr <- dim(ras)[1] # number of columns
  np <- length(pol) # number of polygons
  
  # output variables
  xc <- vector('list', np) # x coordinates
  yc <- vector('list', np) # y coordinates
  id <- vector('list', np) # polygon id coordinates
  pc <- vector('list', np) # pixel cover
  
  # check each polygon
  for (p in 1:np) {
    te <- raster::extent(pol[1,]) # target extent
    tr <- crop(rr, te) # target raster
    r = as.matrix(raster::rasterize(pol[1,], tr, getCover=T)) # rasterize polygon
    ind <- which (r > mpc) # identify usable pixels
    xc[[p]] <- ((ind / nr)*pr) - re[1] # convert pixel coordinates to x
    yc[[p]] <- re[2] - ((ind %% nr)*pr) # convert pixel coordinates to y
    id[[p]] <- replicate(length(ind), p) # assign polygon ID ot coordinates
    pc[[p]] <- r[ind] # assign pixel cover information
  }
  
  # build /return
  df <- data.frame(x=unlist(xc), y=unlist(yc), id=unlist(id), cover=pc, stringsAsFactors=F)
  return(df)
  
}