#' @title poly2sample
#'
#' @description {Converts a raster grid to points.}
#' @param x Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonDataFrame}.
#' @param y A raster object or a numeric element.
#' @importFrom raster raster extent crop rasterToPoints rasterize xyFromCell cellFromXY crs
#' @importFrom sp SpatialPointsDataFrame
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{SpatialPointsDataFrame} with sampled pixels reporting on polygon percent coverage.
#' @details {\emph{poly2Sample} extends on the \code{\link[raster]{rasterize}} function from the raster package making
#' it more efficient over large areas and converting its output into point samples rather than a raster object. For each
#' polygon in (\emph{"x"}), \emph{poly2sample} extracts the overlapping pixels of a reference grid. Then, for each pixel,
#' the function estimates the percentage of it that is covered by the reference polygon. If \emph{y} is a raster object, 
#' the function will use it as a reference  grid. If \emph{y} is a numeric element, the function will build a raster from 
#' the extent of \emph{x} and a resolution equal to \emph{y}.}
#' @examples {
#'
#'  require(raster)
#'
#'  # load example probability image
#'  file <- system.file('extdata', 'probabilities.tif', package="rsMove")
#'  img <- raster(file)
#'
#'  # load area of interest
#'  file <- system.file('extdata', 'roi.shp', package="rsMove")
#'  roi <- shapefile(file)
#'
#'  # extract samples
#'  samples <- poly2sample(roi, img)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------#

poly2sample <- function(x, y) {

#-------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#-------------------------------------------------------------------------------------------------------------------------#

  # check shapefile
  if(is.null('x')) {stop('"x" is missing')}
  if(!class(x)[1]%in%c('SpatialPolygons', 'SpatialPolygonsDataFrame')) {
    stop('"x" is not a valid input')}

  # check/derive reference raster
  if (is.numeric(y)) {y <- extend(raster(extent(x), res=y, crs=crs(x)), y)} else {
    e <- try(extent(y))
    if (class(e) == "try-error") {stop('"y" is not of a valid class')}}

  # check overlap between x and y
  o <- checkOverlap(x, y)
  if (o[1] != 100) {stop('"x" is not contained by "y"')}

#-------------------------------------------------------------------------------------------------------------------------#
# 2. evaluate polygons
#-------------------------------------------------------------------------------------------------------------------------#

  lf <- function(i) {
    r <- extend(crop(y[[1]], extent(x[i,])), 1)
    r <- rasterToPoints(rasterize(x[i,], r, getCover=TRUE))
    return(data.frame(x=r[,1], y=r[,2], c=r[,3], id=i))}

  df0 <- lapply(1:length(x), lf)

#-------------------------------------------------------------------------------------------------------------------------#
# 3. remove duplicated pixels and update pecent count
#-------------------------------------------------------------------------------------------------------------------------#

  df0 <- do.call(rbind, df0)
  df0 <- df0[which(df0$c > 0),]
  
  if (nrow(df0) > 0) {
    
    pp <- cellFromXY(y, df0[,c("x", "y")]) # pixel positions
    up <- unique(pp) # unique positions
    pc <- sapply(up, function(x) {sum(df0$c[which(pp==x)])}) # update percentages
    pc[which(pc>100)] <- 100 # account for miss-calculations (related to e.g. overlapping polygons)
    xy <- xyFromCell(y, up) # convert unique positions to coordinates
    df0 <- data.frame(x=xy[,1], y=xy[,2], cover=pc) # build final data frame
    
    rm(pp, up, pc, xy)
    
    # return output
    return(SpatialPointsDataFrame(df0[,1:2], df0, proj4string=crs(x)))
    
    
  } else {
    
    warning('no samples with a percent coverage >= 1%')
    return(NULL)
    
  }

#------------------------------------------------------------------------------------------------------------------------#

}
