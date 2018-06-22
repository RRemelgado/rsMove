#' @title poly2sample
#'
#' @description {Converts a raster grid to points depending on how much each pixel is covered by a polygon.}
#' @param x Object of class \emph{SpatialPolygons} or \emph{SpatialPolygonDataFrame}.
#' @param y A raster object or a numeric element.
#' @param min.cover Minimum percent a pixel should be covered by a polygon for sampling (1-100). Default is 1.
#' @importFrom raster raster extent crop rasterToPoints rasterize xyFromCell cellFromXY crs
#' @importFrom sp SpatialPointsDataFrame
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{SpatialPointsDataFrame} with sampled pixels reporting on polygon percent coverage.
#' @details {\emph{poly2Sample} extends on the \code{\link[raster]{rasterize}} function from the raster package making
#' it more efficient over large areas and converting its output into point samples rather than a raster object. For each
#' polygon in (\emph{"x"}), \emph{poly2sample} extracts the overlapping pixels of a reference grid. Then, for each pixel,
#' the function estimates the percentage of it that is covered by the reference polygon. Finally, the function extracts
#' coordinate pairs for pixels that has a percent coverage equal to or greater than \emph{min.cover}. If \emph{y} is a
#' raster object, the function will use it as a reference  grid. If \emph{y} is a numeric element, the function will build
#' a raster from the extent of \emph{x} and a resolution equal to \emph{y}.}
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
#'  samples <- poly2sample(roi, img, min.cover=50)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------#

poly2sample <- function(x, y, min.cover=1) {

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

  # check cover value
  if (is.null(min.cover)) {min.cover <- 100}
  if (min.cover < 1 | min.cover > 100) {stop('"min.cover" should be between 0 and 100')}

#-------------------------------------------------------------------------------------------------------------------------#
# 2. evaluate polygons
#-------------------------------------------------------------------------------------------------------------------------#

  lf <- function(i) {
    r <- extend(crop(y, extent(x[i,])), res(y)*2)
    r <- rasterToPoints(rasterize(x[i,], r, getCover=TRUE))
    ind <- which(r[,3] > min.cover) # usable pixels
    if (length(ind) > 0) {return(data.frame(x=r[ind,1], y=r[ind,2], c=r[ind,3]))} else {return(NULL)}}

  df0 <- lapply(1:length(x), lf)

#-------------------------------------------------------------------------------------------------------------------------#
# 3. remove duplicated pixels and update pecent count
#-------------------------------------------------------------------------------------------------------------------------#

  i <- sapply(df0, function(r) {!is.null(r)})

  if (sum(i) > 0) {

    df0 <- do.call(rbind, df0[i])
    pp <- cellFromXY(y, df0[,c("x", "y")]) # pixel positions
    up <- unique(pp) # unique positions
    pc <- sapply(up, function(x) {sum(df0$c[which(pp==x)])}) # update percentages
    pc[which(pc>100)] <- 100 # account for miss-calculations (related to e.g. overlapping polygons)
    xy <- xyFromCell(y, up) # convert unique positions to coordinates
    df0 <- data.frame(x=xy[,1], y=xy[,2], cover=pc) # build final data frame

    rm(pp, up, pc, xy)

    # return output
    if (is.null(nrow(df0))) {
      warning('no pixels with a percent cover greater than 0')
      return(NULL)} else {

        return(SpatialPointsDataFrame(df0[,1:2], df0, proj4string=crs(x)))

      }

  } else {

    warning('no pixels with a percent cover greater than 0')
    return(NULL)

  }

#------------------------------------------------------------------------------------------------------------------------#

}

