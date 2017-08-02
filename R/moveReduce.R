#' @title moveReduce
#'
#' @description Remote sensing based point segmentation.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param edata Object of class \emph{RasterLayer} or \emph{data.frame}.
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates. 
#' @import raster rgdal
#' @seealso \code{\link{sampleMove}} \code{\link{moveSeg}}
#' @return A \emph{list}.
#' @details {SReduces a set of input samples (\emph{xy}) based on their assignment to unique pixels 
#' within a reference raster (\emph{img}). The function looks at consecutive points ordered by time 
#' (\emph{ot}) and aggregates samples if they remain within the same pixel. If the same pixel is 
#' revisited on a later time, that observation is kept as a separate occurrence. For each temporal 
#' segment, the function returns mean x and y coordinates, the start and end timestamps, the mean 
#' timestamp and the elapsed time.}
#' @examples {
#'  
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'  
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(r))
#'  
#'  # observation time
#'  o.time <- strptime(paste0(moveData@data$date, ' ', moveData@data$time), format="%Y/%m/%d %H:%M:%S")
#'  
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(xy=moveData, ot=o.time, img=r)
#'  
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#  

moveReduce <- function(xy=xy, ot=ot, img=img) {
  
#----------------------------------------------------------------------------------------------------------#
# 1. check input variables
#----------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy) # output projection
  
  # sample dates
  if (!is.null(ot)) {
    if (!class(ot)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"ot" is nof of a valid class')}
    if (length(ot)!=length(xy)) {stop('errorr: "xy" and "ot" have different lengths')}}
  
  # environmental data
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    stop('"img" is not a valid raster object')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "edata" have different projections')}

#----------------------------------------------------------------------------------------------------------#
# 2. identify segments
#----------------------------------------------------------------------------------------------------------#
  
  # convert xy to single pixels
  os <- order(ot)
  xy <- xy[os,]
  ot <- ot[os]
  sp <- cellFromXY(img, xy@coords)
  
  rm(os)
  
  # search for segments and return median values
  sp0 <- 1
  li <- 1
  ux <- list() # x coordinates
  uy <- list() # y coordinates
  ft <- list() # observation time (first)
  ut <- list() # observation time (mean)
  lt <- list() # observation time (last)
  et <- list() # elapsed time
  sg <- list() # segment position
  for (r in 2:length(sp)) {
    if (sp[r]!=sp[r-1]) {
      ep <- (r-1)
      ux[[li]] <- mean(xy@coords[sp0:ep,1])
      uy[[li]] <- mean(xy@coords[sp0:ep,2])
      if (!is.null(ot)) {
        ft[[li]] <- ot[sp0]
        lt[[li]] <- ot[ep]
        ut[[li]] <- mean(ot[sp0:ep])
        et[[li]] <- difftime(ot[ep], ot[sp0], units='mins')
      } else {
        ft[[li]] <- NA
        lt[[li]] <- NA
        ut[[li]] <- NA
        et[[li]] <- NA}
      sg[[li]] <- sp[r-1]
      sp0 <- r
      li <- li + 1
    }
  }
  
  # convert to vector
  ux <- unlist(ux)
  uy <- unlist(uy)
  ft <- do.call("c", ft)
  ut <- do.call("c", ut)
  lt <- do.call("c", lt)
  et <- unlist(et)
  sg <- unlist(sg)
  
  rm(sp, sp0, li)
  
#----------------------------------------------------------------------------------------------------------#
# 3. derive single raster
#----------------------------------------------------------------------------------------------------------#
  
  # estimate time sum per cell
  sp <- unique(sg)
  t.sum <- sapply(sp, function(x) {sum(et[which(sg==x)])})
  
  # build raster
  t.sum.r <- rasterize(xyFromCell(img, sp), crop(img, extent(xy)), t.sum)
  
#----------------------------------------------------------------------------------------------------------#
# 4. build output
#----------------------------------------------------------------------------------------------------------#
  
  df <- data.frame(x=ux, y=uy, timestamp=ut, start.time=ft, end.time=lt, elapsed.time=et)
  r.shp <- SpatialPointsDataFrame(df[,1:2], df, proj4string=crs(xy))
  
  return(list(points=r.shp, total.time=t.sum.r))
  
}