#' @title moveReduce
#'
#' @description Pixel based summary of movement data that preserves periodic movements.
#' @param x Object of class \emph{spatVector}.
#' @param y Object of class \emph{spatRaster}.
#' @param z Object of class \emph{Date} or \emph{POSIXct} with the observation time of each element in \emph{x}.
#' @param preserve.revisits Logical. Should the function preserve revisit patterns?
#' @importFrom terra rast crs cellFromXY rasterize geom vect
#' @importFrom plyr ddply summarise
#' @return A \emph{spatVector} object.
#' @details {Translates (\emph{x}) into pixel coordinates within a reference
#' raster (\emph{y}). The function identifies temporal segments corresponding
#' to groups of consecutive observations within the same pixel. In this process,
#' revisits to recorded pixels are preserved. Once the segments are identified,
#' the function derives mean x and y coordinates for each of them and evaluates
#' the time spent within each pixel. The function reports on the start and end
#' timestamps and the elapsed time. If \emph{preserve.revisits} is FALSE, the
#' function will then summarize the output on a pixel level summing the time
#' spent at each pixel. Additionally, if \emph{derive.raster} is TRUE, the
#' function will derive a \emph{RasterLayer} with the same configuration as
#' \emph{y} depicting the the total amount of time spent per pixel.
#' The output of the function consists of a \emph{spatVector}
#' with the reduced sample set.}
#' @examples {
#'
#'  require(terra)
#'
#'  # read raster data
#'  r <- (system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  # read movement data
#'  shortMove <- read.csv(system.file('extdata', 'shortMove.csv', package="rsMove"))
#'
#'  # convert observations to vector
#'  shortMove = vect(shortMove, geom=c("x","y"), crs="EPSG:32632")
#'
#'  # observation time
#'  z <- as.POSIXct(strptime(paste0(shortMove$date, ' ', shortMove$time),
#'  format="%Y/%m/%d %H:%M:%S"))
#'
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(shortMove, r, z, derive.raster=TRUE)
#'
#' }
#' @export

#-----------------------------------------------------------------------------#

moveReduce <- function(x, y, z, preserve.revisits=TRUE) {

  #---------------------------------------------------------------------------#
  # 1. check input variables
  #---------------------------------------------------------------------------#

  # samples
  if (!class(x)[1]%in%c('spatVector')) {stop('"x" is not of a valid class')}

  # sample dates
  if (!class(z)[1]%in%c('Date', 'POSIXct')) {stop('"z" is nof of a valid class')}
  if (length(z)!=length(x)) {stop('"x" and "z" have different lengths')}
  if (sum(is.na(z)) > 0) {stop('please filter missing values in "z"')}

  # environmental data
  if (!class(y)[1]%in%c('spatRaster')) stop('"y" is not a valid raster object')
  if (crs(x)==crs(y)) {stop('"x" and "edata" have different projections')}

  if (!is.logical(preserve.revisits)) stop('"preserve.revisits" not logical')

  #---------------------------------------------------------------------------#
  # 2. identify segments for each temporal window
  #---------------------------------------------------------------------------#

  # combine inputs
  samples = data.frame(geom(x)[,c("x","y")])
  samples$date = as.POSIXct(z)

  rm(x, z)

  # convert x to single pixels
  samples = samples[order(samples$date),]
  samples$cell.id <- cellFromXY(y, samples[,c("x","y")])

  # search for segments and return sample indices
  pd = rle(samples$cell.id)$lengths
  samples$seg.id = vector('numeric', nrow(samples))
  for (p in 1:length(pd)) {
    samples$seg.id[(sum(pd[0:(p-1)])+1):sum(pd[1:p])] <- p}

  rm(pd)

  #---------------------------------------------------------------------------#
  # 3. derive samples
  #---------------------------------------------------------------------------#

  # summarise samples
  samples = ddply(samples, .(seg.id,cell.id), summarise,
                  x=mean(x), y=mean(y),
                  start.time=min(date), end.time=max(date),
                  elapsed.time=difftime(max(date), min(date), units="h"))

  # if preserve.revisits is FALSE, reduce to unique pixels
  if (!preserve.revisits) {
    samples = ddply(samples, .(cell.id), summarise,
                    x=mean(x), y=mean(y), start.time=min(start.time),
                    end.time=max(end.time), elapsed.time=sum(elapsed.time))
  }

  # convert reduce samples to vector before exporting
  return(vect(samples, geom=c("x","y"), crs=crs(y)))

}
