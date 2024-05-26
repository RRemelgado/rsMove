#' @title spaceDir
#'
#' @description Analysis of environmental change in space along a movement track.
#' @param x \emph{SpatVector} object.
#' @param y Object of class \emph{SpatRaster}.
#' @param z Object of class \emph{Date} or \emph{POSIXct} with observation dates for each entry in \emph{x}.
#' @param space.buffer Spatial buffer size expressed in meters.
#' @param time.buffer Temporal buffer size expressed in days.
#' @param fun List of functions to apply to each time step.
#' @param min.count Minimum number of pixels required by \emph{stat.fun}. Default is 2.
#' @importFrom terra extract distance
#' @seealso \code{\link{timeDir}}
#' @details {The function quantifies environmental changes along a spatial
#' gradient defined by GPS tracking dataset. For each GPS observation, the
#' function finds the nearest observations within a given \emph{buffer.size},
#' and uses those observations to apply a user-define list of functions
#' (\emph{fun}) applied to thethe underlying pixel values in \emph{y}. In
#' addition, the function will report on the linear distance traveled
#' between endpoints (in meters) and the associated travel time (in minutes).}
#' @return {A \emph{data.frame} with statistics at each GPS observation.
#' In addition the user defined statistics, the table will contain several
#' entries describing the data used to calculate those statistics:
#' \itemize{
#'  \item{\emph{distance_traveled} - Total distance traveled between all observations}
#'  \item{\emph{elapsed_time} - Time between the first and last observation}
#'  \item{\emph{nr_observations} - nr_observations}
#'  \item{\emph{pixel_value} - Pixel value at the main observation}}}
#'
#' @examples {
#'
#'  require(terra)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  shortMove <- read.csv(system.file('extdata',
#'  'shortMove.csv', package="rsMove"))
#'
#'  # convert observations to vector
#'  shortMove = vect(shortMove, geom=c("x","y"), crs="EPSG:32632")
#'
#'  # observation time
#'  obs.time <- strptime(paste0(shortMove$date, ' ',shortMove$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # construct target functions
#'  functions <- list(
#'     slope=function(i) lm(i~c(1:length(i)))$coefficients[2][[1]],
#'     mean=function(i) mean(i, na.rm=T),
#'     sd=function(i) sd(i, na.rm=T)
#'  )
#'  s.sample <- spaceDir(shortMove, r, obs.time, 30, 1, fun=functions)
#'
#' }
#' @export

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

spaceDir <- function(x, y, z, space.buffer, time.buffer, fun, min.count=2) {

  #---------------------------------------------------------------------------#
  # 1. check variables
  #---------------------------------------------------------------------------#

  # samples
  if (!class(x)%in%c('SpatVector')) stop('"x" is not of a valid class')

  # sample dates
  if (!class(z)[1]%in%c('Date', 'POSIXct')) stop('"z" not of a valid class')
  if (length(z)!=nrow(x)) stop('"x" shorter than "z"')

  # raster
  if (!class(y)[1]=='SpatRaster') stop('"y" is not of a valid class')
  if (crs(x)!=crs(y)) stop('"x" and "y" have different projections')

  # check input metrics
  if (!is.list(fun)) stop('"fun" must be a list of functions') else {
    fun_type = sapply(fun, function(f) class(f))
    if (sum(fun_type == "function") != length(fun)) {
      stop('one or more elements in "fun" not a function')
    }
  }

  # check min.count
  if (!is.numeric(min.count)) stop('"min.count" must be a numeric element')

  #---------------------------------------------------------------------------#
  # 2. evaluate space-time dynamics in movement data
  #---------------------------------------------------------------------------#

  output = do.call(rbind, lapply(1:nrow(x), function(i) {

    # target observations
    ind = which(
      (distance(x[i,], x) < buffer.size) &
        (abs(difftime(z[i],z, units="d")) < time.buffer))

    if (length(ind) > min.count) {

      # extract pixel values
      pixel_values = extract(y, x[ind,], ID=F)[[1]]

      # evaluate extracted pixel values
      tmp = as.data.frame(lapply(fun, function(f) f(pixel_values)[[1]]))
      tmp$distance_travelled=sum(sapply(2:length(ind), function(p) distance(y[(p-1):p])))
      tmp$elapsed_time=difftime(max(z[ind]), min(z[ind]), units="d")
      tmp$nr_observations=length(ind)
      tmp$pixel_value = extract(y, x[i,], mean, ID=F)[[1]]

    } else {

      # return empty data.frame
      tmp = as.data.frame(matrix(NA,1,length(names(fun))+3))
      colnames(tmp) = c(names(fun),
                        "distance_traveled",
                        "elapsed_time",
                        "nr_observations",
                        "pixel_value")

    }

    return(tmp)

  }))

  # export
  return(output)

}
