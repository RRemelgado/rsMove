#' @title timeDir
#'
#' @description Analysis of environmental change in time for a set of coordinate pairs.
#' @param x Object of class \emph{SpatRaster}.
#' @param x.dates Object of class \emph{Date} with observation dates of \emph{y}.
#' @param y Object of class \emph{SpatVector}.
#' @param y.dates Object of class \emph{Date} with observation dates of \emph{y}.
#' @param temporal.buffer two element vector with temporal window size (expressed in days).
#' @param fun List of statistical function to apply to the data.
#' @param min.count Minimum number of samples required by \emph{stat.fun}. Default is 2.
#' @importFrom terra extract geom
#' @importFrom stats median
#' @seealso \code{\link{spaceDir}}
#' @details {This function quantifies environmental changes in time for each
#' GPS entry along a movement track. First, for each point in \emph{y}, the
#' function compares its observation date (\emph{y.dates}) against the dates
#' with environmental data (\emph{x.dates}), and preserves those entries in
#' \emph{x} that fall within the \emph{temporal.buffer}. The user can adjust
#' this window to determine which images are the most important. For example,
#' if one wishes to know how the landscape evolved up to the observation date
#' of the target sample, \emph{temporal.buffer} can be define as, e.g., c(30,0)
#' forcing the function to only consider pixels recorded within the previous 30
#' days. After selecting adequate temporal information for each data point, a
#' list of user-defined statistical metrics are calculated (i.e., \emph{fun}).}
#' @return {A \emph{data.frame} with statistics at each GPS observation.
#' In addition the user defined statistics, the table will contain several
#' entries describing the data used to calculate those statistics:
#' \itemize{
#'  \item{\emph{distance_traveled} - Total distance traveled between all observations}
#'  \item{\emph{elapsed_time} - Time between the first and last observation}
#'  \item{\emph{nr_observations} - nr_observations}
#'  \item{\emph{pixel_value} - Pixel value at the main observation}}}
#' @examples {
#'
#'  require(terra)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)
#'  r.stk <- rast(file)
#'  r.stk <- c(r.stk, r.stk, r.stk) # dummy files for the example
#'
#'  # read movement data
#'  shortMove <- read.csv(system.file('extdata', 'shortMove.csv', package="rsMove"))
#'
#'  # convert observations to vector
#'  shortMove = vect(shortMove, geom=c("x","y"), crs="EPSG:32632")
#'
#'  # raster dates
#'  r.dates <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'
#'  # sample dates
#'  obs.dates <- as.Date(shortMove$date)
#'
#'  # perform directional sampling
#'  functions <- list(
#'     slope=function(x,y) lm(y~x)$coefficients[2][[1]],
#'     mean=function(x,y) mean(y, na.rm=T),
#'     sd=function(x,y) sd(y, na.rm=T)
#'  )
#'  time.env <- timeDir(r.stk, r.dates, shortMove, obs.dates, 3, fun=functions)
#'
#' }
#' @export

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

timeDir <- function(x, x.dates, y, y.dates, temporal.buffer, fun=NULL, min.count=2) {

  #---------------------------------------------------------------------------#
  # 1. check variables
  #---------------------------------------------------------------------------#

  # check raster data
  if (class(x)!=c('SpatRaster')) stop('"x" is not of a valid class')
  if (class(x.dates)!=c('Date')) stop('"x.dates" is not of a valid class')

  # check vector data
  if (class(y)!=c('SpatVector')) stop('"y" is not of a valid class')
  if (class(y.dates)!=c('Date')) stop('"y.dates" is not of a valid class')

  # time information
  if (!is.numeric(temporal.buffer)) stop('"temporal.buffer" us not numeric')

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
  # 2. retrieve environmental data
  #---------------------------------------------------------------------------#

  # retrieve environmental variables
  ind <- which(x.dates%in%seq.Date(min(y.dates-temporal.buffer),
                                   max(y.dates+temporal.buffer), by=1))
  env.data <- extract(x[[ind]], geom(y)[,c("x","y")])
  x.dates = x.dates[ind]

  #---------------------------------------------------------------------------#
  # 3. evaluate changes for each sampling location
  #---------------------------------------------------------------------------#

  statistics = do.call(rbind, lapply(1:nrow(y), function(i) {

    # find observations within temporal window
    ind <- which(
      (x.dates >= (y.dates[i]-temporal.buffer)) &
        (x.dates <= (y.dates[i]+temporal.buffer)) &
        (!is.na(env.data[i,])))

    if (length(ind) >= min.count) {
      tx <- as.numeric(x.dates[ind])
      ty <- as.numeric(env.data[i,ind])
      tmp = as.data.frame(lapply(fun, function(f) f(tx,ty)))
      tmp$start.date = min(x.dates[ind])
      tmp$median.date = median(x.dates[ind])
      tmp$end.date = max(x.dates[ind])
      tmp$nr_steps = length(ind)
      return(tmp)
    } else {
      # return empty data.frame
      tmp = as.data.frame(matrix(NA,1,length(names(fun))+4))
      colnames(tmp) = c(names(fun),
                        "start.date", "median.date",
                        "end.date", "nr_steps")
      return(tmp)

    }
    }))

  return(statistics)

}
