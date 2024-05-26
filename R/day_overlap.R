#' @title day_overlap
#'
#' @description Estimate the number of the day covered by timestamps.
#' @param x Observation time. Object of class \emph{POSIXct}, \emph{POSIXlt}, or \emph{POSIXt}.
#' @details {The function estimates how much of a day is covered
#' by the timestamps in \emph{x}, which correspond to the timestamps
#' in a movement dataset. Note that \emph{x} is not accepted if
#' containing data from more than one day.}
#' @return {A \emph{list} containing:
#'  \itemize{
#'  \item{\emph{day.cover} - Percent of the day with recorded timestamps.}
#'  \item{\emph{day.overlap} - Percent of the day covered between the first and last timestamp.}
#' }
#' }
#' @examples {
#'
#' # load samples
#' multiMove <- read.csv(system.file('extdata', 'multiMove.csv', package="rsMove"))
#'
#' # data-time of species observations
#' times = strptime(multiMove$timestamp,format="%Y-%m-%d %H:%M:%S")
#'
#' # subset to first unique day
#' days = as.Date(times)
#' times = times[which(days == unique(days)[1])]
#'
#' # extract regions
#' day_overlap(date)])
#'
#' }
#' @export

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

day_overlap = function(x) {

  #---------------------------------------------------------------------------#
  # 1. check if input is a valid vector of timestamps
  #---------------------------------------------------------------------------#

  if (!max(class(x) %in% c("POSIXlt","POSIXt","POSIXct"))) stop('"x" is not of a valid class')
  if (length(unique(format(x, "%Y-%m-%d"))) > 1) stop('"x" must refer to a single day')

  #---------------------------------------------------------------------------#
  # 2. estimate proportion of the day covered by timestamps
  #---------------------------------------------------------------------------#

  # convert first timestamp into number of seconds in the day
  t0 = as.numeric(format(x, "%H"))*3600 +
    as.numeric(format(x, "%M"))*60 +
    as.numeric(format(x, "%S"))

  # translate t0 to a hypothetical sequence of seconds (defines max. time cover)
  t1 = seq(min(t0), max(t0), 1)

  # build sequence for a hypothetical complete day
  t2 = seq(0, 86400, 1)

  # estimate percent overlap
  return(
    list(
      day.cover=length(intersect(t0,t2))/length(t2)*100,
      day.overlap=length(intersect(t1,t2))/length(t2)*100
      )
    )

}
