#' @title timeDir
#'
#' @description Analysis of environmental change in time for a set of coordinate pairs.
#' @param xy Object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param obs.dates Object of class \emph{Date} with \emph{xy} observation dates.
#' @param env.data Object of class \emph{RasterStack} or \emph{RasterBrick} or \emph{data.frame}.
#' @param env.dates Object of class \emph{Date} with \emph{env.data} observation dates.
#' @param temporal.buffer two element vector with temporal window size (expressed in days).
#' @param stat.fun Output statistical metric.
#' @param min.count Minimum number of samples required by \emph{stat.fun}. Default is 2.
#' @importFrom raster crs extract
#' @importFrom stats lm
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot geom_point theme guides scale_fill_gradientn scale_size_continuous ylab xlab element_text element_blank geom_hist aes_string
#' @seealso \code{\link{spaceDir}} \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{vector} with a requested statistical metric for each point in \emph{xy} and informative plots.
#' @details {This function quantifies environmental change in time along a movement track. First, for each point in \emph{xy},
#' the function compares its observation date (\emph{obs.dates}) against the acquisition dates (\emph{env.dates}) of \emph{env.data}
#' to select non \emph{NA} timesteps within a predefined temporal window (\emph{temporal.buffer}). The user can adjust this window to
#' determine which images are the most important. For example, if one wishes to know how the landscape evolved up to the observation
#' date of the target sample, \emph{temporal.buffer} can be define as, e.g., c(30,0) forcing the function to only consider pixels recorded
#' within the previous 30 days. After selecting adequate temporal information for each data point, a statistical metric is estimated. This
#' statistical metric is specified by \emph{stat.fun}. By default, the function reports on the slope between the acquisition dates of \emph{env.data}
#' and their corresponding values. When providing a new function, set x for \emph{env.dates} and y for \emph{env.data}. The final output is a list consisting of:
#' \itemize{
#' \item{\emph{stats} - \emph{data.frame} with the estimated statistical metric for each data point.}
#' \item{\emph{hist.plot} - Histogram plot of the requested statistical metric. The bin size is the standard deviation of all estimated values.}
#' \item{\emph{point.plot} - Plot of the \emph{xy} showing the spatial variability of the requested statistical metric.}
#' }}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)
#'  r.stk <- stack(file)
#'  r.stk <- stack(r.stk, r.stk, r.stk) # dummy files for the example
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # raster dates
#'  r.dates <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'
#'  # sample dates
#'  obs.dates <- as.Date(shortMove@data$date)
#'
#'  # perform directional sampling
#'  of <- function(x,y) {lm(y~x)$coefficients[2]}
#'  time.env <- timeDir(r.stk, r.dates, obs.dates, c(30,30), xy=shortMove, stat.fun=of)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

timeDir <- function(env.data, env.dates, obs.dates, temporal.buffer, xy=NULL, stat.fun=NULL, min.count=2) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  # samples
  if (!is.null(xy)) {if (!class(xy)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}}

  # sample dates
  if (!exists('obs.dates')) {stop('"obs.dates" is missing')}
  if (class(obs.dates)[1]!='Date') {stop('"obs.dates" is nof of a valid class')}
  if (length(obs.dates)!=length(xy)) {stop('"xy" and "obs.dates" have different lengths')}

  # environmental data dates
  if (class(env.dates)[1]!='Date') {stop('"env.dates" is nof of a valid class')}

  # environmental data
  if (!class(env.data)[1]%in%c("RasterStack", "RasterBrick", "data.frame")) {stop('"env.data" is not of a valid class')}
  if (class(env.data)[1]%in%c("RasterStack", "RasterBrick")) {
    if (is.null(xy)) {stop('"env.data" is a raster object. Please define "xy"')}
    if (crs(xy)@projargs!=crs(env.data)@projargs) {stop('"xy" and "env.data" have different projections')}
    if (length(env.dates)!=nlayers(env.data)) {stop('"env.data" and "env.dates" have different lengths')}}
  if (class(env.data)[1]=='data.frame') {if (length(env.dates)!=ncol(env.data)) {stop('"env.data" and "env.dates" have different lengths')}}

  # time information
  if (!is.numeric(temporal.buffer)) {stop('"temporal.buffer" us not numeric')}
  if (length(temporal.buffer)!=2) {stop('"temporal.buffer" does not have two elements')}

  # check/define input metrics
  if (is.null(stat.fun)) {stat.fun <- function(x,y) {lm(y~x)$coefficients[2]}} else {
    if(!is.function(stat.fun)) {stop('"stat.fun" is not a valid function')}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. retrieve environmental data
#-------------------------------------------------------------------------------------------------------------------------------#

  if (!is.data.frame(env.data)) {

    # retrieve environmental variables
    ind <- which(env.dates%in%seq.Date(min(obs.dates-temporal.buffer[1]), max(obs.dates+temporal.buffer[2]), by=1))
    env.data <- extract(env.data[[ind]], xy@coords)
    env.dates <- env.dates[ind]

  }

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. apply sampling approach
#-------------------------------------------------------------------------------------------------------------------------------#

  f <- function(i) {
    ind <- which(env.dates >= (obs.dates[i]-temporal.buffer[1]) & env.dates <= (obs.dates[i]+temporal.buffer[2]) & !is.na(env.data[i,]))
    if (length(ind) >= min.count) {
      x <- as.numeric(env.dates[ind])
      y <- as.numeric(env.data[i,ind])
      return(stat.fun(x,y))
    } else {return(NA)}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. query samples
#-------------------------------------------------------------------------------------------------------------------------------#

  df <- data.frame(value=unlist(lapply(1:nrow(env.data), f)))

#-------------------------------------------------------------------------------------------------------------------------------#
# 5. build plot(s)
#-------------------------------------------------------------------------------------------------------------------------------#

  p1 <- ggplot(df, aes_string(x="value")) + theme_bw() +
    geom_histogram(binwidth=sd(df$value, na.rm=TRUE),
                   aes(y=..count../sum(..count..))) +
    ylab('Relative freq. (%)') + xlab("Value")

  # build plot object
  if (!is.null(xy)) {
    cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))
    df0 <- data.frame(x=xy@coords[,1], y=xy@coords[,2], value=df$value)
    p2 <- ggplot(df0) + theme_bw() + xlab('X') + ylab('Y') +
      geom_point(aes_string(x="x", y="y", size="value", fill="value"), color="black", pch=21) +
      scale_size_continuous(guide=FALSE) + guides(col=cr(10)) +
      scale_fill_gradientn(colours=cr(10)) +
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    return(list(stats=df, hist.plot=p1, point.plot=p2))} else {return(list(stats=df, hist.plot=p1))}

}
