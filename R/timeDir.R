#' @title timeDir
#'
#' @description Analysis of environmental change in time for a set of coordinate pairs.
#' @param xy Object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param obs.dates Object of class \emph{Date} with \emph{xy} observation dates.
#' @param img Object of class
#' @param env.data Object of class \emph{RasterStack} or \emph{RasterBrick} or \emph{data.frame}.
#' @param env.dates Object of class \emph{Date} with \emph{env.data} observation dates.
#' @param window.size Temporal moving window size (expressed in days).
#' @param sample.direction One of \emph{forward}, \emph{backward} or \emph{both}. Default is \emph{both}.
#' @param stat.fun Output statistical metric.
#' @param min.count Minimum number of samples required by \emph{stat.fun}. Default is 2.
#' @import raster rgdal
#' @importFrom stats lm
#' @seealso \code{\link{spaceDir}} \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{vector} with a requested statistical metric for each point in \emph{xy}.
#' @details {This function evaluates how environmental conditions change in time along a movement track.
#' First, for each point in \emph{xy}, the function compares its observation date (\emph{obs.dates}) against
#' the acquisition dates (\emph{env.dates}) of \emph{env.data} to select non \emph{NA} timesteps within a
#' pre-defined temporal window (\emph{window.size}). The user can chose to only consider time steps before
#' (\emph{backward}) or after (\emph{forward} the target observation time or look at both directios (\emph{both}).
#' If the latest is chosen, the function applies \emph{window.size} equally to both directions. After selecting
#' adequate temporal information for each data point, a statistical metric is estimated. The statistical metric
#' is provided by (\emph{stat.fun}). By default, the slope is reported from a linear regression between the
#' acquisition times of \emph{env.data} and their corresponding values. When providing a new function, set x
#' for \emph{env.dates} and y for \emph{env.data}.}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'  rsStk <- stack(rsStk, rsStk, rsStk) # dummy files for the example
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(rsStk))
#'
#'  # raster dates
#'  r.dates <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'
#'  # sample dates
#'  obs.dates <- strptime(paste0(moveData@data$date, ' ',moveData@data$time), format="%Y/%m/%d %H:%M:%S")
#'
#'  # perform directional sampling
#'  of <- function(x,y) {lm(y~x)$coefficients[2]}
#'  time.env <- timeDir(xy=moveData, obs.dates=obs.dates, img=rsStk, env.dates=r.dates, window.size=10, sample.direction="backward", stat.fun=of)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

timeDir <- function(xy=NULL, obs.dates=obs.dates, img=NULL, env.data=NULL, env.dates=env.dates, window.size=NULL, sample.direction=NULL, stat.fun=NULL, min.count=2) {

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
  if (is.null(window.size)) {stop('"window.size" is missing')} else {if (!is.numeric(window.size)) {stop('"window.size" us not numeric')}}

  # query type
  if (!is.null(sample.direction)) {
    if (length(sample.direction)>1) {stop('"sample.direction" has too many entries')}
    if (!sample.direction%in%c('forward', 'backward', 'both')) {stop('"sample.direction" is not a valid entry')}
  } else {sample.direction <- 'both'}

  # check/define input metrics
  if (is.null(stat.fun)) {stat.fun <- function(x,y) {lm(y~x)$coefficients[2]}} else {
    if(!is.function(stat.fun)) {stop('"stat.fun" is not a valid function')}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. retrieve environmental data
#-------------------------------------------------------------------------------------------------------------------------------#

  if (!is.data.frame(env.data)) {

    # retrieve environmental variables
    ind <- which(env.dates%in%seq.Date(min(obs.dates-window.size), max(obs.dates+window.size), by=1))
    env.data <- extract(img[[ind]], xy@coords)
    env.dates <- env.dates[ind]

    rm(img)

  }

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. apply sampling approach
#-------------------------------------------------------------------------------------------------------------------------------#

  # backwards sampling
  if (sample.direction=='backward') {
    f <- function(i) {
      ind <- which(env.dates >= (obs.dates[i]-window.size) & env.dates <= obs.dates[i])
      x <- as.numeric(env.dates[ind])
      y <- env.data[i,]
      u <- !is.na(y)
      if (sum(u) >= min.count) {return(stat.fun(x[u],y[u]))} else {return(NA)}}}

  # forward sampling
  if (sample.direction=='forward') {
    f <- function(i) {
      ind <- which(env.dates >= obs.dates & env.dates <= (obs.dates[i]+window.size))
      x <- as.numeric(env.dates[ind])
      y <- env.data[i,]
      u <- !is.na(y)
      if (sum(u) >= min.count) {return(stat.fun(x[u],y[u]))} else {return(NA)}}}

  # Backward-Forward sampling
  if (sample.direction=='both') {
  f <- function(i) {
    ind <- which(env.dates >= (obs.dates[i]-window.size) & env.dates <= (obs.dates[i]+window.size))
    x <- as.numeric(env.dates[ind])
    y <- env.data[i,]
    u <- !is.na(y)
    if (sum(u) >= min.count) {return(stat.fun(x[u],y[u]))} else {return(NA)}}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. query samples
#-------------------------------------------------------------------------------------------------------------------------------#

  df <- data.frame(value=unlist(lapply(1:nrow(env.data), f)))

#-------------------------------------------------------------------------------------------------------------------------------#
# 5. build plot
#-------------------------------------------------------------------------------------------------------------------------------#

  # build plot object
  if (!is.null(xy)) {
    cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))
    df0 <- data.frame(x=xy@coords[,1], y=xy@coords[,2], value=df$value)
    p <- ggplot(df0) + theme_bw() + xlab('X') + ylab('Y') +
      geom_point(aes_string(x="x", y="y", size="value", fill="value"), color="black", pch=21) +
      scale_size_continuous(guide=FALSE) + guides(col=cr(10)) +
      scale_fill_gradientn(colours=cr(10)) +
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    return(list(stats=df, plot=p))} else {return(list(stats=df))}

}
