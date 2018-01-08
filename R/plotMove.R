#' @title plotMove
#'
#' @description {Standardized plotting of environmental and temporal information for a set of coordinate pairs.}
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param size.var Optional. Vector with elapsed time as report by \code{\link{moveReduce}}, \code{\link{sampleMove}} or \code{\link{timeDir}}. Controls the point size.
#' @param fill.var Optional. Vector with environmental information. Controls the fill color.
#' @param var.type One of 'cont' or 'cat'. Defines the type of \emph{fill.var}.
#' @import raster rgdal ggplot2
#' @seealso \code{\link{dataQuery}} \code{\link{moveReduce}}
#' @return A \emph{ggplot} object.
#' @details {This function was designed to extent on other functions such as \code{\link{dataQuery}}, which provides environmental
#' information, and \code{\link{moveReduce}}, which provides information on the time spent per sample. Using these two functions
#' as an example, \emph{plotMove} can represent the relation between the elapsed time and the change in environmental conditions.}
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
#'  time <- strptime(paste0(moveData@data$date, ' ', moveData@data$time), format="%Y/%m/%d %H:%M:%S")
#'
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(xy=moveData, obs.time=time, img=r)
#'
#'  # query data
#'  ov <- extract(r, move.reduce$points)
#'
#'  # plot output
#'  x <- move.reduce$points@coords[,1]
#'  y <- move.reduce$points@coords[,2]
#'  et <- move.reduce$points@data$elapsed.time
#'  op <- plotMove(x=x, y=y, time=et, value=ov[,1], type="cont")
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

plotMove <- function(x=x, y=y, time=NULL, value=NULL, type=NULL) {

#----------------------------------------------------------------------------------------------------------#
# 1. Check input data
#----------------------------------------------------------------------------------------------------------#

  if (length(x)!=length(y)) {stop('"x" and "y" have different lengths')}
  if(!is.null(time)) {
    if (!is.numeric(time)) {stop('"time" is nof of a valid class')}
    if (length(time)!=length(x)) {stop('coordinates and "time" have different lengths')}}
  if (!is.null(value)) {
    if (is.null(type)) {stop('"value" is set. Please specify "type"')}
    if (!type%in%c('cont', 'cat')) {stop('"type" is not a recognized keyword')}
    if (length(value)!=length(x)) {stop('coordinates and "value" have different lengths')}}

  # abort function if no variable is provided
  if (is.null(time) & is.null(value)) {
    warning('neither "time" or "value" were specfied. aborted.')
    return()
  }

#----------------------------------------------------------------------------------------------------------#
# 2. Determine variable limits/breaks
#----------------------------------------------------------------------------------------------------------#

  # time breaks
  if (!is.null(time)) {
    mv <- round(max(time, na.rm=T))
    nc <- nchar(as.character(mv))
    m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
    mv <- mv / m
    tb <- round(mv)
    if (mv > tb) {tb <- (tb+0.2)*m} else {tb <- tb*m}
    tb <- seq(0, tb, m)}

  # color scheme for value
  if (!is.null(value)) {cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))}

#----------------------------------------------------------------------------------------------------------#
# 3. Build plots
#----------------------------------------------------------------------------------------------------------#

  # time and environmental data
  if (!is.null(time) & !is.null(value)) {

    # build data frame
    df <- data.frame(x=x, y=y, time=time, value=value)

    # build plot
    if (type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="time", fill="value"), color="black", pch=21) +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) +
        scale_fill_gradientn(name="Value", colours=cr(10))}
    if (type=="cat") {
      df$value <- factor(df$value)
      p <- ggplot(df, aes_string(x="x", y="y", color="value", size="time")) + theme_bw() + geom_point() +
        guides(col=guide_legend(override.aes=list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) +
        scale_color_discrete(name="Class")}

  }

#----------------------------------------------------------------------------------------------------------#

  # only time
  if (!is.null(time) & is.null(value)) {

    # build data frame
    df <- data.frame(x=x, y=y, time=time)

    # build plot
    p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="time", fill="red"), color="black", pch=21) +
      guides(fill=FALSE, col=guide_legend(override.aes = list(shape=15, size=6))) +
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12))
  }

#----------------------------------------------------------------------------------------------------------#

  # only environmental data
  if (is.null(time) & !is.null(value)) {

    # build data frame
    df <- data.frame(x=x, y=y, value=value)

    # build plot
    if (type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", color="value"), color="black", pch=21) +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) +
        scale_fill_gradientn(name="Value", colours=cr(10))}
    if (type=="cat") {
      df$value <- factor(df$value)
      p <- ggplot(df, aes_string(x="x", y="y", color="value")) + theme_bw() + geom_point() +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name="Elapsed Time", breaks=tb) +
        scale_color_discrete(name="Class")}


  }

    return(p)

}
