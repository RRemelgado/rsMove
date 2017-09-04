#' @title plotMove
#'
#' @description {Standardized plotting of sampled environmental information
#'  and time spent for each sample given a set of coordinate pairs.}
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param obs.time Vector with time length.
#' @param value Vector with environmental data.edata Object of class \emph{RasterLayer} or \emph{data.frame}.
#' @param type One of 'cont' or 'cat'. Defines the type of \emph{value}.
#' @import raster rgdal ggplot2
#' @seealso \code{\link{dataQuery}} \code{\link{moveReduce}}
#' @return A \emph{ggplot} object.
#' @details {This function plots a provided set of x and y coordinates adding information on the elapsed time
#' at each coordinate (\emph{e.time}) and/or a given environmental variable (\emph{value}). If only \emph{time}
#' or \emph{value} are set, the function builds a scatterplot where the size of the points is defined by the input
#' variables. If both are provided, the size of the points is defined by the elapsed time and the raster value is
#' used to build a color scheme for the points. When \emph{value} is provided, the keyword \emph{type} is required.
#' It will influence how the plots are built. The breaks a limits for the point size and colors is define automatically.}
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
#'  obs.time <- strptime(paste0(moveData@data$date, ' ', moveData@data$time), format="%Y/%m/%d %H:%M:%S")
#'
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(xy=moveData, obs.time=obs.time, img=r)
#'
#'  # query data
#'  ov <- dataQuery(xy=move.reduce$points, img=r)
#'
#'  # plot output
#'  x <- move.reduce$points@coords[,1]
#'  y <- move.reduce$points@coords[,2]
#'  et <- move.reduce$points@data$elapsed.time
#'  op <- plotMove(x=x, y=y, obs.time=et, value=ov[,1], type="cont")
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

plotMove <- function(x=x, y=y, obs.time=NULL, value=NULL, type=NULL) {

#----------------------------------------------------------------------------------------------------------#
# 1. Check input data
#----------------------------------------------------------------------------------------------------------#

  if (length(x)!=length(y)) {stop('"x" and "y" have different lengths')}
  if(!is.null(obs.time)) {
    if (!is.numeric(obs.time)) {stop('"obs.time" is nof of a valid class')}
    if (length(obs.time)!=length(x)) {stop('coordinates and "obs.time" have different lengths')}}
  if (!is.null(value)) {
    if (is.null(type)) {stop('"value" is set. Please specify "type"')}
    if (!type%in%c('cont', 'cat')) {stop('"type" is not a recognized keyword')}
    if (length(value)!=length(x)) {stop('coordinates and "value" have different lengths')}}

  # abort function if no variable is provided
  if (is.null(obs.time) & is.null(value)) {
    warning('neither "obs.time" or "value" were specfied. aborted.')
    return()
  }

#----------------------------------------------------------------------------------------------------------#
# 2. Determine variable limits/breaks
#----------------------------------------------------------------------------------------------------------#

  # time breaks
  if (!is.null(obs.time)) {
    mv <- round(max(obs.time, na.rm=T))
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
  if (!is.null(obs.time) & !is.null(value)) {

    # build data frame
    df <- data.frame(x=x, y=y, time=obs.time, value=value)

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
  if (!is.null(obs.time) & is.null(value)) {

    # build data frame
    df <- data.frame(x=x, y=y, time=obs.time)

    # build plot
    p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="time", fill="red"), color="black", pch=21) +
      guides(fill=FALSE, col=guide_legend(override.aes = list(shape=15, size=6))) +
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12))
  }

#----------------------------------------------------------------------------------------------------------#

  # only environmental data
  if (is.null(obs.time) & !is.null(value)) {

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
