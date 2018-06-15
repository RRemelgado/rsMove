#' @title plotMove
#'
#' @description {Standardized plotting of environmental and temporal information for a set of coordinate pairs.}
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param size.var Optional. Vector with elapsed time as report by \code{\link{moveReduce}}, \code{\link{sampleMove}} or \code{\link{timeDir}}. Controls the point size.
#' @param fill.var Optional. Vector with environmental information. Controls the fill color.
#' @param var.type One of 'cont' or 'cat'. Defines the type of \emph{fill.var}.
#' @importFrom ggplot2 ggplot aes_string theme geom_bar scale_fill_gradientn xlab ylab theme_bw geom_point guides
#' scale_size_continuous scale_color_discrete scale_fill_gradientn scale_size_continuous guide_legend element_text element_blank
#' @importFrom grDevices colorRampPalette
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
#'  r <- raster(system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # observation time
#'  time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time), format="%Y/%m/%d %H:%M:%S")
#'
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(shortMove, time, r)
#'
#'  # query data
#'  ov <- extract(r, move.reduce$points)
#'
#'  # plot output
#'  x <- move.reduce$points@data$x
#'  y <- move.reduce$points@data$y
#'  et <- move.reduce$points@data$`Elapsed time (minutes)`
#'  op <- plotMove(x, y, size.var=et, fill.var=ov, var.type="cont")
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

plotMove <- function(x, y, size.var=NULL, fill.var=NULL, var.type=NULL) {

#----------------------------------------------------------------------------------------------------------#
# 1. Check input data
#----------------------------------------------------------------------------------------------------------#

  if (length(x)!=length(y)) {stop('"x" and "y" have different lengths')}
  if(!is.null(size.var)) {
    if (!is.numeric(size.var)) {stop('"size.var" is nof of a valid class')}
    if (length(size.var)!=length(x)) {stop('coordinates and "size.var" have different lengths')}}
  if (!is.null(fill.var)) {
    if (is.null(var.type)) {stop('"fill.var" is set. Please specify "var.type"')}
    if (!var.type%in%c('cont', 'cat')) {stop('"var.type" is not a recognized keyword')}
    if (length(fill.var)!=length(x)) {stop('coordinates and "fill.var" have different lengths')}}

  # abort function if no variable is provided
  if (is.null(size.var) & is.null(fill.var)) {
    warning('neither "size.var" or "fill.var" were specfied. aborted.')
    return()
  }

#----------------------------------------------------------------------------------------------------------#
# 2. Determine variable limits/breaks
#----------------------------------------------------------------------------------------------------------#

  # time breaks
  if (!is.null(size.var)) {
    mv <- round(max(size.var, na.rm=T))
    nc <- nchar(as.character(mv))
    m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
    mv <- mv / m
    tb <- round(mv)
    if (mv > tb) {tb <- (tb+0.2)*m} else {tb <- tb*m}
    tb <- seq(0, tb, m)}

  # color scheme for value
  if (!is.null(fill.var)) {cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))}

#----------------------------------------------------------------------------------------------------------#
# 3. Build plots
#----------------------------------------------------------------------------------------------------------#

  # time and environmental data
  if (!is.null(size.var) & !is.null(fill.var)) {

    # build data frame
    df <- data.frame(x=x, y=y, time=size.var, value=fill.var)

    # build plot
    if (var.type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="time", fill="value"), color="black", pch=21) +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) +
        scale_fill_gradientn(name="Value", colours=cr(10))}
    if (var.type=="cat") {
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
  if (!is.null(size.var) & is.null(fill.var)) {

    # build data frame
    df <- data.frame(x=x, y=y, time=size.var)

    # build plot
    p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="time"), color="black", pch=21) +
      guides(fill=FALSE, col=guide_legend(override.aes = list(shape=15, size=6))) +
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12))
  }

#----------------------------------------------------------------------------------------------------------#

  # only environmental data
  if (is.null(size.var) & !is.null(fill.var)) {

    # build data frame
    df <- data.frame(x=x, y=y, value=fill.var)

    # build plot
    if (var.type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", color="value"), color="black", pch=21) +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) +
        scale_fill_gradientn(name="Value", colours=cr(10))}
    if (var.type=="cat") {
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
