#' @title plotMove
#'
#' @description {Standardized plotting of environmental and temporal information for a set of coordinate pairs.}
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param size.var Optional. Controls the point size.
#' @param fill.var Optional. Controls the fill color.
#' @param var.type One of 'cont' or 'cat'. Defines the type of \emph{fill.var}.
#' @param var.names Character vector with names for \emph{size.var} and \emph{fill.var} to ve added to the plot.
#' @importFrom ggplot2 ggplot aes_string theme geom_bar scale_fill_gradientn xlab ylab theme_bw geom_point guides
#' scale_size_continuous scale_color_discrete scale_fill_gradientn scale_size_continuous guide_legend element_text element_blank
#' @importFrom grDevices colorRampPalette
#' @seealso \code{\link{dataQuery}} \code{\link{moveReduce}}
#' @return A \emph{ggplot} object.
#' @details {This function was designed to extent on other functions such as \code{\link{dataQuery}}, which provides environmental
#' information, and \code{\link{moveReduce}}, which provides information on e.g. the time spent per sample. Using these two functions
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
#'  move.reduce <- moveReduce(shortMove, r, time)
#'
#'  # query data
#'  ov <- extract(r, move.reduce$points)
#'
#'  # plot output
#'  x <- move.reduce$points@data$x
#'  y <- move.reduce$points@data$y
#'  et <- move.reduce$points@data$elapsed.time
#'  op <- plotMove(x, y, size.var=et, fill.var=ov, var.type="cont")
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

plotMove <- function(x, y, size.var=NULL, fill.var=NULL, var.type=NULL, var.names=NULL) {

  #----------------------------------------------------------------------------------------------------------#
  # 1. Check input data
  #----------------------------------------------------------------------------------------------------------#

  # check coordinates
  if (!is.numeric(x) | !is.numeric(y)) {stop('"x" and/or "y" not numeric')}
  if (length(x)!=length(y)) {stop('"x" and "y" have different lengths')}

  # check size variable
  if(!is.null(size.var)) {
    if (!is.numeric(size.var)) {stop('"size.var" is nof of a valid class')}
    if (length(size.var)!=length(x)) {stop('coordinates and "size.var" have different lengths')}}

  # check fill variable
  if (!is.null(fill.var)) {
    if (is.null(var.type)) {stop('"fill.var" is set. Please specify "var.type"')}
    if (!var.type%in%c('cont', 'cat')) {stop('"var.type" is not a recognized keyword')}
    if (length(fill.var)!=length(x)) {stop('coordinates and "fill.var" have different lengths')}}

  # check variable names
  if (is.null(var.names)) {
    if (var.type == 'cont') {var.names <- c('size', 'value')}
    if (var.type == 'cat') {var.names <- c('size', 'class')}
  } else {
    if (!is.character(var.names)) {stop('"var.names" is not a character object')}
    if (length(var.names) != 2) {stop('"var.names" should have two elements')}
  }

  # abort function if no variable is provided
  if (is.null(size.var) & is.null(fill.var)) {
    warning('neither "size.var" or "fill.var" were specfied. aborted.')
    return()
  }

  #----------------------------------------------------------------------------------------------------------#
  # 2. Determine variable limits/breaks
  #----------------------------------------------------------------------------------------------------------#

  # size breaks
  if (!is.null(size.var)) {tb <- seq(min(size.var), max(size.var), sd(size.var))}

  # color scheme for fill
  if (!is.null(fill.var)) {cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))}

  #----------------------------------------------------------------------------------------------------------#
  # 3. Build plots
  #----------------------------------------------------------------------------------------------------------#

  # size and environmental data
  if (!is.null(size.var) & !is.null(fill.var)) {

    # build data frame
    df <- data.frame(x=x, y=y, size=size.var, fill=fill.var)

    # build plot
    if (var.type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="size", fill="fill"), color="black", pch=21) +
        guides(size=guide_legend(order=1, override.aes=list(shape=1)), col=guide_legend(order=0, override.aes=list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name=var.names[1], breaks=tb, range=c(2,12), labels=sprintf("%0.2f", tb)) +
        scale_fill_gradientn(name=var.names[2], colours=cr(10))}
    if (var.type=="cat") {
      df$fill <- factor(df$fill)
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="size", color="fill")) +
        guides(size=guide_legend(order=1, override.aes=list(shape=1)), col=guide_legend(order=0, override.aes=list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name=var.names[1], breaks=tb, range=c(2,12), labels=sprintf("%0.2f", tb)) +
        scale_color_discrete(name=var.names[2])}

  }

  #----------------------------------------------------------------------------------------------------------#

  # only size
  if (!is.null(size.var) & is.null(fill.var)) {

    # build data frame
    df <- data.frame(x=x, y=y, size=size.var)

    # build plot
    p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", size="size"), color="black", pch=21) +
      guides(fill=FALSE, col=guide_legend(override.aes = list(shape=15, size=6))) +
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      scale_size_continuous(name=var.names[1], breaks=tb, range=c(2,12), labels=sprintf("%0.2f", tb))
  }

  #----------------------------------------------------------------------------------------------------------#

  # only environmental data
  if (is.null(size.var) & !is.null(fill.var)) {

    # build data frame
    df <- data.frame(x=x, y=y, fill=fill.var)

    # build plot
    if (var.type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes_string(x="x", y="y", color="fill"), color="black", pch=21) +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name=var.names[1], breaks=tb, range=c(2,12)) +
        scale_fill_gradientn(name=var.names[2], colours=cr(10))}
    if (var.type=="cat") {
      df$fill <- factor(df$fill)
      p <- ggplot(df, aes_string(x="x", y="y", color="fill")) + theme_bw() + geom_point() +
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) +
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(),
              panel.grid.minor=element_blank()) +
        scale_size_continuous(name=var.names[1], breaks=tb, labels=sprintf("%0.2f", tb)) +
        scale_color_discrete(name=var.names[2])}


  }

  return(p)

}

