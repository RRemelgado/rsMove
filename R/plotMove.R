#' @title plotMove
#'
#' @description Plots per pixel time and value.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param o.time Vector with time length.
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
#'  o.time <- strptime(paste0(moveData@data$date, ' ', moveData@data$time), format="%Y/%m/%d %H:%M:%S")
#'  
#'  # reduce amount of samples
#'  move.reduce <- moveReduce(xy=moveData, ot=o.time, img=r)
#'  
#'  # query data
#'  ov <- dataQuery(xy=move.reduce, img=r)
#'  
#'  plot output
#'  op <- plotMove(x=move.reduce@coords[,1], y=move.reduce@coords[,2], 
#'  o.time=move.reduce@data$time, value=ov, type="cont")
#'  
#' }
#' @export

#----------------------------------------------------------------------------------------------------------#

plotMove <- function(x=x, y=y, o.time=NULL, value=NULL, type=NULL) {
  
#----------------------------------------------------------------------------------------------------------#
# 1. Check input data
#----------------------------------------------------------------------------------------------------------#
  
  if (length(x)!=length(y)) {stop('"x" and "y" have different lengths')}
  if(!is.null(o.time)) {
    if (!is.numeric(o.time)) {stop('"o.time" is nof of a valid class')}
    if (length(o.time)!=length(x)) {stop('coordinates and "o.time" have different lengths')}}
  if (!is.null(value)) {
    if (is.null(type)) {stop('"value" is set. Please specify "type"')}
    if (!type%in%c('cont', 'cat')) {stop('"type" is not a recognized keyword')}
    if (length(value)!=length(x)) {stop('coordinates and "value" have different lengths')}}
  
  # abort function if no variable is provided
  if (is.null(o.time) & is.null(value)) {
    warning('neither "o.time" or "value" were specfied. aborted.')
    return()
  }
  
#----------------------------------------------------------------------------------------------------------#
# 2. Determine variable limits/breaks
#----------------------------------------------------------------------------------------------------------#  
  
  # time breaks
  if (!is.null(o.time)) {
    mv <- round(max(o.time, na.rm=T))
    nc <- nchar(as.character(mv))
    m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
    mv <- mv / m
    tb <- round(mv)
    if (mv > tb) {tb <- (tb+0.2)*m} else {tb <- tb*m}
    tb <- seq(0, tb, m)}
  
  # value limit
  if (!is.null(value)) {
    if (type=='cont') {
      mv <- round(max(value, na.rm=T))
      nc <- nchar(as.character(mv))
      m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
      mv <- mv / m
      vl <- round(mv)
      if (mv > vl) {vl <- (vl+0.2)*m} else {tb <- vl*m}}}
  
#----------------------------------------------------------------------------------------------------------#
# 3. Build plots
#----------------------------------------------------------------------------------------------------------#  
  
  # time and environmental data
  if (!is.null(o.time) & !is.null(value)) {
    
    # build data frame
    df <- data.frame(x=x, y=y, time=o.time, value=value)
    
    # build plot
    if (type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes(x=x, y=y, size=time, fill=value), color="black", pch=21) + 
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) + 
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank()) + 
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) + 
        scale_fill_gradientn(name="Value", breaks=c(0, (vl/2), vl), limits=c(0,vl))}
    if (type=="cat") {
      df$value <- factor(df$value)
      p <- ggplot(df, aes(x=x, y=y, color=value, size=time)) + theme_bw() + geom_point() + 
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) + 
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank()) + 
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) + 
        scale_color_discrete(name="Class")}
    
  }

#----------------------------------------------------------------------------------------------------------#
  
  # only time
  if (!is.null(o.time) & is.null(value)) {
    
    # build data frame
    df <- data.frame(x=x, y=y, time=o.time)
    
    # build plot
    p <- ggplot(df) + theme_bw() + geom_point(aes(x=x, y=y, size=time, fill="red"), color="black", pch=21) + 
      guides(fill=FALSE, col=guide_legend(override.aes = list(shape=15, size=6))) + 
      theme(legend.text=element_text(size=10), panel.grid.major=element_blank(), 
            panel.grid.minor=element_blank()) + 
      scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12))
  }
  
#----------------------------------------------------------------------------------------------------------#
    
  # only environmental data
  if (is.null(o.time) & !is.null(value)) {
    
    # build data frame
    df <- data.frame(x=x, y=y, value=value)
    
    # build plot
    if (type=="cont") {
      p <- ggplot(df) + theme_bw() + geom_point(aes(x=x, y=y, size=time, fill=value), color="black", pch=21) + 
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) + 
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank()) + 
        scale_size_continuous(name="Elapsed Time", breaks=tb, range=c(2,12)) + 
        scale_fill_gradientn(name="Value", breaks=c(0, (vl/2), vl), limits=c(0,vl))}
    if (type=="cat") {
      df$value <- factor(df$value)
      ggplot(df, aes(x=x, y=y, color=factor(value))) + theme_bw() + geom_point() + 
        guides(col=guide_legend(override.aes = list(shape=15, size=6))) + 
        theme(legend.text=element_text(size=10), panel.grid.major=element_blank(), 
              panel.grid.minor=element_blank()) + 
        scale_size_continuous(name="Elapsed Time", breaks=tb) + 
        scale_color_discrete(name="Class")}
    
    
  }
  
    return(p)
    
}
