#' @title timeDir
#'
#' @description Temporal directional raster analysis.
#' @param xy Object of class "SpatialPoints" or "SpatialPointsDataFrame".
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterStack} or \emph{RasterBrick}.
#' @param edata Object of class \emph{data frame}.
#' @param rt Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{img} observation dates.
#' @param mws Moving window size (expressed in days).
#' @param dir One of \emph{fwd}, \emph{bwd} or \emph{both}. Default is \emph{both}.
#' @param fun Summary function.
#' @import raster rgdal
#' @importFrom stats lm
#' @seealso \code{\link{spaceDir}} \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{vector}.
#' @details {This function evaluates how do environmental conditions change in time 
#' along a movement track. First, it compares the observation times (\emph{ot}) 
#' of \emph{xy} against the acquisition times (\emph{rt}) of \emph{img} to search for relevant 
#' information within a pre-defined temporal window (\emph{mws}). The user can chose to 
#' only consider time steps before (\emph{bwd}) or after (\emph{fwd} the target observation 
#' time or look at both directios (\emph{both}). If the latest is chosen, the function 
#' applies \emph{mws} equally to both directions. After selecting adequate temporal 
#' information, a statistical metric (\emph{fun}) is used to summarize the selected 
#' time steps. By default, the slope will be used. The slope is estimated from 
#' a linear regression estimated between the acquisition times of \emph{img} and their 
#' corresponding values. When providing a new function, set x for time and y for 
#' the raster values. The output reports on:
#' \itemize{
#'  \item{\emph{x} - mean x coordinates}
#'  \item{\emph{y} - mean y coordinates}
#'  \item{\emph{timestamp} - mean observation time}
#'  \item{\emph{pixel.time} - elapsed time within a pixel for a given segment}
#'  \item{\emph{stat}: statistical metric}}
#' If \emph{edata} is provided, \emph{img} will only be used as a reference grid as \emph{edata} 
#' will contain the environmental data with each column representing a different variable. Otherwise, 
#' this data will be retrieved from \emph{img}.}
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
#'  rd <- seq.Date(as.Date("2013-08-01"), as.Date("2013-08-09"), 1)
#'  
#'  # sample dates
#'  td <- as.Date(moveData@data$date)
#'  
#'  # perform directional sampling
#'  of <- function(x,y) {lm(y~x)$coefficients[2]}
#'  t.sample <- timeDir(xy=moveData, ot=td, img=rsStk, rt=rd, mws=10, dir="bwd", fun=of)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

timeDir <- function(xy=xy, ot=ot, img=img, edata=NULL, rt=rt, mws=NULL, dir=NULL, fun=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  
  # sample dates
  if (!exists('ot')) {stop('"ot" is missing')}
  if (!class(ot)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"ot" is nof of a valid class')}
  if (length(ot)!=length(xy)) {stop('errorr: "xy" and "ot" have different lengths')}
  
  # raster
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}   
  
  # environmental data
  if (!is.null(edata)) {
    if (class(edata)[1]!='data.frame') {stop('"edata" provided but not a data frame')}
    if (nrow(edata)!=length(xy)) {stop('number of elements in "xy" and "edata" do not match')}}
  
  # raster dates
  if (!exists('rt')) {stop('"rt" is missing')}
  if (!class(rt)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"rt" is nof of a valid class')}
  if (length(rt)!=nlayers(img)) {stop('errorr: "img" and "rt" have different lengths')}
  
  # time information
  if (is.null(mws)) {stop('"mws" is missing')} else {
     if (!is.numeric(mws)) {stop('"mws" us not numeric')}}
  
  # query type
  if (!is.null(dir)) {
    if (length(dir)>1) {stop('"dir" has too many entries')}
    if (!dir%in%c('fwd', 'bwd', 'both')) {stop('"dir" is not a valid entry')}
  } else {dir <- 'both'}
  
  # check/define input metrics
  if (is.null(fun)) {fun <- function(x,y) {lm(y~x)$coefficients[2]}} else {
    if(!is.function(fun)) {stop('"fun" is not a valid function')}}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 2. retrieve environmental data
#-------------------------------------------------------------------------------------------------------------------------------#
  
  if (is.null(edata)) {
    
    # retrieve environmental variables
    ot <- as.Date(ot)
    rt <- as.Date(rt)
    ind <- which(rt%in%seq.Date(min(ot-mws), max(ot+mws), by=1))
    edata <- extract(img[[ind]], xy@coords)
    rt <- rt[ind]
    
    rm(img)
    
  }
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 3. apply sampling approach
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # backwards sampling
  if (dir=='bwd') {
    f <- function(i) {
      ind <- which(rt >= (ot[i]-mws) & rt <= ot[i])
      x <- as.numeric(rt[ind])
      y <- edata[i,]
      u <- !is.na(y)
      if (sum(u)>1) {return(fun(x[u],y[u]))} else {return(NA)}}}
  
  # forward sampling    
  if (dir=='fwd') {
    f <- function(i) {
      ind <- which(rt >= ot & rt <= (ot[i]+mws))
      x <- as.numeric(rt[ind])
      y <- edata[i,]
      u <- !is.na(y)
      if (sum(u)>1) {return(fun(x[u],y[u]))} else {return(NA)}}}
    
  # Backward-Forward sampling
  if (dir=='both') {
  f <- function(i) {
    ind <- which(rt >= (ot[i]-mws) & rt <= (ot[i]+mws))
    x <- as.numeric(rt[ind])
    y <- edata[i,]
    u <- !is.na(y)
    if (sum(u)>1) {return(fun(x[u],y[u]))} else {return(NA)}}}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 4. query samples
#-------------------------------------------------------------------------------------------------------------------------------#
  
  df <- data.frame(value=unlist(lapply(1:length(xy), f)))
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 5. build plot
#-------------------------------------------------------------------------------------------------------------------------------#
  
  
  # determine y scale range
  mv <- max(df$value)
  nc <- nchar(as.character(mv))
  m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
  mv <- mv / m
  yr <- round(mv)
  if (mv > yr) {yr <- (yr+0.2)*m} else {yr <- yr*m}
  
  # build plot object
  p <- ggplot(df, aes(y=value)) + geom_bar(stat='identity') +  ylim(0,yr) + 
    xlab('') + theme(axis.text.x=element_blank())
  
  
  return(list(stats=df, plot=p))
 
}