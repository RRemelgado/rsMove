#' @title spaceDir
#'
#' @description Spatial directional raster analysis.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterLayer}
#' @param edata Object of class \emph{vector}.
#' @param dir One of \emph{fwd}, \emph{bwd} or \emph{both}. Default is \emph{both}.
#' @param fun Summary function.
#' @import raster sp rgdal
#' @importFrom stats lm
#' @seealso \code{\link{timeDir}} \code{\link{dataQuery}} \code{\link{imgInt}
#' @return A \emph{list}.
#' @details {This function evaluates how do environmental conditions change in space 
#' along a movement track. Before an output is returned, the function looks for segments 
#' of consecutive samples that appear within the same pixel and provides mean values for 
#' their coordinates and observation times as well as the corresponding pixel value and 
#' the elapsed time during the segment (in minutes). this step avoids the replication of 
#' samples while preserving observation related to revisits to the same pixels. Then, for 
#' each observation, the function looks at it imediate neighbors to define spatial segments 
#' over which a statistical metric will be estimated. All the pixels between the two endpoints 
#' of the segment are considered in this process. By defaut, the slope of the linear regression 
#' between the first and the last observation is returned. The user can use the argument \emph{dir} #
#' to prompt the function to focus at previous (\emph{bwd}) or following (\emph{fwd}) observations or to 
#' look in both directions (\emph{both}). The output consists of a list containing two shapefiles, one 
#' \emph{SpatialPointsDataFrame} (\emph{$endpoints}) with the endpoints of each segment and a 
#' \emph{SpatialLinesDataFrame} (\emph{$segments}) representing each segment. The data frames provided 
#' with these shapefiles reports on:
#' \itemize{
#'  \item{\emph{x} - mean x coordinates}
#'  \item{\emph{y} - mean y coordinates}
#'  \item{\emph{timestamp} - mean observation time}
#'  \item{\emph{pixel.time} - elapsed time within a pixel for a given segment}
#'  \item{\emph{travel.distance} - distance between a segment endpoints}
#'  \item{\emph{travel.time} - elapsed time between a segment endpoints}
#'  \item{\emph{stat}: statistical metric}}
#' If \emph{edata} is provided, \emph{img} will only be used as a reference grid as \emph{edata} will 
#' contain the environmental data. Otherwise, this will be retrieved from \emph{img}.}
#' @note "xy" should be provided with a cartesian coordinate system (e.g. UTM).
#' @examples {
#'  
#'  require(raster)
#'  
#'  # read movement data
#'  file <- system.file('extdata', 'konstanz_20130804.shp', package="rsMove")
#'  moveData <- shapefile(file)
#'  
#'  # observation time
#'  td <- as.Date(moveData@data$date)
#'  
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'  
#'  # perform directional sampling
#'  of <- function(x,y) {lm(y~x)$coefficients[2]}
#'  s.sample <- spaceDir(xy=moveData, ot=td, img=r, dir="bwd", fun=of)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

spaceDir <- function(xy=xy, ot=ot, img=img, edata=NULL, dir=dir, fun=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy)
  
  # sample dates
  if (!exists('ot')) {stop('"ot" is missing')}
  if (!class(ot)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"ot" is nof of a valid class')}
  if (length(ot)!=length(xy)) {stop('errorr: "xy" and "ot" have different lengths')}
  
  # raster
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]=='RasterLayer') {stop('"img" is not of a valid class')}
  if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}   
  
  # query direction
  if (!is.null(dir)) {
    if (length(dir)>1) {stop('"dir" has too many entries')}
    if (!dir%in%c('bwd', 'fwd', 'both')) {stop('"dir" is not a valid entry')}
  } else {dir <- 'both'}
  
  # environmental data
  if (!is.null(edata)) {
    if (class(edata)[1]!='vector') {stop('"edata" provided but not a vector')}
    if (length(edata)!=length(xy)) {stop('number of elements in "xy" and "edata" do not match')}}
  
  # check/define input metrics
  if (is.null(fun)) {fun <- function(x,y) {lm(y~x)$coefficients[2]}} else {
    if(!is.function(fun)) {stop('"fun" is not a valid function')}}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 2. summarize unique pixel segments
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # convert xy to single pixels
  ext <- extent(img)
  pxr <- res(img)
  nr <- dim(img)[1]
  sp <- sp <- ((ext[4]-xy@coords[,2])%/%pxr) + nr * (((xy@coords[,1]-ext[1])%/%pxr)-1)
  
  # search for segments and return median values
  sp0 <- 1
  li <- 1
  ux <- list() # x coordinates
  uy <- list() # y coordinates
  ut <- list() # observation time
  et <- list() # elapsed time
  sg <- list() # segment position
  for (r in 2:length(sp)) {
    if (sp[r]!=sp[r-1]) {
      ep <- (r-1)
      ux[[li]] <- mean(xy@coords[sp0:ep,1])
      uy[[li]] <- mean(xy@coords[sp0:ep,2])
      ut[[li]] <- mean(ot[sp0:ep])
      et[[li]] <- difftime(ot[ep], ot[sp0], units='mins')
      sg[[li]] <- li
      sp0 <- r
      li <- li + 1
    }
  }
  
  # convert to vector
  ux <- unlist(ux)
  uy <- unlist(uy)
  ut <- do.call("c", ut)
  et <- unlist(et)
  sg <- unlist(sg)
  
  rm(ext, sp, sp0, li, ot)

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. select pixels between consecutive points
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # backward sampling
  if (dir=='bwd') {
    f1 <- function(i) {
      si <- i-1
      ei <- i
      d0 <- sqrt((ux[si]-ux[ei])^2 + (uy[si]-uy[ei])^2)
      t0 <- difftime(ut[si], ut[ei], units='mins')
      x0 <- ux[si:ei]
      y0 <- uy[si:ei]
      if((d0 > (pxr[1]*2))) {
        dx <- x0[1]-x0[2]
        dy <- y0[1]-y0[2]
        if (abs(dx)>abs(dy)) {
          m <- lm(y0~x0)$coefficients
          if (dx>0) {cm<- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dx)/pxr[1]))*pxr[1])}
          x0 <- (x0[1]+cm)
          y0 <- m[1] + x0 * m[2]
        } else {
          m <- lm(x0~y0)$coefficients
          if (dy>0) {cm<- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dy)/pxr[1]))*pxr[1])}
          y0 <- (y0[1]+cm)
          x0 <- m[1] + y0 * m[2]}}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), sg[i])))}
    op <- lapply(2:length(sg), f1)}
  
  # forward sampling
  if (dir=='fwd') {
    f1 <- function(i) {
      si <- i
      ei <- i+1
      d0 <- sqrt((ux[si]-ux[ei])^2 + (uy[si]-uy[ei])^2)
      t0 <- difftime(ut[si], ut[ei], units='mins')
      x0 <- ux[si:ei]
      y0 <- uy[si:ei]
      if((d0 > (pxr[1]*2))) {
        dx <- x0[1]-x0[2]
        dy <- y0[1]-y0[2]
        if (abs(dx)>abs(dy)) {
          m <- lm(y0~x0)$coefficients
          if (dx>0) {cm<- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dx)/pxr[1]))*pxr[1])}
          x0 <- (x0[1]+cm)
          y0 <- m[1] + x0 * m[2]
        } else {
          m <- lm(x0~y0)$coefficients
          if (dy>0) {cm<- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dy)/pxr[1]))*pxr[1])}
          y0 <- (y0[1]+cm)
          x0 <- m[1] + y0 * m[2]}}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), sg[i])))}
    op <- lapply(1:(length(sg)-1), f1)}
  
  # backward-forward sampling
  if (dir=='both') {
    f1 <- function(i) {
      si <- i-1
      ei <- i+1
      d0 <- sqrt((ux[si]-ux[ei])^2 + (uy[si]-uy[ei])^2)
      t0 <- difftime(ut[si], ut[ei], units='mins')
      x0 <- ux[si:ei]
      y0 <- uy[si:ei]
      if((d0 > (pxr[1]*2))) {
        dx <- x0[1]-x0[2]
        dy <- y0[1]-y0[2]
        if (abs(dx)>abs(dy)) {
          m <- lm(y0~x0)$coefficients
          if (dx>0) {cm<- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dx)/pxr[1]))*pxr[1])}
          x0 <- (x0[1]+cm)
          y0 <- m[1] + x0 * m[2]
        } else {
          m <- lm(x0~y0)$coefficients
          if (dy>0) {cm<- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dy)/pxr[1]))*pxr[1])}
          y0 <- (y0[1]+cm)
          x0 <- m[1] + y0 * m[2]}}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), sg[i])))}
    op <- lapply(2:(length(sg)-1), f1)}
    
  # retrieve x, y and p
  xc <- unlist(sapply(op, function(i){i$x}))
  yc <- unlist(sapply(op, function(i){i$y}))
  pd <- unlist(sapply(op, function(i){i$d}))
  td <- unlist(sapply(op, function(i){i$t}))
  us <- unlist(sapply(op, function(i){i$p}))
    
  rm(op)
    
#-------------------------------------------------------------------------------------------------------------------------------#
# 4. query samples
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # retrieve environmental variables
  if (class(edata)[1]!='data.frame') {edata <- extract(img, cbind(xc,yc))}
  
  # query function
  f2 <- function(i) {
    ind <- us==i
    u <- !is.na(edata[ind])
    if (sum(u)>1) {return(as.numeric(fun(1:length(u),as.numeric(edata[ind[u]]))))} else {return(NA)}}
  
  # apply user provided functon
  if(dir=='bwd') {ov <- sapply(2:length(sg), f2)}
  if(dir=='fwd') {ov <- sapply(1:(length(sg)-1), f2)}
  if(dir=='both') {ov <- sapply(2:(length(sg)-1), f2)}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 5. derive output
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # build dara table
  ind <- which(sg %in%us) # used positions
  df <- data.frame(x=ux[ind], y=uy[ind], timestamp=ut[ind], pixel.time=et[ind], travel.distance=pd, travel.time=td, stat=ov, sid=sg[ind])
  
  # build segment endpoint shapefile
  p.shp <- SpatialPointsDataFrame(df[,1:2], df, proj4string=rProj)
  
  # build segment path shapefile
  f <- function(x) {
    loc <- which(us==sg[ind[x]])
    return(Lines(list(Line(cbind(xc[loc],yc[loc]))), x))}
  l.shp = SpatialLinesDataFrame(SpatialLines(lapply(1:length(ind), f), proj4string=rProj), df)
  
  # output
  return(list(endpoints=p.shp, segments=l.shp))
  
}