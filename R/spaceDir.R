#' @title spaceDir
#'
#' @description Spatial directional raster analysis.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param img Object of class \emph{RasterLayer}.
#' @param dir One of \emph{fwd}, \emph{bwd} or \emph{both}. Default is \emph{both}.
#' @param type One of 'cont' or 'cat'. Defines which type of variable is in use.
#' @param dm One of 'm' or 'deg' specifying the projection unit. Default is 'm'.
#' @param fun Summary function.
#' @import raster sp rgdal
#' @importFrom stats lm
#' @seealso \code{\link{timeDir}} \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{list}.
#' @details {This function evaluates how do environmental conditions change in space 
#' along a movement track. Before an output is returned, the function looks for segments 
#' of consecutive samples that appear within the same pixel and provides mean values for 
#' their coordinates and observation times as well as the corresponding pixel value and 
#' the elapsed time during the segment (in minutes). this step avoids the replication of 
#' samples while preserving observation related to revisits to the same pixels. Then, for 
#' each observation, the function looks at it imediate neighbors to define spatial segments. 
#' All the pixels between the two endpoints of the segment are considered in this process. 
#' The user can use the argument \emph{dir} to prompt the function to focus on previous 
#' (\emph{bwd}) or following (\emph{fwd}) observations or to look in both directions (\emph{both}) 
#' when defining the segments. The output consists of a list containing two shapefiles, one 
#' \emph{SpatialPointsDataFrame} (\emph{$endpoints}) with the endpoints of each segment and a 
#' \emph{SpatialLinesDataFrame} (\emph{$segments}) representing each segment. The data frames provided 
#' with these shapefiles reports on:
#' \itemize{
#'  \item{\emph{x} - mean x coordinates}
#'  \item{\emph{y} - mean y coordinates}
#'  \item{\emph{timestamp} - mean observation time}
#'  \item{\emph{pixel.time} - elapsed time within a pixel for a given segment}
#'  \item{\emph{travel.distance} - distance between a segment endpoints}
#'  \item{\emph{travel.time} - elapsed time between a segment endpoints}}
#'  The additionaly information added to the data frame depends on the type of data. If \emph{type} is set 
#'  to \emph{'cont'}, the function will assume the input raster is a continuous variable as estimate a 
#'  statistical metric for each segment. If none is provided throught the keyword \emph{fun}, the slope 
#'  will be returned by default and assigned to the column \emph{'stat'} within the output data frame. If 
#'  \emph{type} is set to '\emph{'cat'}, the function will assume the data is categorical. As a result,
#'  the keyword \emph{fun} is ignored. Instead, the function will return the proportion of pixels occupied 
#'  by each class for each segment. In addition, the dominant class (\emph{'main'}) and the shannon index 
#'  (\emph{'shannon'}) will be returned for each segment within the output data frame.}
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
#'  td <- as.Date(moveData@data$date)
#'  
#'  # perform directional sampling
#'  of <- function(x) {lm(1:lengh(x)~x)$coefficients[2]}
#'  s.sample <- spaceDir(xy=moveData, ot=td, img=r, dir="bwd", type='cont', fun=of)
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

spaceDir <- function(xy=xy, ot=ot, img=img, dir=dir, type=type, dm='m', fun=NULL) {
  
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
  
  # variable type
  if (is.null(type)) {stop('"type" is missing')} else {
    if (!type%in%c('cont', 'cat')) {stop('"type" is not a recognized keyword')}}
  
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
      t0 <- difftime(ut[ei], ut[si], units='mins')
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
      if (dm=='deg') {
        x00 <- ux[si:ei]*pi/180
        y00 <- uy[si:ei]*pi/180
        sf <- function(o) {
          xDiff <- abs(x00[(o-1)]-x00[o])
          yDiff <- abs(y00[(o-1)]-y00[o])
          aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y00[o]) * cos(x00[o]) * sin(xDiff/2.) * sin(xDiff/2.)
          cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
          return(6371000 * cCoef)}
        d0 <- sum(sapply(2:length(x00), sf))}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), sg[i])))}
    op <- lapply(2:length(sg), f1)}
  
  # forward sampling
  if (dir=='fwd') {
    f1 <- function(i) {
      si <- i
      ei <- i+1
      d0 <- sqrt((ux[si]-ux[ei])^2 + (uy[si]-uy[ei])^2)
      t0 <- difftime(ut[ei], ut[si], units='mins')
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
      if (dm=='deg') {
        x00 <- ux[si:ei]*pi/180
        y00 <- uy[si:ei]*pi/180
        sf <- function(o) {
          xDiff <- abs(x00[(o-1)]-x00[o])
          yDiff <- abs(y00[(o-1)]-y00[o])
          aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y00[o]) * cos(x00[o]) * sin(xDiff/2.) * sin(xDiff/2.)
          cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
          return(6371000 * cCoef)}
        d0 <- sum(sapply(2:length(x00), sf))}
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
        x00 <- vector('list', 2)
        y00 <- vector('list', 2)
        for (o in 2:3) {
          dx <- x0[(o-1)]-x0[o]
          dy <- y0[(o-1)]-y0[o]
          if (abs(dx)>abs(dy)) {
            m <- lm(y0~x0)$coefficients
            if (dx>0) {cm <- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm <- ((0:round(abs(dx)/pxr[1]))*pxr[1])}
            tx <- (x0[(o-1)]+cm)
            ty <- m[1] + x0[(o-1):o] * m[2]
          } else {
            m <- lm(x0~y0)$coefficients
            if (dy>0) {cm <- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm <- ((0:round(abs(dy)/pxr[1]))*pxr[1])}
            y0 <- (y0[(o-1)]+cm)
            x0 <- m[1] + y0[(o-1):o] * m[2]}
            if (o==3) {
              x00[(o-1)] <- tx[2:length(tx)]
              y00[(o-1)] <- ty[2:length(ty)]
            } else {
              x00[(o-1)] <- tx
              y00[(o-1)] <- ty}}
        x0 <- unlist(x00)
        y0 <- unlist(y00)}
      if (dm=='deg') {
        x00 <- ux[si:ei]*pi/180
        y00 <- uy[si:ei]*pi/180
        sf <- function(o) {
          xDiff <- abs(x00[(o-1)]-x00[o])
          yDiff <- abs(y00[(o-1)]-y00[o])
          aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y00[o]) * cos(x00[o]) * sin(xDiff/2.) * sin(xDiff/2.)
          cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
          return(6371000 * cCoef)}
        d0 <- sum(sapply(1:length(x00), sf))}
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
  edata <- extract(img, cbind(xc,yc))
  
  
  if (type=='cont') {
    
    # query function
    f2 <- function(i) {
      ind <- which(us==i)
      u <- which(!is.na(edata[ind]))
      if (sum(u)>1) {return(as.numeric(fun(as.numeric(edata[ind[u]]))))} else {return(NA)}}
    
    # apply user provided functon
    if(dir=='bwd') {ov <- sapply(2:length(sg), f2)}
    if(dir=='fwd') {ov <- sapply(1:(length(sg)-1), f2)}
    if(dir=='both') {ov <- sapply(2:(length(sg)-1), f2)}
    ov <- data.frame(stat=ov)
  
  }
  
  if (type=='cat') {
    
    uc <- unique(img) # unique classes
    
    # function to sumarize class composition
    f2 <- function(i) {
      ind <- which(us==i)
      return(sapply(uc, function(y) {sum(edata[ind]==uc[y], na.rm=T)}))}
    
    if(dir=='bwd') {ov <- do.call(rbind, lapply(2:length(sg), f2))}
    if(dir=='fwd') {ov <- do.call(rbind, lapply(1:(length(sg)-1), f2))}
    if(dir=='both') {ov <- do.call(rbind, lapply(2:(length(sg)-1), f2))}
    
    # derive additional statistics
    dc <- apply(ov, 1, function(x) {uc[which(x==max(x))]})
    si <- apply(ov, 1, function(x) {-sum((x/sum(x))*log(x/sum(x)))})
    
    # update output table
    ov <- data.frame(ov, dc, si, stringsAsFactors=F)
    colnames(ov) <- c(as.character(uc), 'main', 'shannon')
    
  }
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 5. derive output
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # build dara table
  ind <- which(sg %in%us) # used positions
  df <- data.frame(x=ux[ind], y=uy[ind], timestamp=ut[ind], pixel.time=et[ind], travel.distance=pd, travel.time=td, sid=sg[ind])
  df <- cbind(df, ov)
  
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