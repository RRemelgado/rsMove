#' @title moveSeg
#'
#' @description Remote sensing based point segmentation
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param edata Object of class \emph{RasterLayer} or \emph{data.frame}.
#' @param type Raster data type. One of \emph{cont} (continues) or \emph{cat} (for categorical).
#' @param o.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param threshold Change threshold.
#' @param b.size Buffer size expressed in the map units.
#' @param s.fun Output summary function. Default is mean.
#' @import raster rgdal ggplot2
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{list}.
#' @details {Segmentation of a point shapefile based on the spatial variability
#' of a raster dataset. When the \emph{type} is set to \emph{'cont'}, the \emph{edata} 
#' is assumed to be continuous. Then, the function determines the percentual
#' change between each pair of two consecutive coordinate pairs. If this change is
#' above a predifined \emph{threshold}, a new pointer is added and the previous
#' sequence of samples is labeled as a unique segment.
#' If \emph{method} is set as \emph{'cont'}, the function assumes the raster data is a 
#' continuous variable andwill require the user to define \emph{theshold} which indicates 
#' when the difference between consecutive points should be considered a change. If \emph{b.size} 
#' is set and \emph{edata} is a raster, the function will use a buffer when extracting the values for 
#' \emph{xy} and use \emph{s.fun} to summarize the values for each point.
#' If \emph{type} is \emph{'cat'}, the \emph{edata} will be assumed to be categorical in which case 
#' \emph{threshold} and {s.fun} will be ignored. If \emph{type} is 'cat' and \emph{b.size} is provided 
#' the function will return the dominant value for each point. 
#' The output consists of a list containing a \emph{SpatialPointsDataFrame} (\emph{$points}) reporting 
#' on thesegment ID (\emph{sid}) associated to each sample and a data frame (\emph{$report}) with the 
#' amount of points in each region and the \emph{edata} value. If \emph{o.time} is provided, 
#' the function also estimate the elapsed time within each segment.}
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
#'  # perform directional sampling
#'  seg <- moveSeg(xy=moveData, o.time=o.time, edata=r, type="cont", threshold=0.1)
#'
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------#

moveSeg <- function(xy=xy, edata=edata, type='cont', o.time=NULL, b.size=NULL, threshold=NULL, s.fun=NULL) {
  
  #---------------------------------------------------------------------------------------------------------------------#
  # 1. check input variables
  #---------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy) # output projection
  
  # sample dates
  if (!is.null(o.time)) {
    if (!class(o.time)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"o.time" is nof of a valid class')}
    if (length(o.time)!=length(xy)) {stop('errorr: "xy" and "o.time" have different lengths')}}
  
  # environmental data
  if (class(edata)[1]=='RasterLayer') {
    if (crs(xy)@projargs!=crs(edata)@projargs) {stop('"xy" and "edata" have different projections')}
    prd <- TRUE
  } else {
    if (!class(edata)[1]%in%c('data.frame')) {stop('"edata" is neither a raster or a data frame')}
    if (nrow(edata)!=length(xy)) {stop('number of elements in "xy" and "edata" do not match')}
    prd=FALSE}
  
  # check threshold
  if (type=='cont') {
    if (is.null(threshold)) {stop('"type" is set to "cont". Please define "threshold"')}
    if (!is.numeric(threshold)) {stop('"threshold" is not numeric')}}
  if (type=='cat') {threshold <- 1}
  
  # check query type
  if (!type%in%c('cont', 'cat')) {stop('"type" is not  avalid keyword')}
  if (type=='cont') {
    if (!is.null(s.fun)) {if (!is.function(s.fun)){stop('"s.fun" is not a valid keyword or function')}}
    if (is.null(s.fun)) {s.fun <- function(x) {return(mean(x, na.rm=T))}}}
  
  #---------------------------------------------------------------------------------------------------------------------#
  # 2. query data
  #---------------------------------------------------------------------------------------------------------------------#
  
  if (prd) {
    
    # apply buffer if required
    if (!is.null(b.size)) {
      
      # average samples within buffer
      if (type=='cont') {edata <- extract(edata, xy@coords, buffer=b.size, fun=s.fun, na.rm=T)}
      
      # determine main class within the buffer
      if (type=='cat') {
        
        # dilate samples
        tmp <- lapply(1:length(xy), function(x) {
          ind <- raster(extent((xy@coords[x,1]-b.size), (xy@coords[x,1]+b.size),
                               (xy@coords[x,2]-b.size), (xy@coords[x,2]+b.size)), crs=rProj)
          ind <- xyFromCell(ind, 1:ncell(ind))
          return(list(c=ind, s=replicate(nrow(ind), x)))})
        si <- unlist(lapply(tmp, function(x) {x$s}))
        tmp <- do.call(rbind, lapply(tmp, function(x) {x$c}))
        
        # extract values
        edata0 <- extract(edata, tmp)
        
        # sumarize data (extract dominant class)
        edata <- sapply(1:length(xy), function(x) {
          ind <- which(si==x)
          r0 <- as.vector(edata0[ind[!duplicated(cellFromXY(edata, tmp[ind,1:2]))]])
          uc <- unique(r0)
          uc <- uc[!is.na(uc)]
          if (length(uc)>0) {
            count <- sapply(uc, function(x) {sum(r0==x)})
            return(uc[which(count==max(count))[1]])
          } else {return(NA)}})
        
        rm(tmp, si)
        
      }
      
      # simple query
    } else {edata <- extract(edata, xy@coords)}}
  
  #---------------------------------------------------------------------------------------------------------------------#
  # 3. identify segments
  #---------------------------------------------------------------------------------------------------------------------#
  
  # search for segments
  r0 <- 1
  li <- 1
  id <- list() # segment id
  rv <- list() # segment value
  
  for (r in 2:length(xy)) {
    diff <- abs(edata[r]-edata[(r-1)])
    if (!is.na(diff)) {
      if (diff >= threshold) {
        ep <- r-1
        rv[[li]] <- mean(edata[c(r0:ep)])
        id[[li]] <- replicate(length(c(r0:ep)), li)
        r0 <- r
        li <- li + 1
        if (r==length(xy)) {
          id[[li]] <- li
          rv[[li]] <- edata[r]}
      } else {if (r==length(xy)) {
        ep <- r
        rv[[li]] <- mean(edata[c(r0:ep)])
        id[[li]] <- replicate(length(c(r0:ep)), li)}}}}
  rv <- unlist(rv)
  id <- unlist(id)
  
  #---------------------------------------------------------------------------------------------------------------------#
  # 4. derive statistics
  #---------------------------------------------------------------------------------------------------------------------#
  
  # update original shapefile
  p.shp <- SpatialPointsDataFrame(xy@coords, data.frame(sid=id), proj4string=rProj)
  
  # build region report
  uid <- sort(unique(id))
  if (!is.null(o.time)) {
    f <- function(x) {
      ind <- which(id==x)
      et <- difftime(o.time[ind[length(ind)]], o.time[ind[1]], units="mins")
      np <- length(ind)
      return(list(time=as.numeric(et), count=np))}
    sstat <- lapply(uid, f)
    df <- data.frame(sid=uid, count=sapply(sstat, function(x) {x$count}),
                     time=sapply(sstat, function(x) {x$time}), value=rv)
  } else {df <- data.frame(sid=uid, count=sapply(uid, function(x){sum(id==x)}))}
  
  
  #---------------------------------------------------------------------------------------------------------------------#
  # 5. build plot
  #---------------------------------------------------------------------------------------------------------------------#
  
  if (type=='cont') {
    
    df$sid <- factor(df$sid, levels=unique(df$sid))
    
    # plot with time
    if (!is.null(o.time)) {
      
      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="time", fill="value")) + geom_bar(stat="identity") +
        xlab('') + ylab('Time Spent (minutes)') +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))
      
      # plot without time
    } else {
      
      # determine y range scale range
      mv <- max(df$count)
      nc <- nchar(as.character(mv))
      m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
      mv <- mv / m
      yr <- round(mv)
      if (mv > yr) {yr <- (yr+0.2)*m} else {yr <- yr*m}
      
      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="time", fill="count")) + geom_bar(stat="identity") +
        xlab('') + ylab('Segment length') +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}
  
  if (type=='cat') {
    
    df$sid <- factor(df$sid, levels=unique(df$sid))
    df$value <- factor(df$value, levels=unique(df$value))
    
    # plot with time
    if (!is.null(o.time)) {
      
      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="time", fill="value")) + geom_bar(stat="identity") +
        xlab('') + ylab('Time Spent (minutes)') + scale_fill_discrete(name="Class") +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))
      
      # plot without time
    } else {
      
      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="time", fill="count")) + geom_bar(stat="identity") +
        xlab('') + ylab('Segment length') + scale_fill_discrete(name="Class") +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}
  
  #---------------------------------------------------------------------------------------------------------------------#
  # 6. return output
  #---------------------------------------------------------------------------------------------------------------------#
  
  return(list(points=p.shp, report=df, plot=p))
  
}
