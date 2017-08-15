#' @title moveSeg
#'
#' @description Remote sensing based point segmentation
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param edata Object of class \emph{RasterLayer} or \emph{data.frame}.
#' @param type Raster data type. One of \emph{cont} (continues) or \emph{cat} (for categorical).
#' @param ot Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates. 
#' @param threshold Change threshold.
#' @param b.size Buffer size expressed in the map units.
#' @param s.fun Output summary function. Default is mean.
#' @import raster rgdal ggplot2
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{list}.
#' @details {Segmentation of a point shapefile based on the spatial variability 
#' of a raster dataset. When the \emph{method} is set to \emph{'cont'}, the raster 
#' data is assumed to be continuous. Then, the function determines the percentual 
#' change between each pair of two consecutive coordinate pairs. If this change is 
#' above a predifined \emph{threshold}, a new pointer is added and the previous 
#' sequence of samples is labeled as a unique segment. If \emph{method} is set as 
#' \emph{'cat'}, the function assumes the raster data is categorical ignoring the 
#' \emph{theshold} and \emph{s.fun} keywords. In this case, a new segment is identified 
#' if the any change is observed between two consecutife points. The output consists of 
#' a list containing a \emph{SpatialPointsDataFrame} (\emph{$points}) reporting on the 
#' segment ID (\emph{sid}) associated to each sample and a data frame (\emph{$report}) 
#' with the amount of points in each region and the value returned by \emph{s.fun}. If 
#' \emph{ot} is provided, the function also provides the elapsed time within each segment. 
#' If \emph{fun} is set by the user, the provided function will be used to summarize the 
#' raster values at each segment. Also, if \emph{edata} is a \emph{RasterStack} or a 
#' \emph{RasterBrick}, \emph{r.fun} is used to reduce the multi-layered object to a single 
#' layer. he user can either use a Principal Component Analysis (PCA) analysis by setting 
#' \emph{r.fun} to \emph{pca} or provide a function. If \emph{pca} is selected, the first 
#' Principal Conponent (PC) \emph{threshold} will be set automatically to the standard 
#' deviation of the first PC. This is the default. If \emph{method} is \emph{cat} the 
#' function will assume the data is categorical and wil define segments whenever a change occurres.}
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
#'  seg <- moveSeg(xy=moveData, ot=o.time, edata=r, type="cont", threshold=0.1)
#'  
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------#

moveSeg <- function(xy=xy, edata=edata, type='cont', ot=NULL, b.size=NULL, threshold=NULL, s.fun=NULL) {
  
#---------------------------------------------------------------------------------------------------------------------#
# 1. check input variables  
#---------------------------------------------------------------------------------------------------------------------#
  
  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy) # output projection
  
  # sample dates
  if (!is.null(ot)) {
    if (!class(ot)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"ot" is nof of a valid class')}
    if (length(ot)!=length(xy)) {stop('errorr: "xy" and "ot" have different lengths')}}
  
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
    if (!is.null(threshold)) {stop('"type" is set to "cat". Please define "threshold"')}
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
      if (type=='cont') {edata <- extract(img, xy@coords, buffer=b.size, fun=s.fun, na.rm=T)}
      
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
        edata <- extract(img, tmp)
        
        # sumarize data (extract dominant class)
        edata <- sapply(1:length(xy), function(x) {
          ind <- which(si==x)
          r0 <- as.vector(edata[ind[!duplicated(cellFromXY(img, tmp[ind,]))]])
          uc <- unique(r0)
          uc <- uc[!is.na(uc)]
          if (length(uc)>0) {
            count <- sapply(uc, function(x) {sum(r0==x)})
            return(uc[which(count==max(count))[1]])
          } else {return(NA)}})
          
          rm(tmp, si)
        
      }
      
      # simple query
    } else {edata <- extract(img, xy@coords)}}
  
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
  if (!is.null(ot)) {
    f <- function(x) {
      ind <- which(id==x)
      et <- difftime(ot[ind[length(ind)]], ot[ind[1]], units="mins")
      np <- length(ind)
      return(list(time=as.numeric(et), count=np))}
    sstat <- lapply(uid, f)
    df <- data.frame(sid=uid, count=sapply(sstat, function(x) {x$count}), 
          time=sapply(sstat, function(x) {x$time}), value=rv)
  } else {df <- data.frame(sid=uid, count=sapply(uid, function(x){sum(id==x)}))}
  
  # estimate relative time
  df$rel.time <- time / diffTime(max(ot), min(ot), units='mins') * 100
  
#---------------------------------------------------------------------------------------------------------------------#
# 5. build plot
#---------------------------------------------------------------------------------------------------------------------#
  
  if (type=='cont') {
    
    # plot with time
    if (!is.null(ot)) {
      
      # buid plot object
      p <- ggplot(df, aes(x=as.factor(sid), y=time, fill=as.factor(value))) + geom_bar(stat="identity") +
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
      p <- ggplot(df, aes(x=as.factor(sid), y=time, fill=count)) + geom_bar(stat="identity") +
        xlab('') + ylab('Segment length') +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12), 
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}
  
  if (type=='cat') {
    
    # plot with time
    if (!is.null(ot)) {
      
      # buid plot object
      p <- ggplot(df, aes(x=as.factor(sid), y=time, fill=as.factor(value))) + geom_bar(stat="identity") +
        xlab('') + ylab('Time Spent (minutes)') + scale_fill_discrete(name="Class") + 
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12), 
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))
    
    # plot without time  
    } else {
      
      # buid plot object
      p <- ggplot(df, aes(x=as.factor(sid), y=time, fill=count)) + geom_bar(stat="identity") +
        xlab('') + ylab('Segment length') + scale_fill_discrete(name="Class") +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12), 
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}
  
#---------------------------------------------------------------------------------------------------------------------#
# 6. return output
#---------------------------------------------------------------------------------------------------------------------#
  
  return(list(points=p.shp, report=df, plot=p))
  
}