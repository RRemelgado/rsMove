#' @title specVar
#'
#' @description Evaluates how a change in raster resolution impacst the availability of data points.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param pxr vector of target resolutions.
#' @param p.res Should the output be ploted on screen? Default is TRUE.
#' @import ggplot2 sp rgdal grDevices
#' @return A \emph{list}.
#' @details {Given a vector of pixel resolutions (\emph{pxr}), the function determines 
#' the number of unique pixels and unique pixel groups. Additionaly, for each pixel, 
#' the function returns the corresponding pixel indices per resolution showing which 
#' samples would be grouped. The function returns a data frame (\emph{$stats}) and a 
#' plot (\emph{$plot}) with the statistics per resolution as well as a data frame with 
#' the pixel indices per resolution (\emph{$indices}).}
#' @seealso \code{\link{tMoveRes}}
#' @examples {
#'  
#'  require(raster)
#'  
#'  # read movement data
#'  file <- system.file('extdata', 'konstanz_20130804.shp', package="rsMove")
#'  moveData <- shapefile(file)
#'  
#'  # test function for 5, 10 20 and 30 m
#'  a.res <- sMoveRes(xy=moveData, pxr=c(5, 10, 20, 30))
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

specVar <- function(img=img, xy=NULL, p.res=T) {
  
#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#
  
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (!is.numeric(pxr)) {stop('"pxr" is not numeric')}
  if (!is.vector(pxr)) {stop('"pxr" is not a vector')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. extract raster parameters and define scaling factor
#---------------------------------------------------------------------------------------------------------------------#
    
  # evaluate each resolution
  ext <- extent(img) # reference extent
  rproj <- crs(img) # reference projection
  ixy <- xyFromCell(img, 1:ncell(img)) # original xy
  
  # aggregation factor
  af <- cbind(pxr / res(img)[1], pxr / res(img)[2])
  cc <- (as.integer(af)) == af
  if (min(cc)==0) {stop('one or more elements in "pxr" is not a multiple of "img" resolution')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 3. extract statistics
#---------------------------------------------------------------------------------------------------------------------#
  
  # output variable
  out <- vector('list', length(pxr))
  
  # loop through each resolution
  for (p in 1:length(pxr)) {    
    
    # resample data to and from a higher resolution
    tmp <- aggregate(img, fact=af[p,], fun=mean, na.rm=T)
    rv <- crop(disaggregate(tmp, fact=af[p,]), img)
    
    # extract difference vales
    rv <- getValues(img - rv)
      
    # convert coordinates to pixel positions
    sp <- cellFromXY(tmp, ixy)
    
    # derive RMSE for each pixel
    up <- unique(sp) # unique pixel positions
    rv <- sapply(up, function(x) {ind <- which(sp==x & !is.na(x))
    if (length(ind)>0) {1/sum(abs(rv[ind]))*sum(abs(rv[ind])*abs(rv[ind]))} else {return(NA)}})
    
    # assign values to original raster
    orv <- sapply(sp, function(x) {rv[which(up==x)]})
    
    # if xy is provided, extract values
    if (!is.null(xy)) {
      rv <- setValues(tmp, rv)
      rv <- extract(rv, xy@coords)}
    
    # add output to list
    out[[p]] <- list(value=rv, resolution=replicate(length(rv), pxr[p]), original=orv)
    
    # remove temporary data from memory
    rm(tmp, sp, rv, up, orv)
    
  }

#---------------------------------------------------------------------------------------------------------------------#
# 4. derive output statistics
#---------------------------------------------------------------------------------------------------------------------#
  
  # retrieve per-pixel values and determine optimal resolution
  opr <- data.frame(lapply(out, function(x){x$original}))
  opr <- setValues(img, apply(opr, 1, function(x) {max(pxr[which(x==min(x))])}))
  
  # derive raster stats
  uv <- unique(opr)
  nc <- ncell(opr)
  cp <- data.frame(Resolution=uv, Percentage=sapply(uv, function(x) {cellStats(opr==x, sum)/nc *100}))
  
  # build output data frame
  if (!is.null(xy)) {
    odf <- as.data.frame(do.call(cbind, lapply(out, function(x) {x$value})))
    colnames(odf) <- c(as.character(pxr))
    # optimal resolution per per sample / pixel
    osr <- as.numeric(apply(odf, 1, function(x) {max(pxr[which(x==min(x))])}))
    
    # sort data frame (ggplot format)
    out <- lapply(1:length(pxr), function(x) {
      pr <- factor(replicate(length(odf[,x]), pxr[x]), 
                   levels=as.character(pxr))
      return(data.frame(MAE=odf[,x], Resolution=pr))})
    out <- do.call(rbind, out) 
    
  } else {
    
    # sort data frame (ggplot format)
    out <- data.frame(MAE=unlist(lapply(out, function(x) {x$value})), 
                      resolution=unlist(lapply(out, function(x) {x$resolution})))}
  
#---------------------------------------------------------------------------------------------------------------------#
# 5. plot output
#---------------------------------------------------------------------------------------------------------------------#

  # determine y scale range
  mv = max(out$MAE)
  if (mv < 100) {
    mv <- mv / 10
    yr <- round(mv*2)/2
    if (mv > yr) {yr <- (yr+0.5)*10} else {yr = yr*10}}
  if (mv >= 100) {
    mv <- mv / 100
    yr <- round(mv*20)/20
    if (mv > yr) {yr <- (yr+0.5)*100} else {yr <- yr*100}}
  
  # build plot object
  p <- ggplot(out, aes(x=Resolution, y=MAE)) + geom_boxplot() + ylim(0,yr)
  
  if (p.res) {p} # plot on screen
  
  # return data frame and plot
  if (!is.null(xy)) {
    return(list(mae=odf, sample.optimal=osr, pixel.optimal=opr, pixel.optimal.stats=cp, plot=p))
  } else {return(list(mae=out, pixel.optimal=opr, pixel.optimal.stats=cp, plot=p))}
  
}