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

specVar <- function(img=img, p.res=T) {
  
#---------------------------------------------------------------------------------------------------------------------#
  #  1. check inpur variables
  #---------------------------------------------------------------------------------------------------------------------#
  
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (!is.numeric(pxr)) {stop('"pxr" is not numeric')}
  if (!is.vector(pxr)) {stop('"pxr" is not a vector')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. convert pixels to coordinates and extract raster info
#---------------------------------------------------------------------------------------------------------------------#
    
  # evaluate each resolution
  ext <- extent(img) # reference extent
  rproj <- crs(img) # reference projection
  xy <- xyFromCell(img, 1:ncell(img)) # original xy

#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#
  
  # output variable
  out <- vector('list', length(pxr))
  
  # loop through each resolution
  for (p in 1:length(pxr)) {    
    
    # resample data to and from a higher resolution
    tmp <- resample(resample(img, raster(ext, res=pxr[p], crs=rproj)), img)
    
    # extract difference vales
    rv <- getValues(img - tmp)
    
    # convert coordinates to pixel positions
    nc <- round((ext[2]-ext[1]) / pxr[p]) + 1 # number of columns
    nr <- round((ext[4]-ext[3]) / pxr[p]) + 1 # number of rows
    sp <- (round((ext[4]-xy[,2])/pxr[p])+1) + nr * round((xy[,1]-ext[1])/pxr[p])
    
    # remove NA values
    cc <- complete.cases(rv)
    rv <- rv[cc]
    sp <- sp[cc]
    
    # derive RMSE for each pixel
    up <- unique(sp) # unique pixel positions
    rv <- sapply(up, function(x) {ind <- which(sp==x)
    sqrt((sum(rv[ind]^2) / length(ind)))})
    
    # add output to list
    out[[p]] <- list(value=rv, resolution=replicate(length(rv), pxr[p]))
    
    # remove temporary data from memory
    rm(tmp, nc, nr, sp, cc, rv, up)
    
  }
  
  # final data frame
  odf <- data.frame(RMSE=unlist(sapply(out, function(x){x$value})), 
                    Resolution=factor(unlist(sapply(out, function(x){x$resolution})), levels=as.character(pxr)))
  
  # remove temporary data from memory
  rm(xy, ext, rproj, out)
  
#---------------------------------------------------------------------------------------------------------------------#
# 4. plot output
#---------------------------------------------------------------------------------------------------------------------#
  
  # determine y scale range
  mv = max(odf$RMSE)
  if (mv < 100) {
    mv <- mv / 10
    yr <- round(mv*2)/2
    if (mv > yr) {yr <- (yr+0.5)*10} else {yr = yr*10}}
  if (mv >= 100) {
    mv <- mv / 100
    yr <- round(mv*20)/20
    if (mv > yr) {yr <- (yr+0.5)*100} else {yr <- yr*100}}
  
  # build plot object
  p <- ggplot(odf, aes(x=Resolution, y=RMSE)) + geom_boxplot() + ylim(0,yr)
  
  if (p.res) {p} # plot on screen
  
  # return data frame and plot
  return(list(stats=odf, plot=p))
  
}