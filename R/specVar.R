#' @title specVar
#'
#' @description {Tool to support the selection of adequate satellite spatial resoltuon. Evaluates
#' how the spectral variability within a pixel change with the change in spatial resolution.}
#' @param img Object of class \emph{RasterLayer}.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param pixel.res vector of target resolutions.
#' @import ggplot2 sp raster rgdal grDevices
#' @return A \emph{list}.
#' @details {Given a raster object (\emph{img}), the function determines how degrading its spatial resolution impacts
#' our ability to perceive the complexity of the landscape. For each pixel resolution given by \emph{pixel.res}, The
#' function resamples \emph{img} and estimates the Mean Absolute Error (MAE) for each pixel based on the differences
#' between the original and aggregated values. If a point shapefile is provided (\emph{xy}), the function will only
#' report on the values that overlap with it. The output of the function consists of:
#' \itemize{
#'  \item{\emph{mae} - MAE for all pixels or for the points within \emph{xy}.}
#'  \item{\emph{pixel.optimal} - raster reporting on the resolution with the lowest MAE for each pixel.}
#'  \item{\emph{pixel.optimal.stats} - Proportion of samples per resolution derived from \emph{pixel.optimal}.}
#'  \item{\emph{plot} - boxplots of the variability of the MAE per resolution.}
#'  \item{\emph{sample.optial} - optimal resolution per sample (available if \emph{xy} is provided).}}}
#' @seealso \code{\link{tMoveRes}} \code{\link{sMoveRes}}
#' @examples \dontrun{
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
#'  # apply function
#'  s.var <- specVar(img=r, xy=moveData, pixel.res=60)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

specVar <- function(img=img, xy=NULL, pixel.res=pixel.res) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  # raster/shapefile
  if (!exists('img')) {stop('"img" is missing')}
  if (!class(img)[1]=='RasterLayer') {stop('"img" is not of a valid class')}
  if (!is.null(xy)) {
    if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
    if (crs(xy)@projargs!=crs(img)@projargs) {stop('"xy" and "img" have different projections')}}

  # pixel resolution
  if (!is.numeric(pixel.res)) {stop('"pixel.res" is not numeric')}
  if (!is.vector(pixel.res)) {stop('"pixel.res" is not a vector')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. extract raster parameters and define scaling factor
#---------------------------------------------------------------------------------------------------------------------#

  # evaluate each resolution
  ext <- extent(img) # reference extent
  rproj <- crs(img) # reference projection
  ixy <- xyFromCell(img, 1:ncell(img)) # original xy

  # aggregation factor
  af <- cbind(pixel.res / res(img)[1], pixel.res / res(img)[2])
  cc <- (as.integer(af)) == af
  if (min(cc)==0) {stop('one or more elements in "pixel.res" is not a multiple of "img" resolution')}

#---------------------------------------------------------------------------------------------------------------------#
# 3. extract statistics
#---------------------------------------------------------------------------------------------------------------------#

  # output variable
  out <- vector('list', length(pixel.res))

  # loop through each resolution
  for (p in 1:length(pixel.res)) {

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
    out[[p]] <- list(value=rv, resolution=replicate(length(rv), pixel.res[p]), original=orv)

    # remove temporary data from memory
    rm(tmp, sp, rv, up, orv)

  }

#---------------------------------------------------------------------------------------------------------------------#
# 4. derive output statistics
#---------------------------------------------------------------------------------------------------------------------#

  # retrieve per-pixel values and determine optimal resolution
  opr <- data.frame(lapply(out, function(x){x$original}))
  opr <- setValues(img, apply(opr, 1, function(x) {max(pixel.res[which(x==min(x))])}))

  # derive raster stats
  uv <- unique(opr)
  nc <- ncell(opr)
  cp <- data.frame(Resolution=uv, Percentage=sapply(uv, function(x) {cellStats(opr==x, sum)/nc *100}))

  # build output data frame
  if (!is.null(xy)) {
    odf <- as.data.frame(do.call(cbind, lapply(out, function(x) {x$value})))
    colnames(odf) <- c(as.character(pixel.res))
    # optimal resolution per per sample / pixel
    osr <- as.numeric(apply(odf, 1, function(x) {max(pixel.res[which(x==min(x))])}))

    # sort data frame (ggplot format)
    out <- lapply(1:length(pixel.res), function(x) {
      pr <- factor(replicate(length(odf[,x]), pixel.res[x]),
                   levels=as.character(pixel.res))
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
  p <- ggplot(out, aes_string(x="Resolution", y="MAE")) + geom_boxplot() + ylim(0,yr)

  # return data frame and plot
  if (!is.null(xy)) {
    return(list(mae=odf, sample.optimal=osr, pixel.optimal=opr, pixel.optimal.stats=cp, plot=p))
  } else {return(list(mae=out, pixel.optimal=opr, pixel.optimal.stats=cp, plot=p))}

}
