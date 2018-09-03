#' @title specVar
#'
#' @description {Tool to support the selection of adequate satellite spatial resolution. Evaluates
#' how the spectral variability within a pixel change with the change in spatial resolution.}
#' @param x Object of class \emph{RasterLayer}.
#' @param y Spatial resolution (unit depends on the spatial projection).
#' @importFrom raster extent xyFromCell crs aggregate crop disaggregate getValues setValues res extract
#' @importFrom ggplot2 ggplot aes_string geom_boxplot ylim geom_histogram xlim
#' @return A \emph{list}.
#' @details {Given a raster object (\emph{x}), the function determines how degrading its spatial resolution
#' impacts our ability to perceive the complexity of the landscape. For the pixel resolution given by \emph{y},
#' The function resamples \emph{x} and estimates the Mean Absolute Percentage Error (MAPE) for each pixel. The
#' MAPE is estimated as \eqn{100 / n * sum(abs(O - A / O)} where \emph{O} are the original value in \emph{x},
#' \emph{A} the aggregated value in the aggregated image and \emph{n} the number of non-NA pixels in the original
#' image The output of the function consists of:
#' \itemize{
#'  \item{\emph{mape} - MAPE raster.}
#'  \item{\emph{plot} - Histogram of \emph{mape}.}}}
#' @seealso \code{\link{tMoveRes}} \code{\link{sMoveRes}}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  # apply function
#'  s.var <- specVar(r, 60)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

specVar <- function(x, y) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  # raster
  if (!class(x)[1]=='RasterLayer') {stop('"x" is not of a valid class')}

  # pixel resolution
  if (!is.numeric(y)) {stop('"y" is not numeric')}
  if (length(y) > 1) {stop('"y" is not a vector')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. extract raster parameters and define scaling factor
#---------------------------------------------------------------------------------------------------------------------#

  # evaluate each resolution
  ext <- extent(x) # reference extent
  rproj <- crs(x) # reference projection
  ixy <- xyFromCell(x, 1:ncell(x)) # original xy
  orv <- getValues(x)

  # aggregation factor
  af <- cbind(y / res(x)[1], y / res(x)[2])
  cc <- (as.integer(af)) == af
  if (min(cc)==0) {warning('one or more elements in "y" is not a multiple of "x" resolution')}

#---------------------------------------------------------------------------------------------------------------------#
# 3. derive Mean Absolute Percentage Error (MAPE)
#---------------------------------------------------------------------------------------------------------------------#

  # MAPE function
  mape <- function(x) {
    ind <- !is.na(x)
    if (length(ind) > 1) {100 / length(x[ind]) * abs(sum(x[ind]) / length(x[ind]))}}

  # resample data to and from a higher resolution
  tmp <- aggregate(x, fact=af, fun=mean, na.rm=TRUE)

  # convert coordinates to pixel positions
  sp <- cellFromXY(tmp, ixy)

  # derive MAPE for each pixel
  out <- sapply(1:ncell(tmp), function(x) {
    ind <- which(sp==x)
    ind <- ind[which(!is.na(orv[ind]))]
    if (length(ind) > 1) {
      return(100 / length(ind) * sum(abs((orv[ind]-tmp[x])/orv[ind])))} else {return(NA)}})

  # retrieve per-pixel values and determine optimal resolution
  tmp <- setValues(tmp, out)

  # remove temporary data from memory
  rm(x, ixy, orv, sp)

#---------------------------------------------------------------------------------------------------------------------#
# 4. derive plot
#---------------------------------------------------------------------------------------------------------------------#

  # determine y scale range
  mv = max(out, na.rm=TRUE)
  if (mv < 100) {
    mv <- mv / 10
    yr <- round(mv*2)/2
    if (mv > yr) {yr <- (yr+0.5)*10} else {yr = yr*10}}
  if (mv >= 100) {
    mv <- mv / 100
    yr <- round(mv*20)/20
    if (mv > yr) {yr <- (yr+0.5)*100} else {yr <- yr*100}}

  # build plot object
  out <- data.frame(MAPE=out)
  p <- ggplot(out, aes_string(x="MAPE")) + theme_bw() +
    geom_histogram(binwidth=sd(out$MAPE, na.rm=TRUE)) +
    ylab('Relative freq. (%)') + xlab("Value") + ylim(0,100)

  #---------------------------------------------------------------------------------------------------------------------#
  # 5. return output
  #---------------------------------------------------------------------------------------------------------------------#

  # return data frame and plot
  return(list(mape=tmp, plot=p))

}
