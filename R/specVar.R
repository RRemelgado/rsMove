#' @title specVar
#'
#' @description {Tool to support the selection of adequate satellite spatial resoltuon. Evaluates
#' how the spectral variability within a pixel change with the change in spatial resolution.}
#' @param img Object of class \emph{RasterLayer}.
#' @param pixel.res Spatial resolution (unit depends on the spatial projection).
#' @importFrom raster extent xyFromCell crs aggregate crop disaggregate getValues setValues res extract
#' @importFrom ggplot2 ggplot aes_string geom_boxplot ylim
#' @return A \emph{list}.
#' @details {Given a raster object (\emph{img}), the function determines how degrading its spatial resolution impacts
#' our ability to perceive the complexity of the landscape. For the pixel resolution given by \emph{pixel.res}, The
#' function resamples \emph{img} and estimates the Mean Absolute Percentage Error (MAPE) for each pixel. The MAPE is
#' estimated as \eqn{100 / n * sum(abs(O - A / O)} where \emph{O} are the original value in \emph{img}, \emph{A} the
#' aggregated value in the aggregated image and \emph{n} the number of non-NA pixels in the original image The
#' output of the function consists of:
#' \itemize{
#'  \item{\emph{mape} - MAPE raster.}
#'  \item{\emph{plot} - Histogram of \emph{mape}.}}}
#' @seealso \code{\link{tMoveRes}} \code{\link{sMoveRes}}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'
#'  # apply function
#'  s.var <- specVar(img=r, pixel.res=60)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

specVar <- function(img=img, xy=NULL, pixel.res=pixel.res) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  # raster
  if (!class(img)[1]=='RasterLayer') {stop('"img" is not of a valid class')}

  # pixel resolution
  if (!is.numeric(pixel.res)) {stop('"pixel.res" is not numeric')}
  if (length(pixel.res) > 1) {stop('"pixel.res" is not a vector')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. extract raster parameters and define scaling factor
#---------------------------------------------------------------------------------------------------------------------#

  # evaluate each resolution
  ext <- extent(img) # reference extent
  rproj <- crs(img) # reference projection
  ixy <- xyFromCell(img, 1:ncell(img)) # original xy
  orv <- getValues(img)

  # aggregation factor
  af <- cbind(pixel.res / res(img)[1], pixel.res / res(img)[2])
  cc <- (as.integer(af)) == af
  if (min(cc)==0) {warning('one or more elements in "pixel.res" is not a multiple of "img" resolution')}

#---------------------------------------------------------------------------------------------------------------------#
# 3. derive Mean Absolute Percentage Error (MAPE)
#---------------------------------------------------------------------------------------------------------------------#

  # MAPE function
  mape <- function(x) {
    ind <- !is.na(x)
    if (length(ind) > 1) {100 / length(x[ind]) * abs(sum(x[ind]) / length(x[ind]))}}

  # resample data to and from a higher resolution
  tmp <- aggregate(img, fact=af, fun=mean, na.rm=TRUE)

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
  rm(img, ixy, orv, sp)

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
  p <- ggplot(out, aes_string('\nMAPE')) + theme_bw() +
    geom_histogram(binwidth=1) +
    xlim(0, 100) + ylab('Pixel Frequency\n')

  #---------------------------------------------------------------------------------------------------------------------#
  # 5. return output
  #---------------------------------------------------------------------------------------------------------------------#

  # return data frame and plot
  return(list(mape=tmp, plot=p))

}
