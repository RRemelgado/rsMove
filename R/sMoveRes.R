#' @title sMoveRes
#'
#' @description {Tool to support the selection of an adequate satellite spatial resolution. Evaluates how the change
#' in spatial resolution changes the amount of samples and sample regions based on a set of coordinate pairs.}
#' @param x Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param y vector of spatial resolutions (unit depends on spatial projection).
#' @importFrom raster extent
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar
#' @return A \emph{list}.
#' @details {Given a vector of pixel resolutions (\emph{y}), the function determines the number of unique pixels
#' and unique pixel regions after their temporal aggregation. For each spatial resolution, the function starts by converting
#' \emph{x} to unique pixel coordinates and labels them based on their spatial aggregation. Then, the function counts the number
#' of samples and sample regions. The output of the function consists of:
#' \itemize{
#'  \item{\emph{stats} - Summary statistics reporting on the number of unique samples and sample regions per spatial resolution.}
#'  \item{\emph{plot} - Plot representing the change in number of samples and sample regions per spatial resolution.}}}
#' @seealso \code{\link{tMoveRes}} \code{\link{specVar}}
#' @examples {
#'
#'  require(raster)
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # test function for 5, 10 20 and 30 m
#'  a.res <- sMoveRes(shortMove, c(5, 10, 20, 30))
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

sMoveRes <- function(x, y) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!class(x)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"x" is not of a valid class')}
  if (!is.numeric(y)) {stop('"y" is not numeric')}
  if (!is.vector(y)) {stop('"y" is not a vector')}

  # evaluate each resolution
  out <- list() # output variable

#---------------------------------------------------------------------------------------------------------------------#
# 2. find unique sample regions
#---------------------------------------------------------------------------------------------------------------------#

  out <-do.call(rbind, lapply(y, function(r) {

    # reference raster (extend to avoid missing samples along the borders)
    ext <- extend(raster(extent(x), res=r, crs=crs(x)), c(2,2))

    # cell positions of x
    sp <- cellFromXY(ext, x)
    up <- unique(sp)

    # build connected component image
    regions <- clump(rasterize(x, ext, 1))

    # output data frame with statistics
    return(data.frame(nr.pixels=length(up), nr.regions=length(unique(extract(regions, up))), resolution=r))

  }))

#---------------------------------------------------------------------------------------------------------------------#
# 3. plot output
#---------------------------------------------------------------------------------------------------------------------#

  # determine fill scale range
  mv = max(out$nr.regions)
  if (mv < 100) {
    mv <- mv / 10
    fr <- round(mv*2)/2
    if (mv > fr) {fr <- (fr+0.5)*10} else {fr <- fr*10}
  }
  if (mv >= 100) {
    mv <- mv / 100
    fr <- round(mv*20)/20
    if (mv > fr) {fr <- (fr+0.5)*100} else {fr <- fr*100}}

  # determine yscale range
  mv <- max(out$nr.pixels)
  if (mv < 100) {
    mv <- mv / 10
    yr <- round(mv*2)/2
    if (mv > yr) {yr <- (yr+0.5)*10} else {yr <- yr*10}}
  if (mv >= 100) {
    mv <- mv / 100
    yr <- round(mv*20)/20
    if (mv > yr) {yr <- (yr+0.5)*100} else {yr <- yr*100}}

  # make color palette
  cr <- colorRampPalette(c("khaki2", "forestgreen"))

  # build plot object
  out$resolution <- factor(out$resolution, levels=y)
  p <- ggplot(out, aes_string(x="pixel.res", y="nr.pixels", fill="nr.regions")) +
    theme_bw() + scale_fill_gradientn(colors=cr(10), breaks=c(0.0, (fr/2), fr),
    limits=c(0,fr), name="Nr. Regions\n") + xlab("\nResolution (m)") +
    ylab("Nr. Pixels\n") + geom_bar(width=0.7, stat="identity") +
    theme(axis.text.x=element_text(size=12),
          axis.title.x =element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14)) + ylim(0,yr)

  # return data frame and plot
  return(list(stats=out, plot=p))

}
