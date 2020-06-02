#' @title tMoveRes
#'
#' @description {Tool to support the selection of an adequate satellite temporal resolution. It evaluates how the change in temporal
#' resolution changes the amount of samples and sample regions based on a set of coordinate pairs and their observation dates.}
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.date Object of class \emph{Date} with \emph{xy} observation dates.
#' @param time.res Vector of temporal resolutions (expressed in days).
#' @param pixel.res Spatial resolution (unit depends on spatial projection).
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar element_text
#' @importFrom raster raster extent extend cellFromXY crs
#' @importFrom utils download.file
#' @importFrom grDevices colorRampPalette
#' @return A \emph{list} object reporting on the amount and distribution of unique pixels and connected pixel regions per temporal resolution.
#' @details {Given a base spatial resolution (\emph{pixel.res} and a vector of temporal resolutions (\emph{time.res}), the function determines
#' the number of unique pixels and unique pixel regions after their temporal aggregation. For each temporal resolution, the function starts by
#' converting \emph{xy} to unique pixel coordinates and labels them based on their spatial aggregation. Then, the function counts the number of
#' samples and sample regions. The output of the function consists of:
#' \itemize{
#'  \item{\emph{stats} - Summary statistics reporting on the number of temporal widows, unique samples and unique sample regions per temporal resolution.}
#'  \item{\emph{plot} - Plot representing the change in number of samples and sample regions per temporal resolution.}}}
#' @seealso \code{\link{sMoveRes}} \code{\link{specVar}}
#' @examples {
#'
#'  require(raster)
#'
#'  # reference data
#'  data(longMove)
#'  longMove <- longMove[1:1000,]
#'  # test function for intervals of 1, 8 and 16 days (e.g. of MODIS data)
#'  obs.date <- as.Date(longMove@data$timestamp)
#'  a.res <- tMoveRes(longMove, obs.date, c(1,8,16), 0.1)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

tMoveRes <- function(xy, obs.date, time.res, pixel.res) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (length(pixel.res)>1) {stop('"pixel.res" has more than one element')}
  if (!is.numeric(time.res)) {stop('"time.res" is not numeric')}
  if (class(obs.date)!="Date") {stop('"obs.date" is not of class "Date"')}
  if(sum(is.na(obs.date)) > 0) {stop('please filter missing values in "obs.date"')}
  rp <- crs(xy) # reference projection

  #---------------------------------------------------------------------------------------------------------------------#
  # 2. determine pixel aggregations
  #---------------------------------------------------------------------------------------------------------------------#

  st <- min(obs.date) # start time
  et <- max(obs.date) # end time

  out <- do.call(rbind, lapply(time.res, function(r) {

    nw <- round(as.numeric((et - st)) / r + 1) # number of temporal windows

    tmp <- do.call(rbind, lapply(1:nw, function(w) {

      loc <- which(obs.date >= (st+r*(w-1)) & obs.date <= (((st+r)+(r*w))-1)) # reference samples
      if (length(loc) > 0) {
        ext <- raster(extend(extent(xy[loc, ]), c(pixel.res, pixel.res)), res = pixel.res, crs = rp)
        sp <- cellFromXY(ext, xy[loc, ])
        up <- unique(sp)
        ext[up] <- 1
        regions <- clump(ext)
        return(data.frame(nr.pixels = length(up), nr.regions = cellStats(regions, max, na.rm = TRUE)))
      } else {
        return(data.frame(nr.pixels=0, nr.regions =0))
      }

    }))

    #  estimate final count of pixels/regions
    return(as.data.frame(t(apply(tmp, 2, sum))))

  }))

  row.names(out) <- NULL

#---------------------------------------------------------------------------------------------------------------------#
# 3. plot output
#---------------------------------------------------------------------------------------------------------------------#

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
  out$resolution <- factor(time.res, levels=sort(time.res))
  p <- ggplot(out, aes_string(x="resolution", y="nr.pixels", fill="nr.regions")) + theme_bw() +
    scale_fill_gradientn(colors=cr(10), name="Nr. Regions\n") + xlab("\nResolution (days)") +
    ylab("Nr. Pixels\n") + geom_bar(width=0.7, stat = "identity") +
    theme(axis.text.x=element_text(size=12),
          axis.title.x =element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14)) + ylim(0,yr)

  # return data frame and plot
  return(list(stats=out, plot=p))

}

