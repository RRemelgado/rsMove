#' @title tMoveRes
#'
#' @description Provides historical information on cloud cover.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param o.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param t.res Temporal resolution.
#' @param s.res Spatial resolution.
#' @import ggplot2 sp rgdal grDevices
#' @importFrom utils download.file
#' @return A \emph{list}.
#' @details {Given a vector of temporal resolutions (\emph{t.res}), the function determines
#' the number of unique pixels and unique pixel groups after their temporal agggation. The
#' function returns the corresponding pixel indices per resolution showing which
#' samples would be grouped (\emph{$indices}). The function returns a data frame (\emph{$stats})
#' and a plot (\emph{$plot}) with the statistics per temporal resolution.}
#' @seealso \code{\link{sMoveRes}}
#' @examples {
#'
#'  require(raster)
#'
#'  # reference data
#'  sprj <- crs("+proj=longlat +ellps=WGS84 +no_defs")
#'  moveData <- read.csv(system.file('extdata', 'latlon_example.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,2:3], moveData, proj4string=sprj)
#'
#'  # test function for 5, 10 20 and 30 m
#'  obs.date <- as.Date(moveData@data$timestamp)
#'  a.res <- tMoveRes(xy=moveData, o.time=obs.date, t.res=c(1,8,16), s.res=0.01)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

tMoveRes <- function(xy=xy, o.time=o.time, t.res=t.res, s.res=s.res) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (length(s.res)>1) {stop('"s.res" has more than one element')}
  if (!is.numeric(t.res)) {stop('"t.res" is not numeric')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#

  ext <- extend(raster(extent(xy), res=s.res, crs=crs(xy)), c(2,2)) # raster extent
  rd <- dim(ext)# raster dimensions
  nr <- rd[1] # number of rows
  nc <- rd[2] # number of columns
  sp <- cellFromXY(ext, xy) # unique pixels
  up <- unique(sp) # unique pixel positions

#---------------------------------------------------------------------------------------------------------------------#
# 3. find unique sample regions
#---------------------------------------------------------------------------------------------------------------------#

  # evaluate pixel connectivity
  regions <- matrix(0, nr, nc)
  for (r in 1:length(up)) {
    rp <- ((up[r]-1) %% nr)+1
    cp <- ((up[r]-1) %/% nr)+1
    if (cp > 1) {sc<-cp-1} else {sc<-cp}
    if (cp < nc) {ec<-cp+1} else {ec<-cp}
    if (rp > 1) {sr<-rp-1} else {sr<-rp}
    if (rp < nr) {er<-rp+1} else {er<-rp}
    if (max(regions[sr:er,sc:ec])>0) {
      uv <- unique(regions[sr:er,sc:ec])
      uv <- uv[which(uv > 0)]
      mv <- min(uv)
      regions[rp,cp] <- mv
      for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
    } else {regions[rp,cp] <- max(regions)+1}
  }

  # update region ID's
  uregions <- sapply(up, function(x) {regions[x]})

  rm(regions)

#---------------------------------------------------------------------------------------------------------------------#
# 4. determine pixel aggregations
#---------------------------------------------------------------------------------------------------------------------#

  st <- min(o.time) # start time
  et <- max(o.time) # end time

  out <- list() # output variable
  for (r in 1:length(t.res)) {

    id <- 0 # reference sample ID
    ind <- vector('numeric', length(xy)) # position index
    nw <- as.numeric(((et - st) / t.res[r]) + 1) # number of temporal windows
    sc <- nr <- vector('numeric', nw) # number of regions

    for (w in 1:nw) {

      # locate pixels within the temporal window
      loc1 <- which(o.time >= (st+(t.res[r]*(w-1))) &
                     o.time <= ((st+t.res[r])+(t.res[r]*(w-1))))

      # quantify unique samples
      upr <- unique(sp[loc1])
      for (p in 1:length(upr)) {
        id <- id + 1
        loc2 <- which(sp[loc1]==upr[p])
        ind[loc1[loc2]] <- id}

      # derive statistics
      nr[w] <- length(unique(uregions[up%in%upr])) # number of regions
      sc[w] <- length(upr) # number of samples

    }

    # update output
    out[[r]] <- list(indices=ind, regions=sum(nr), count=sum(sc), regions.window=nr, count.window=sc)

  }

  # output data frame with statistics
  out1 <- data.frame(n.pixels=sapply(out, function(x) {x$count}),
                     n.regions=sapply(out, function(x) {x$regions}),
                     resolution=t.res)

  # output data frame with sample indices
  out2 <- do.call(cbind, lapply(out, function(x) {x$indices}))
  colnames(out2) <- as.character(t.res)

  # count per window
  out3 <- lapply(out, function(x) {list(regions.window=x$regions.window, count.window=x$count.window)})

  #---------------------------------------------------------------------------------------------------------------------#
  # 4. plot output
  #---------------------------------------------------------------------------------------------------------------------#

  # determine yscale range
  mv <- max(out1$n.pixels)
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
  p <- ggplot(out1, aes(x=factor(t.res), y=n.pixels, fill=n.regions)) + theme_bw() +
    scale_fill_gradientn(colors=cr(10), name="Nr. Regions\n") + xlab("\nResolution (days)") +
    ylab("Nr. Pixels\n") + geom_bar(width=0.7, stat = "identity") +
    theme(axis.text.x=element_text(size=12),
          axis.title.x =element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14)) + ylim(0,yr)

  # return data frame and plot
  return(list(stats=out1, plot=p, indices=out2, window.stats=out3))

}
