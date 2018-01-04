#' @title tMoveRes
#'
#' @description {Tool to support the selection of an adequate satellite temporal resoltuon. It evaluates how the change in temporal
#' resolution changes the amount of samples and sample regions based on a set of coordinate pairs and their observation dates.}
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.time Object of class \emph{Date} with \emph{xy} observation dates.
#' @param time.res Vector of temporal resolutions (expressed in days).
#' @param pixel.res Spatial resolution (unit depends on spatial projection).
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar
#' @importFrom raster raster extent extend cellFromXY crs dim
#' @importFrom utils download.file
#' @importFrom grDevices colorRampPalette
#' @return A \emph{list} object reporting on the amount and distribution of unique pixels and connected pixel regions per temporal resolution.
#' @details {Given a base spatial resolution (\emph{pixel.res} and a vector of temporal resolutions (\emph{time.res}), the function determines
#' the number of unique pixels and unique pixel regions after their temporal agggregation. For each temporal resolution, the function starts by
#' converting \emph{xy} to unique pixel coordinates and labels them based on their spatial aggregation. Then, the function counts the number of
#' samples and sample regions. The output of the function consists of:
#' \itemize{
#'  \item{\emph{stats} - Summarity statistics reporting on the number of temporal widows, unique samples and unique sample regions per temporal resolution.}
#'  \item{\emph{plot} - Plot representing the change in number of samples and sample regions per temporal resolution.}
#'  \item{\emph{indices} - Indices for each sample in \emph{xy} based on their spatial aggregation within each temporal resolution.}}}
#' @seealso \code{\link{sMoveRes}} \code{\link{specVar}}
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
#'  a.res <- tMoveRes(xy=moveData, obs.time=obs.date, time.res=c(1,8,16), pixel.res=0.01)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

tMoveRes <- function(xy=xy, obs.time=obs.time, time.res=time.res, pixel.res=pixel.res) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (length(pixel.res)>1) {stop('"pixel.res" has more than one element')}
  if (!is.numeric(time.res)) {stop('"time.res" is not numeric')}
  if (!class(obs.time)!="Date") {stop('"obs.time" is not of class "Date"')}

#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#

  ext <- extend(raster(extent(xy), res=pixel.res, crs=crs(xy)), c(2,2)) # raster extent
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

  st <- min(obs.time) # start time
  et <- max(obs.time) # end time

  out <- list() # output variable
  for (r in 1:length(time.res)) {

    id <- 0 # reference sample ID
    ind <- vector('numeric', length(xy)) # position index
    nw <- as.numeric(((et - st) / time.res[r]) + 1) # number of temporal windows
    sc <- nr <- vector('numeric', nw) # number of regions

    for (w in 1:nw) {

      # locate pixels within the temporal window
      loc1 <- which(obs.time >= (st+(time.res[r]*(w-1))) & obs.time <= ((st+time.res[r])+(time.res[r]*(w-1))))

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
    out[[r]] <- list(indices=ind, regions=sum(nr), count=sum(sc), window.count=nw)

  }

  # output data frame with statistics
  out1 <- data.frame(n.pixels=sapply(out, function(x) {x$count}),
                     n.regions=sapply(out, function(x) {x$regions}),
                     n.windows=window.count, resolution=time.res)

  # output data frame with sample indices
  out2 <- do.call(cbind, lapply(out, function(x) {x$indices}))
  colnames(out2) <- as.character(time.res)

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
  p <- ggplot(out1, aes(x=factor(time.res), y=n.pixels, fill=n.regions)) + theme_bw() +
    scale_fill_gradientn(colors=cr(10), name="Nr. Regions\n") + xlab("\nResolution (days)") +
    ylab("Nr. Pixels\n") + geom_bar(width=0.7, stat = "identity") +
    theme(axis.text.x=element_text(size=12),
          axis.title.x =element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14)) + ylim(0,yr)

  # return data frame and plot
  return(list(stats=out1, plot=p, indices=out2))

}
