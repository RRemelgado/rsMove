#' @title sMoveRes
#'
#' @description {Tool to support the selection of an adequate satellite spatial resoltuon. Evaluates how the change
#' in spatial resolution changes the amount of samples and sample regions based on a set of coordinate pairs.}
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param pixel.res vector of spatial resolutions (unit depends on spatial projection).
#' @importFrom raster extent
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar
#' @return A \emph{list}.
#' @details {Given a vector of pixel resolutions (\emph{pixel.res}), the function determines the number of unique pixels
#' and unique pixel regions after their temporal agggregation. For each spatial resolution, the function starts by converting
#' \emph{xy} to unique pixel coordinates and labels them based on their spatial aggregation. Then, the function counts the number
#' of samples and sample regions. The output of the function consists of:
#' \itemize{
#'  \item{\emph{stats} - Summarity statistics reporting on the number of unique samples and sample regions per spatial resolution.}
#'  \item{\emph{plot} - Plot representing the change in number of samples and sample regions per spatial resolution.}
#'  \item{\emph{indices} - Indices for each sample in \emph{xy} based on their spatial aggregation within each spatial resolution.}}}
#' @seealso \code{\link{tMoveRes}} \code{\link{specVar}}
#' @examples {
#'
#'  require(raster)
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData)
#'
#'  # test function for 5, 10 20 and 30 m
#'  a.res <- sMoveRes(xy=moveData, pixel.res=c(5, 10, 20, 30))
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

sMoveRes <- function(xy=xy, pixel.res=pixel.res) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (!is.numeric(pixel.res)) {stop('"pixel.res" is not numeric')}
  if (!is.vector(pixel.res)) {stop('"pixel.res" is not a vector')}

  # evaluate each resolution
  ext <- extent(xy) # reference extent
  out <- list() # output variable
  for (p in 1:length(pixel.res)) {

#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#

    nc <- round((ext[2]-ext[1]) / pixel.res[p]) + 1 # number of columns
    nr <- round((ext[4]-ext[3]) / pixel.res[p]) + 1 # number of rows
    sp <- (round((ext[4]-xy@coords[,2])/pixel.res[p])+1) + nr * round((xy@coords[,1]-ext[1])/pixel.res[p]) # convert coordinates to pixel positions
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

    # update output
    uv <- unique(regions[which(regions>0)])
    out[[p]] <- list(count=length(up), regions=length(uv), indices=sp)

  }

  # output data frame with statistics
  out1 <- data.frame(n.pixels=sapply(out, function(x) {x$count}),
                    n.regions=sapply(out, function(x) {x$regions}))
  row.names(out1) <- as.character(pixel.res)

  # output data frame with sample indices
  out <- lapply(out, function(x) {x$indices})
  out2 <- do.call(cbind, lapply(out, data.frame, stringsAsFactors=FALSE))
  colnames(out2) <- as.character(pixel.res)

#---------------------------------------------------------------------------------------------------------------------#
# 4. plot output
#---------------------------------------------------------------------------------------------------------------------#

  # determine fill scale range
  mv = max(out1$n.regions)
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
  out1$pixel.res <- factor(pixel.res, levels=unique(pixel.res))
  p <- ggplot(out1, aes_string(x="pixel.res", y="n.pixels", fill="n.regions")) +
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
  return(list(stats=out1, plot=p, indices=out2))

}
