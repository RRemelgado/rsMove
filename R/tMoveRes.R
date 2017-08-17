#' @title tMoveRes
#'
#' @description Provides historical information on cloud cover.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param o.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param t.res Temporal resolution.
#' @param s.res Spatial resolution.
#' @param p.res Should the output be ploted on screen? Default is TRUE.
#' @import ggplot2 sp rgdal grDevices
#' @importFrom utils download.file
#' @return A \emph{list}.
#' @details {Given a vector of temporal resolutions (\emph{t.res}), the function determines 
#' the number of unique pixels and unique pixel groups after their temporal agggation. The 
#' function returns the corresponding pixel indices per resolution showing which 
#' samples would be grouped (\emph{$indices}). The function returns a data frame (\emph{$stats}) 
#' and a plot (\emph{$plot}) with the statistics per temporal resolution.}
#' @seealso \code{\link{sMoveRes}}
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
#'  # test function for 5, 10 20 and 30 m
#'  a.res <- tMoveRes(xy=moveData, dpath='.')
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

tMoveRes <- function(xy=xy, o.time=o.time, t.res=t.res, s.res=s.res, p.res=T) {
  
#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#
  
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rr <- crs(xy) # reference projection
  if (is.na(rr@projargs)) {stop('"xy" does not have a valid projection')}
  if (length(s.res)>1) {stop('"s.res" has more than one element')}
  if (!is.logical(p.res)) {stop('"p.res" is not a logical argument')}
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. determine grid coordinates for given pixels
#---------------------------------------------------------------------------------------------------------------------#
  
  ext <- extent(xy) # reference extent
  nc <- round((ext[2]-ext[1]) / s.res) + 1 # number of columns
  nr <- round((ext[4]-ext[3]) / s.res) + 1 # number of rows
  sp <- (round((ext[4]-xy@coords[,2])/s.res)+1) + nr * round((xy@coords[,1]-ext[1])/s.res) # convert coordinates to pixel positions
  
#---------------------------------------------------------------------------------------------------------------------#
# 3. find unique sample regions
#---------------------------------------------------------------------------------------------------------------------#

  up <- unique(sp) # unique pixel positions
  
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
    nr <- vector('numeric', length(nw)) # number of regions
    
    for (w in 1:length(nw)) {
      
      # locate pixels within the temporal window
      loc1 <- which(o.time >= (st+(t.res[r]*(w-1))) & 
                     o.time <= ((st+t.res[r])+(t.res[r]*(w-1))))
       
      # quantify unique samples
      upr <- unique(sp[loc1])
      for (p in 1:length(upr)) {
        id <- id + 1
        loc2 <- which(sp[loc1]==upr[p])
        ind[loc1[loc2]] <- id}
      
      # number of unique regions
      nr[w] <- length(unique(uregions[up%in%upr]))}
    
    # update output
    out[[r]] <- list(count=max(ind), indices=ind, regions=sum(nr))
    
  }
  
  # output data frame with statistics
  out1 <- data.frame(n.pixels=sapply(out, function(x) {x$count}), 
                     n.regions=sapply(out, function(x) {x$regions}))
  row.names(out1) <- as.character(t.res)
  
  # output data frame with sample indices
  out <- lapply(out, function(x) {x$indices})
  out2 <- do.call(cbind, lapply(out, data.frame, stringsAsFactors=FALSE))
  colnames(out2) <- as.character(t.res)
  
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
  p <- ggplot(out1, aes(x=factor(t.res), y=n.pixels, fill=n.regions)) + theme_bw() + 
    scale_fill_gradientn(colors=cr(10), breaks=c(0.0, (fr/2), fr), 
                         limits=c(0,fr), name="Nr. Regions\n") + xlab("\nResolution (days)") + 
    ylab("Nr. Pixels\n") + geom_bar(width=0.7, stat = "identity") + 
    theme(axis.text.x=element_text(size=12), 
          axis.title.x =element_text(size=14), 
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12), 
          legend.title=element_text(size=14)) + ylim(0,yr)
  
  if (p.res) {p} # plot on screen
  
  # return data frame and plot
  return(list(stats=out1, plot=p, indices=out2))
  
}