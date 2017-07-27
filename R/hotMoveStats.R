#' @title hotMoveStats
#'
#' @description Provides statistics for the output of hotMove.
#' @param rid List object as provided by \emph{hotMove()}.
#' @param aid Optional. Unique identifiers.
#' @param time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct}.
#' @param tUnit Time unit for stats. Default is \emph{days}. See \code{\link[base]{difftime}} for additional keywords.
#' @param method Method used to estimate polygon area. 
#' @return A \emph{data frame}.
#' @import sp rgdal
#' @details {This functions analysis the attributes of sample regions define by hotMove().
#' For each sample region, the function returns the amount of samples (\emph{tns}. If a vector 
#' of unique identifiers is provided (\emph{aid}) the number of unique identifiers observed within 
#' each region is also reported. If temporal information is provided (\emph{time}) the function 
#' identifies unique temporal segment corresponding to periods of consecutive days with 
#' observations. For each segment, the function reports on the amount of segment s(\emph{nts}) 
#' as well as the minimum (\emph{mnt}), maximum (\emph{mxt}) and mean (\emph{avt}), ) of the time segments 
#' and the total amount of time that they amount to (\emph{tts}). If \emph{rid} contains polygons for 
#' each region, the function also reports on the area of convex polygons. In this case, the 
#' user can use the \emph{method} keyword to specify how the polygons should be handled. If the 
#' polygons are in lat-lon, method \emph{deg} can be used to re-project each polygon to its 
#' corresponding UTZ zone before retrieving the area. This is the default, if the polygons 
#' are in a cartesian coordinate system use \emph{m}.}
#' @seealso \code{\link{hotMove}}
#' @examples {
#' 
#' require(sp)
#' 
#' # reference data
#' sprj <- CRS("+proj=longlat +ellps=WGS84 +no_defs")
#' moveData <- read.csv(system.file('extdata', 'latlon_example.csv', package="rsMove"))
#' moveData <- SpatialPointsDataFrame(moveData[,2:3], moveData, proj4string=sprj)
#' 
#' # extract regions
#' hm <- hotMove(xy=moveData, pxr=0.1, shp=TRUE)
#' 
#' # plot shapefile (color by region)
#' plot(hm$polygons, col=hm$indices)
#' 
#' # add new information to original shapefile
#' moveData@data <- cbind(moveData@data, hm$indices)
#' 
#' # derive statistics
#' hm.stats <- hotMoveStats(rid=hm, time=as.Date(moveData@data$timestamp))
#' 
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

hotMoveStats <- function(rid=rid, time=NULL, tUnit=NULL, aid=NULL, method='deg') {

  # check if input is valid
  if (!exists('rid')) {stop('rid not provided')}
  if (!is.null(time)) {if (!class(time)[1]%in%c("Date", "POSIXlt", "POSIXt")) {stop('"time" is not a valid "POSIXlt"/POSIXt" object')}}
  if (!is.null(aid)) {if (length(rid$indices)!=length(aid)) {stop('"aid" not the same length as "rid"')}}
  if(!is.null(time)) {if (length(rid$indices)!=length(time)) {stop('"time" not the same length as "rid"')}}
  if (is.null(rid$indices)) {stop('provided list is not valid ("indices" keyword missing)')}
  if (!(method %in% c('m', 'deg')) & !is.null(rid$polygons)) {stop(paste0('shapefile provided in "rid" but method ', method, ' not valid (choose between "m" and "deg")'))}
  
  # output varibles
  ur <- sort(unique(rid$indices)) # unique regions
  nr <- length(ur) # number of unique regions
  ura <- matrix(0, nr) # area of unique regions
  tns <- matrix(0, nr) # number of samples per region
  nui <- matrix(0, nr) # number of individuals per region
  mnt <- matrix(0, nr) # smallest time segment per region
  avt <- matrix(0, nr) # average of time segments per region
  mxt <- matrix(0, nr) # largest time segment per region
  tts <- matrix(0, nr) # sum of recorded time
  nts <- matrix(0, nr) # number of time segments
  
  # evaluate each region separately
  for (r in 1:length(ur)) {

    # extract base stats
    ind <- which(rid$indices==ur[r])
    tns[r] = length(ind)
    if (!is.null(aid)) {nui[r] = length(unique(aid[ind]))} else {nui[r]<-NA}
  
    # identify unique temporal segments and count number of days
    if (!is.null(time)) {
      ts <- list()
      sp <- 1
      st <- sort(unique(time[ind]))
      if (length(st) > 1) {
        for (t in 2:length(st)) {
          diff <- as.numeric(difftime(st[t], st[(t-1)], units='days'))
          if (diff > 1) {
            ts[[(length(ts)+1)]] <- as.numeric(difftime(st[(t-1)], st[sp], units='days')) + 1
            sp <- t
          }}
        ts <- unlist(ts)
        if(!is.null(ts)) {
          mnt[r] <- min(ts)
          avt[r] <- mean(ts)
          mxt[r] <- max(ts)
          tts[r] <- sum(ts)
          nts[r] <- length(ts)
        } else {
          tts[r] <- as.numeric(difftime(max(st), min(st), units='days')) + 1
          mnt[r] <- tts[r]
          avt[r] <- tts[r]
          mxt[r] <- tts[r]
          nts[r] <- 1}
        
        rm(ts, st, diff, sp)
      } else {
        mnt[r] <- 1
        mxt[r] <- 1
        tts[r] <- 1
        nts[r] <- 1}
    } else {
      mnt <- NA
      mat <- NA
      mxt <- NA
      tts <- NA
      nts <- NA
    }

    # if provided, estimate polygon area
    if (!is.null(rid$polygons)) {

      # if the data is in degrees, project it to UTM
      if (method=='deg') {

        # find target UTM zone
        sc <- rid$polygons[r]@polygons[[1]]@Polygons[[1]]@coords[,2] # polygon vertices
        zone <- sapply(sc, function(x) {if (x > 0) {return(round(x/6.)+31)} else {return(round((180+x)/6)+1)}})
        uv = unique(zone)
        count <- sapply(uv, function(x){sum(zone==x)}) # zone code counts
        zone <- zone[which(count==max(count))] # dominant zone/orientation
        refProj <- sp::CRS(paste0('+proj=utm +zone=', zone[1], ' +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

        # project polygon and estimate area
        shp <- sp::spTransform(rid$polygons[r], refProj)
        ura[r] <- shp@polygons[[1]]@Polygons[[1]]@area * 0.000001

        rm(sc, zone, count, refProj)

      } else {ura[r] <- rid$polygons[r]@polygons[[1]]@Polygons[[1]]@area * 0.000001}

    } else {ura[r] <- NA}

  }

  # build/return data frame with statistics
  return(data.frame(rid=ur, tns=tns, nui=nui, mnt=mnt, avt=avt, mxt=mxt, tts=tts, nts=nts, ura=ura))

}