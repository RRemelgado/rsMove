#' @title hotMoveStats
#'
#' @description Provides statistics for the output of hotMove.
#' @param rid Output of hotMove or an object with the same structure (i.e. list containing indices under the keyword $indices and region polygons of the class "SpatialPolygon" under the keyword "polygons")
#' @param aid Optional. Animal unique id. Use to evaluate the number of individuals per region.
#' @param time Optional. Valid "POSIXlt" or "POSIXct" type vector with the time of observation for each record.
#' @param method Method used to estimate polygon area. Required if rid contains a shapefile. If "deg" the function assumes that the used projectn is in degrees. If so, the program projects each polygon to the closest UTM zone so that an area can be estimated. Default is 'deg'.
#' @return Statistics for a sample regions derived with hotMove(). The function provides information on: Total number of samples (tns), numbdf of unique individuals (nui), time of the smallest segment (mnt), average time per segment (avt), time of largest segment (mxt), number of time segments (nts) and area of the polygon (ura). Missing information is reported as NA
#' @import sp
#' @seealso \code{\link{hotMove}},  \code{\link{CRS}}
#' @examples \dontrun{
#' # reference data
#' whiteStork <- data(WhiteStork_2013.rdata)
#' 
#' # lat/lon WGS84
#' outProj <-  CRS("+init=epsg:4978")
#' 
#' # determine regions
#' hm <- hotMove(x=whiteStork$lon, y=whiteStork=lon, res=0.01, cs=outProj) 
#' 
#' # extract statistics
#' hotMoveStats(rid=hm, time=whiteStork$time, aid=whiteStork$id, method='m')
#' }

#-------------------------------------------------------------------------------------------------------------------------------#

hotMoveStats <- function(rid=rid, time=NULL, aid=NULL, method='deg') {

  # check if input is valid
  if (!exists(rid)) {stop('error: rid not provided')}
  if (!is.null(time)) {if (max(c("POSIXlt", "POSIXt") %in% class(a))!=1) {stop('error: "time" is not a valid "POSIXlt"/POSIXt" object')}}
  if (!is.null(aid)) {if (length(rid$indices)!=length(aid)) {stop('error: "aid" not the same length as "rid"')}}
  if(!is.null(time)) {if (length(rid$indices)!=length(aid)) {stop('error: "time" not the same length as "rid"')}}
  if (is.null(rid$indices)) {stop('error: provided list is not valid ("indices" keyword missing)')}
  if (!(method %in% c('m', 'deg')) & !is.null(rid$polygons)) {stop(paste0('error: shapefile provided in "rid" but method ', method, ' not valid (choose between "m" and "deg")'))}

  # output varibles
  ur <- unique(rid$indices) # unique regions
  nr <- length(ur) # number of unique regions
  ura <- matrix(0, nr) # area of unique regions
  tns <- matrix(0, nr) # number of samples per region
  nui <- matrix(0, nr) # number of individuals per region
  mnt <- matrix(0, nr) # small time segment per region
  avt <- matrix(0, nr) # average of time segments per region
  mxt <- matrix(0, nr) # largest time segment per region
  nts <- matrix(0, nr) # number of tme intervals
  # evaluate each region separately
  for (r in 1:ur) {

    # extract base stats
    ind <- which(rid$indices==ur[r])
    tns[r] = length(ind)
    if (!is.null(aid)) {nui[r] = length(unique(aid[ind]))} else {nui[r]<-NA}

    # identify unique temporal segments and count number of days
    if (!is.null(time)) {
      ts <- list()
      sp <- 1
      st <- sort(time[ind])
      for (t in 2:length(st)) {ts[]}
      diff <- as.numeric(difftime(st[(t-1)], st[t], units='days'))
      if (diff > 1) {
        ts[[(length(ts)+1)]] <- as.numeric(difftime(st[sp], st[(t-1)], units='days'))
        sp <- t
      }
      ts <- unlist(ts)
      mnt[r] <- min(ts)
      mat[r] <- mean(ts)
      mxt[r] <- max(ts)
      nts[r] <- length(ts)
      
      rm(ts, st, diff, sp)

    } else {
      mnt <- NA
      mat <- NA
      mxt <- NA
      nts <- NA
    }

    # if provided, estimate polygon area
    if (!is.null(rid$polygons)) {

      # if the data is in degrees, project it to UTM
      if (method=='deg') {

        # find target UTM zone
        sc <- rid$polygons[r]@polygons[[1]]@Polygons[[1]]@coords[,2] # polygon vertices
        zc <- vector('character', length(sc))
        zone <- sapply(sc, function(x) {if (x > 0) {return(round(x/6.)+31)} else {return(round((180+x)/6)+1)}})
        uv = unique(zone)
        count <- vector('numeric', length(uv)) # zone code counts
        for (u in 1:length(uv)) {count[u] <- length(zone==uv[u])}
        zone <- zone[which(count==max(count))] # dominant zone/orientation
        refProj <- sp::CRS(paste0('+proj=utm +zone=', zone, ' +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'))

        # project polygon and estimate area
        shp <- sp::spTransform(rid$polygons[r], refProj)
        ura[r] <- shp@polygons[[1]]@Polygons[[1]]@area

        rm(sc, zc, zone, count, refProj)

      } else {ura[r] <- rid$polygons[r]@polygons[[1]]@Polygons[[1]]@area}

    } else {ura[r] <- NA}

  }

  # build/return data frame with statistics
  return(data.frame(tns, nui, mnt, avt, mxt, nts, ura, stringsAsFactors=F))

}
