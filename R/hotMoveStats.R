#' @title hotMoveStats
#'
#' @description Provides statistics for the output of hotMove.
#' @param rid List object as provided by \emph{hotMove()}.
#' @param aid Optional. Unique identifiers.
#' @param o.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct}.
#' @param tUnit Time unit for stats. Default is \emph{days}. See \code{\link[base]{difftime}} for additional keywords.
#' @param method Method used to estimate polygon area. 
#' @return A \emph{data frame}.
#' @import sp rgdal ggplot2
#' @details {This functions analysis the attributes of sample regions define by hotMove(). Alternatively, 
#' the user can keep \emph{rid} as NULL. This case, all information will be assumed as part of one region. 
#' For each sample region, the function returns the amount of samples (\emph{tns}. If a vector of unique 
#' identifiers is provided (\emph{aid}) the number of unique identifiers observed within each region is also 
#' reported. If temporal information is provided (\emph{time}) the function identifies unique temporal segment 
#' corresponding to periods of consecutive days with observations. For each segment, the function reports on the 
#' amount of segments (\emph{nts}) as well as the minimum (\emph{mnt}), maximum (\emph{mxt}) and mean (\emph{avt}) 
#' of the time segments and the total amount of time that they amount to (\emph{tts}). For each region, the function 
#' will also report on the start and end of each temporal segment (\emph{$temporal.segments}) and will provide the 
#' sample indices for associated to each segment (\emph{$segment.indices}). If \emph{rid} contains polygons for each 
#' region, the function also reports on the area of convex polygons. In this case, the user can use the \emph{method} 
#' keyword to specify how the polygons should be handled. If the polygons are in lat-lon, method \emph{deg} can be used 
#' to re-project each polygon to its corresponding UTZ zone before retrieving the area. This is the default, if the polygons 
#' are in a cartesian coordinate system use \emph{m}.}
#' @seealso \code{\link{hotMove}}
#' @examples {
#' 
#' require(raster)
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
#' hm.region.stats <- hotMoveStats(rid=hm, o.time=as.Date(moveData@data$timestamp))
#' hm.time.stats <- hotMoveStats(o.time=as.Date(moveData@data$timestamp))
#' 
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

hotMoveStats <- function(rid=NULL, o.time=NULL, tUnit=NULL, aid=NULL, method='deg') {

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  if(is.null(rid) & is.null(o.time) & is.null(aid)) {stop('no information provided ("rid", "o.time" and "aid" are NULL)')}
  if (!is.null(o.time)) {
    if (!class(o.time)[1]%in%c("Date", "POSIXlt", "POSIXt")) {stop('"o.time" is not a valid "POSIXlt"/POSIXt" object')}
    pt <- TRUE} else {pt <- FALSE}
  if(!is.null(o.time) & !is.null(aid)) {if (length(aid)!=length(o.time)) {stop('"o.time" and "aid" have different lengths')}}
  if (!(method %in% c('m', 'deg')) & !is.null(rid$polygons)) {stop(paste0('shapefile provided in "rid" but method ', method, ' not valid (choose between "m" and "deg")'))}
  if (!is.null(rid)) {
    if (!is.list(rid)) {stop('"rid" is not a list')}
    if (min(names(rid)%in%c('indices', 'polygons'))==0) {stop('names in "rid" are not valid (requires "indices" and/or "polygons"')}
    if (!is.null(o.time)) {if (length(rid$indices)!=length(o.time)) {stop('"o.time" and "rid" have different lengths')}}
    if (!is.null(aid)) {if (length(rid$indices)!=length(aid)) {stop('"aid" and "rid" have different lengths')}}
    ur <- sort(unique(rid$indices)) # unique regions
    nr <- length(ur) # number of unique regions
    pr <- TRUE
  } else {
    ur <- 1
    nr <- 1
    pr <- FALSE}
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 2. define output stats
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  ura <- matrix(0, nr) # area of unique regions
  tns <- matrix(0, nr) # number of samples per region
  nui <- matrix(0, nr) # number of individuals per region
  mnt <- matrix(0, nr) # smallest o.time segment per region
  avt <- matrix(0, nr) # average of o.time segments per region
  mxt <- matrix(0, nr) # largest o.time segment per region
  tts <- matrix(0, nr) # sum of recorded o.time
  nts <- matrix(0, nr) # number of o.time segments
  
  if (pt) {
    ss1 <- list() # temporal segment stats
    ss2 <- list() # temporal segment indices
  }

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3. evaluate each region
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    
  # evaluate each region separately
  for (r in 1:length(ur)) {

    # extract base stats
    if (pr) {ind <- which(rid$indices==ur[r])} else {ind <- 1:length(o.time)}
    if (!is.null(aid)) {nui[r] = length(unique(aid[ind]))} else {nui[r]<-NA}
    tns[r] = length(ind)
    
    # identify unique temporal segments and count number of days
    if (!is.null(o.time)) {
      ts0 <- list()
      sp <- 1
      st <- sort(unique(o.time[ind]))
      if (length(st) > 1) {
        for (t in 2:length(st)) {
          diff <- as.numeric(difftime(st[t], st[(t-1)], units=tUnit))
          if (diff > 1) {
            ts0[[(length(ts0)+1)]] <- as.numeric(difftime(st[(t-1)], st[sp], units=tUnit)) + 1
            ss1[[length(ss1)+1]] <- data.frame(start=st[sp], end=st[(t-1)], id=ur[r], count=length(sp:(t-1)))
            ss2[[length(ss2)+1]] <- which(o.time >= ss1[[length(ss1)]]$start & 
                                            o.time <= ss1[[length(ss1)]]$end & 
                                            rid$indices==ur[r])
            sp <- t
          }}
        ts0 <- unlist(ts0)
        if(!is.null(ts0)) {
          mnt[r] <- min(ts0)
          avt[r] <- mean(ts0)
          mxt[r] <- max(ts0)
          tts[r] <- sum(ts0)
          nts[r] <- length(ts0)
        } else {
          tts[r] <- as.numeric(difftime(max(st), min(st), units=tUnit)) + 1
          mnt[r] <- tts[r]
          avt[r] <- tts[r]
          mxt[r] <- tts[r]
          nts[r] <- 1
          ss1[[length(ss1)+1]] <- data.frame(start=min(st), end=max(st), id=ur[r], count=length(st))
          ss2[[length(ss2)+1]] <- which(o.time >= ss1[[length(ss1)]]$start & 
                                          o.time <= ss1[[length(ss1)]]$end & 
                                          rid$indices==ur[r])
        }
        rm(ts0, st, diff, sp)
      } else {
        mnt[r] <- 1
        mxt[r] <- 1
        tts[r] <- 1
        nts[r] <- 1
        ss1[[length(ss1)+1]] <- data.frame(start=min(st), end=max(st), id=ur[r], count=length(st))
        ss2[[length(ss2)+1]] <- which(o.time >= ss1[[length(ss1)]]$start & 
                                        o.time <= ss1[[length(ss1)]]$end & 
                                        rid$indices==ur[r])}
    } else {
      mnt <- NA
      mat <- NA
      mxt <- NA
      tts <- NA
      nts <- NA
    }

    if (pr) {
      
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
  
      } else {ura[r] <- NA}} else {ura[r] <- NA}

  }
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 4. build output data frames
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  # build data frame with statistics
  df <- data.frame(rid=ur, tns=tns, nui=nui, mnt=mnt, avt=avt, mxt=mxt, tts=tts, nts=nts, ura=ura)
  
  # temporal segments
  if (pt) {time.seg <- do.call(rbind, ss1)}
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 5. build output data frames
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  # color scheme
  if (pr) {
    cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))
    p1 <- ggplot(df, aes(x=factor(rid, levels=unique(rid)), y=tts, fill=tns)) + theme_bw() + 
      geom_bar(stat="identity") + xlab("\nRegion ID") + ylab("Number of Days") + 
      scale_fill_gradientn(name="Nr. Samples", colours=cr(10))}

  if (pt & !pr) {
    p2 <- ggplot(time.seg, aes(x=start, y=count)) + theme_bw() + geom_bar(stat="identity") + 
      xlab("Start time") + ylab("Sample Count")}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 6. derive output
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  # return output
  if (pt & pr) {return(list(stats=df, plot=p1, temporal.segments=time.seg, segment.indices=ss2))}
  if (pt & !pr) {return(list(stats=df, plot=p2, temporal.segments=time.seg, segment.indices=ss2))}
  if (!pt) {return(list(stats=df))}
  
}