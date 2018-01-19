#' @title hotMoveStats
#'
#' @description {Segmentation and statistical analysis of the time spent by an animal within a geographical region.}
#' @param region.id Region unique identifiers. Vector of class \emph{numeric}.
#' @param individual.id Individual identifier. Vector of class \emph{character}.
#' @param obs.time Observation time. Object of class \emph{Date}.
#' @return A list containing statistical information for each region (\emph{region.stats}) and for each temporal segment (\emph{temporal.segment.stats}) and sample indices for each segment (temporal.segment.indices)
#' @importFrom sp spTransform CRS
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_fill_gradientn xlab ylab theme_bw
#' @importFrom grDevices colorRampPalette
#' @details {For each unique region defined by \emph{region.id}, the function identifies unique temporal segments
#' defined as periods of consecutive days with observations. Then, for each region, the function uses the identified segments
#' to repoprt on the minimum, maximum and mean time spent as well as the total amount of time spent within the region.
#' Moreover, the function provides a detailed report of each segment and informs on the corresponding sample indices. If
#' \emph{individual.id} is specified, the function will in addition count the number of individuals found within each region
#' and within each temporal segment.}
#' @seealso \code{\link{hotMove}}
#' @examples {
#'
#' require(raster)
#'
#' # reference data
#' data(longMove)
#'
#' # extract regions
#' hm <- hotMove(xy=longMove, pixel.res=0.1, return.shp=TRUE)
#'
#' # plot shapefile (color by region)
#' plot(hm$polygons)
#'
#' # add new information to original shapefile
#' longMove@data <- cbind(longMove@data, hm$indices)
#'
#' # derive statistics
#' hm.region.stats <- hotMoveStats(region.id=hm$indices,
#' obs.time=as.Date(longMove@data$timestamp))
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

hotMoveStats <- function(region.id=region.id, obs.time=obs.time, individual.id=NULL) {

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  if(is.null(region.id) & is.null(obs.time) & is.null(individual.id)) {stop('no information provided ("region.id", "obs.time" and "individual.id" are NULL)')}
  if (!class(region.id)%in%c("numeric", "character")) {stop('"region.id" is not of a valid class')}
  if (class(obs.time)[1]!="Date") {stop('"obs.time" is not of a valid class')}
  if (!is.null(individual.id)) {if (!class(individual.id)%in%c("numeric", "character")) {stop('"individual.id" is not of a valid class')}}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 2. define output stats
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  ur <- sort(unique(region.id)) # unique regions
  nr <- length(ur) # number of unique regions

  ura <- matrix(0, nr) # area of unique regions
  tns <- matrix(0, nr) # number of samples per region
  nui <- matrix(0, nr) # number of individuals per region
  mnt <- matrix(0, nr) # smallest obs.time segment per region
  avt <- matrix(0, nr) # average of obs.time segments per region
  mxt <- matrix(0, nr) # largest obs.time segment per region
  tts <- matrix(0, nr) # sum of recorded obs.time
  nts <- matrix(0, nr) # number of obs.time segments
  ss1 <- list() # temporal segment stats
  ss2 <- list() # temporal segment indices

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3. evaluate each region
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  # evaluate each region separately
  for (r in 1:length(ur)) {

    # extract base stats
    ind <- which(region.id==ur[r])
    if (!is.null(individual.id)) {nui[r] = length(unique(individual.id[ind]))} else {nui[r]<-NA}
    tns[r] = length(ind)

    # identify unique temporal segments and count number of days
    ts0 <- list()
    sp <- 1
    si <- order(obs.time[ind])
    st <- obs.time[ind[si]]
    ui <- individual.id[ind[si]]
    if (length(st) > 1) {
      for (t in 2:length(st)) {
        diff <- as.numeric(difftime(st[t], st[(t-1)], units="days")) + 1
        if (diff > 1) {
          ts0[[(length(ts0)+1)]] <- as.numeric(difftime(st[(t-1)], st[sp], units="days")) + 1
          ss1[[length(ss1)+1]] <- data.frame(start=st[sp], end=st[(t-1)], id=ur[r], count=length(sp:(t-1)), individuals=length(unique(ui[sp:(t-1)])))
          ss2[[length(ss2)+1]] <- which(obs.time >= ss1[[length(ss1)]]$start & obs.time <= ss1[[length(ss1)]]$end & region.id==ur[r])
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
        tts[r] <- as.numeric(difftime(max(st), min(st), units="days")) + 1
        mnt[r] <- tts[r]
        avt[r] <- tts[r]
        mxt[r] <- tts[r]
        nts[r] <- 1
        ss1[[length(ss1)+1]] <- data.frame(start=min(st), end=max(st), id=ur[r], count=length(st), individuals=length(unique(ui)))
        ss2[[length(ss2)+1]] <- which(obs.time >= ss1[[length(ss1)]]$start & obs.time <= ss1[[length(ss1)]]$end & region.id==ur[r])
      }
      rm(ts0, st, diff, sp)
    } else {
      mnt[r] <- 1
      mxt[r] <- 1
      tts[r] <- 1
      nts[r] <- 1
      ss1[[length(ss1)+1]] <- data.frame(start=min(st), end=max(st), id=ur[r], count=length(st), individuals=length(unique(ui)))
      ss2[[length(ss2)+1]] <- which(obs.time >= ss1[[length(ss1)]]$start & obs.time <= ss1[[length(ss1)]]$end & region.id==ur[r])}

  }

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 4. build output data frames
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  # build data frame with statistics
  df <- data.frame(region.id=ur, tns=tns, nui=nui, mnt=mnt, avt=avt, mxt=mxt, tts=tts, nts=nts)

  # temporal segments
  time.seg <- do.call(rbind, ss1)
  colnames(time.seg) <- c("Start Time", "End Time", "Region ID", "Nr Samples", "Nr Individuals")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 5. build output data frames
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))
  df$region.id <- factor(df$region.id, levels=unique(df$region.id))
  p <- ggplot(df, aes(x=region.id, y=nts, fill=tns)) +
    theme_bw() + geom_bar(stat="identity") + xlab("\nRegion ID") +
    ylab("Number of Days") + scale_fill_gradientn(name="Nr. Samples", colours=cr(10))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 6. derive output
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  # return output
  colnames(df) <- c("Region ID", "Nr Samples", "Nr Individuals", "Min Time", "Average Time", "Max Time", "Total Time", "Nr Segments")
  return(list(region.stats=df, plot=p, temporal.segment.stats=time.seg, temporal.segment.indices=ss2))

}
