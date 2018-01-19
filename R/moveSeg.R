#' @title moveSeg
#'
#' @description Pixel based segmentation of movement data using environmental data.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param env.data Object of class \emph{RasterLayer} or \emph{data.frame}.
#' @param data.type Raster data data.type. One of \emph{cont} (continous) or \emph{cat} (for categorical).
#' @param obs.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{xy} observation dates.
#' @param threshold Change threshold. Required if \emph{data.type} is set to \emph{cat}.
#' @param summary.fun Summary function used to summarize the values within each segment when \emph{method} is \emph{cont}. Default is mean.
#' @param buffer.size Spatial buffer size applied around each segment (unit depends on spatial projection).
#' @param smooth.fun Smoothing function applyed with \emph{buffer.size} when \emph{method} is \emph{cont}. Default is mean.
#' @importFrom raster extract crs
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar scale_fill_discrete
#' @import raster rgdal ggplot2
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}} \code{\link{timeDir}} \code{\link{spaceDir}}
#' @return A \emph{list}.
#' @details {This function identifies segments of comparable environmental conditions along the movement track given by \emph{xy}.
#' Looking at consective data points, the function queries \emph{env.data} and proceeds to identify a new segment if \emph{threshold}
#' is exceeded. Then, for each segment, the function summarizes \emph{env.data} using \emph{summary.fun} and reports on the amount
#' of points found within it. Moreover, if \emph{obs.time} is set, the function reports on the start and end timestamps and the elapsed time.
#' If \emph{method} is set as \emph{'cont'}, the function assumes the raster data is a continuous variable. This will require the user to
#' define \emph{theshold} which indicates when the difference between consecutive points should be considered a change.
#' In order to smooth the extracted values the user can specify \emph{buffer.size}. This will prompt the function to summarize the values around
#' each sample in \emph{xy} using a metric define by\emph{smooth.fun}. However, if \emph{data.type} is set to \emph{cat} \emph{smooth.fun} is ignored.
#' In this case, the function will report on the majority value within the buffer. The output of this function consists of:
#'\itemize{
#'  \item{\emph{indices} - Vector reporting on the segment identifiers associated to each sample in \emph{xy}.}
#'  \item{\emph{stats} - Statistical information for each segment reporting on the conrresponding environmental and temporal information.}
#'  \item{\emph{plot} - plot of \emph{stats} showing the variability of environmental conditions and time spent per segment.}}}
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', 'landCover.tif', package="rsMove"))
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # observation time
#'  obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # perform directional sampling
#'  seg <- moveSeg(xy=shortMove, obs.time=obs.time, env.data=r, data.type="cat")
#'
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------#

moveSeg <- function(xy=xy, env.data=env.data, data.type='cont', threshold=threshold,
                    obs.time=NULL, summary.fun=NULL, buffer.size=NULL, smooth.fun=NULL) {

  #---------------------------------------------------------------------------------------------------------------------#
  # 1. check input variables
  #---------------------------------------------------------------------------------------------------------------------#

  # samples
  if (!exists('xy')) {stop('"xy" is missing')}
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rProj <- crs(xy) # output projection

  # sample dates
  if (!is.null(obs.time)) {
    if (!class(obs.time)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"obs.time" is nof of a valid class')}
    if (length(obs.time)!=length(xy)) {stop('errorr: "xy" and "obs.time" have different lengths')}
    processTime <- TRUE
  } else {processTime <- FALSE}

  # environmental data
  if (class(env.data)[1]=='RasterLayer') {
    if (crs(xy)@projargs!=crs(env.data)@projargs) {stop('"xy" and "env.data" have different projections')}
    prd <- TRUE
  } else {
    if (!class(env.data)[1]%in%c('data.frame')) {stop('"env.data" is neither a raster or a data frame')}
    if (nrow(env.data)!=length(xy)) {stop('number of elements in "xy" and "env.data" do not match')}
    prd=FALSE}

  # check threshold
  if (data.type=='cont') {
    if (is.null(threshold)) {stop('"data.type" is set to "cont". Please define "threshold"')}
    if (!is.numeric(threshold)) {stop('"threshold" is not numeric')}}
  if (data.type=='cat') {threshold <- 1}

  # check query data.type
  if (!data.type%in%c('cont', 'cat')) {stop('"data.type" is not  avalid keyword')}
  if (data.type=='cont') {
    if (!is.null(smooth.fun)) {if (!is.function(smooth.fun)){stop('"smooth.fun" is not a valid keyword or function')}}
    if (is.null(smooth.fun)) {smooth.fun <- function(x) {return(mean(x, na.rm=T))}}
    if (!is.null(summary.fun)) {if (!is.function(smooth.fun)){stop('"summary.fun" is not a valid keyword or function')}}
    if (is.null(summary.fun)) {summary.fun <- function(x) {return(mean(x, na.rm=T))}}}
  if (data.type=='cat') {summary.fun <- function(x) {return(mean(x, na.rm=T))}}

  #---------------------------------------------------------------------------------------------------------------------#
  # 2. query data
  #---------------------------------------------------------------------------------------------------------------------#

  if (prd) {

    # apply buffer if required
    if (!is.null(buffer.size)) {

      # average samples within buffer
      if (data.type=='cont') {env.data <- extract(env.data, xy@coords, buffer=buffer.size, fun=smooth.fun, na.rm=T)}

      # determine main class within the buffer
      if (data.type=='cat') {

        # dilate samples
        tmp <- lapply(1:length(xy), function(x) {
          ind <- raster(extent((xy@coords[x,1]-buffer.size), (xy@coords[x,1]+buffer.size),
                               (xy@coords[x,2]-buffer.size), (xy@coords[x,2]+buffer.size)), crs=rProj)
          ind <- xyFromCell(ind, 1:ncell(ind))
          return(list(c=ind, s=replicate(nrow(ind), x)))})
        si <- unlist(lapply(tmp, function(x) {x$s}))
        tmp <- do.call(rbind, lapply(tmp, function(x) {x$c}))

        # extract values
        edata0 <- extract(env.data, tmp)

        # sumarize data (extract dominant class)
        env.data <- sapply(1:length(xy), function(x) {
          ind <- which(si==x)
          r0 <- as.vector(edata0[ind[!duplicated(cellFromXY(env.data, tmp[ind,1:2]))]])
          uc <- unique(r0)
          uc <- uc[!is.na(uc)]
          if (length(uc)>0) {
            count <- sapply(uc, function(x) {sum(r0==x)})
            return(uc[which(count==max(count))[1]])
          } else {return(NA)}})

        rm(tmp, si, edata0)

      }

      # simple query
    } else {env.data <- extract(env.data, xy@coords)}}

  #---------------------------------------------------------------------------------------------------------------------#
  # 3. identify segments
  #---------------------------------------------------------------------------------------------------------------------#

  # search for segments
  r0 <- 1
  li <- 1
  id <- list() # segment id
  rv <- list() # segment value

  for (r in 2:length(xy)) {
    diff <- abs(env.data[r]-env.data[(r-1)])
    if (!is.na(diff)) {
      if (diff >= threshold) {
        ep <- r-1
        rv[[li]] <- summary.fun(env.data[c(r0:ep)])
        id[[li]] <- replicate(length(c(r0:ep)), li)
        r0 <- r
        li <- li + 1
        if (r==length(xy)) {
          id[[li]] <- li
          rv[[li]] <- env.data[r]}
      } else {if (r==length(xy)) {
        ep <- r
        rv[[li]] <- summary.fun(env.data[c(r0:ep)])
        id[[li]] <- replicate(length(c(r0:ep)), li)}}}}
  rv <- unlist(rv)
  id <- unlist(id)

  #---------------------------------------------------------------------------------------------------------------------#
  # 4. derive segment statistics
  #---------------------------------------------------------------------------------------------------------------------#

  # build region report
  uid <- sort(unique(id))
  df <- do.call(rbind, lapply(uid, function(x) {
    ind <- which(id==x)
    if (processTime) {
      tt <- difftime(obs.time[ind[length(ind)]], obs.time[ind[1]], units="mins")
      st <- min(obs.time[ind])
      et <- max(obs.time[ind])
    } else {
      tt <- NA
      st <- NA
      et <- NA}
    return(data.frame(sid=x, count=length(ind), start=st, end=et,
                      elapsed=as.numeric(tt), value=rv[x], stringsAsFactors=FALSE))}))

  #---------------------------------------------------------------------------------------------------------------------#
  # 5. build plot
  #---------------------------------------------------------------------------------------------------------------------#

  if (data.type=='cont') {

    df$sid <- factor(df$sid, levels=sort(unique(df$sid)))

    # plot with time
    if (!is.null(obs.time)) {

      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="time", fill="value")) + geom_bar(stat="identity") +
        xlab('') + ylab('Time Spent (minutes)') +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))

      # plot without time
    } else {

      # determine y range scale range
      mv <- max(df$count)
      nc <- nchar(as.character(mv))
      m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
      mv <- mv / m
      yr <- round(mv)
      if (mv > yr) {yr <- (yr+0.2)*m} else {yr <- yr*m}

      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="count", fill="value")) + geom_bar(stat="identity") +
        xlab('\nSegment ID') + ylab('Sample count\n') +
        theme(axis.text.x=element_text(size=10, hjust=1), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}

  if (data.type=='cat') {

    df$sid <- factor(df$sid, levels=sort(unique(df$sid)))
    df$value <- factor(df$value, levels=sort(unique(df$value)))

    # plot with time
    if (!is.null(obs.time)) {

      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="elapsed", fill="value")) + geom_bar(stat="identity") +
        xlab('\nSegment ID') + ylab('Time Spent (minutes)\n') + scale_fill_discrete(name="Class") +
        theme(axis.text.x=element_text(size=10, hjust=1), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))

      # plot without time
    } else {

      # buid plot object
      p <- ggplot(df, aes_string(x="sid", y="count", fill="value")) + geom_bar(stat="identity") +
        xlab('\nSegment ID') + ylab('Sample count\n') + scale_fill_discrete(name="Class") +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size=10, hjust=1), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}

  #---------------------------------------------------------------------------------------------------------------------#
  # 6. return output
  #---------------------------------------------------------------------------------------------------------------------#

  colnames(df) <- c("Segment ID", "Nr. Samples", "Timestamp (start)", "Timestamp (end)", "Total time (minutes)", "Value")
  return(list(indices=id, stats=df, plot=p))

}
