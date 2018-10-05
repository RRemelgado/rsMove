#' @title moveSeg
#'
#' @description Pixel based segmentation of movement data using environmental data.
#' @param x Object of class \emph{RasterLayer} or \emph{data.frame}.
#' @param z Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct}.
#' @param y Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param data.type Raster data data.type. One of \emph{cont} (continuous) or \emph{cat} (for categorical).
#' @param threshold Change threshold. Required if \emph{data.type} is set to \emph{cat}.
#' @param summary.fun Summary function used to summarize the values within each segment when \emph{method} is \emph{cont}. Default is mean.
#' @param buffer.size Spatial buffer size applied around each segment (unit depends on spatial projection).
#' @param smooth.fun Smoothing function applied with \emph{buffer.size} when \emph{method} is \emph{cont}. Default is mean.
#' @importFrom raster extract crs
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar scale_fill_discrete element_blank element_text
#' @seealso \code{\link{dataQuery}} \code{\link{imgInt}} \code{\link{timeDir}} \code{\link{spaceDir}}
#' @return A \emph{list}.
#' @details {This function identifies segments of comparable environmental conditions along a movement track given by \emph{y}.
#' Looking at consecutive data points, the function queries \emph{x} and proceeds to identify a new segment if \emph{threshold}
#' is exceeded. Then, for each segment, the function summarizes \emph{x} using \emph{summary.fun} and reports on the amount of
#' points found within it. Moreover, if \emph{z} is set, the function reports on the start and end timestamps and the elapsed
#' time for each segment. If \emph{data.type} is set as \emph{'cont'}, the function assumes \emph{y} is a continuous variable. This
#' will require the user to define \emph{threshold} which indicates when the difference between consecutive points should be
#' considered as a change. If \emph{data.type} is set as \emph{'cat'}, then the function will ignore \emph{threshold} and map
#' a change every time a change in value occurs. The user might choose to smooth the extracted values using \emph{buffer.size}.
#' This will prompt the function to summarize the values arounde ach sample in \emph{y} using a metric define by\emph{smooth.fun}.
#' However, if \emph{data.type} is set to \emph{cat}, \emph{smooth.fun} is ignored. In this case, the function will report on the
#' majority value within the spatial buffer. The output of this function consists of:
#'\itemize{
#'  \item{\emph{segment.id} - Vector reporting on the segment identifiers associated to each sample in \emph{y}.}
#'  \item{\emph{segment.stats} - Statistical information for each segment reporting on the corresponding environmental and temporal information.}
#'  \item{\emph{segment.plot} - plot of \emph{stats} showing the variability of environmental conditions and time spent per segment.}}}
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
#'  z <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # perform directional sampling
#'  seg <- moveSeg(r, z, shortMove, data.type="cat")
#'
#' }
#' @export

#---------------------------------------------------------------------------------------------------------------------#

moveSeg <- function(x, z, y, data.type='cont', threshold=NULL, summary.fun=NULL, buffer.size=NULL, smooth.fun=NULL) {

  #---------------------------------------------------------------------------------------------------------------------#
  # 1. check input variables
  #---------------------------------------------------------------------------------------------------------------------#

  # sample dates
  if (!class(z)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"z" is nof of a valid class')}
  if (length(z)!=length(y)) {stop('errorr: "y" and "z" have different lengths')}

  # environmental data
  if (class(x)[1]=='RasterLayer') {
    if (crs(y)@projargs!=crs(x)@projargs) {stop('"y" and "x" have different projections')}
    if (missing(y)) {stop('"x" is a "RasterLayer", "y" is required')}
    if (!class(y)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"y" is not of a valid class')}
    rProj <- crs(y) # output projection
    prd <- TRUE
  } else {
    if (!class(x)[1]%in%c('RasterLayer', 'data.frame')) {stop('"x" is neither a raster or a data frame')}
    if (nrow(x)!=length(y)) {stop('number of elements in "y" and "x" do not match')}
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
  if (data.type=='cat') {summary.fun <- function(x) {return(x[1])}}

  # check buffer size
  if (!is.null(buffer.size)) {
    if (!is.numeric(buffer.size)) {stop('"buffer.size" is not numeric')}
    if (length(buffer.size) > 1) {stop('"buffer.size" has more than 1 element')}
  }

  #---------------------------------------------------------------------------------------------------------------------#
  # 2. query data
  #---------------------------------------------------------------------------------------------------------------------#

  if (prd) {

    # apply buffer if required
    if (!is.null(buffer.size)) {

      # average samples within buffer
      if (data.type=='cont') {x <- extract(x, y@coords, buffer=buffer.size, fun=smooth.fun, na.rm=T)}

      # determine main class within the buffer
      if (data.type=='cat') {

        # dilate samples
        tmp <- lapply(1:length(y), function(p) {
          ind <- raster(extent((y@coords[p,1]-buffer.size), (y@coords[p,1]+buffer.size),
                               (y@coords[p,2]-buffer.size), (y@coords[p,2]+buffer.size)), crs=rProj)
          ind <- xyFromCell(ind, 1:ncell(ind))
          return(list(c=ind, s=replicate(nrow(ind), p)))})
        si <- unlist(lapply(tmp, function(p) {p$s}))
        tmp <- do.call(rbind, lapply(tmp, function(p) {p$c}))

        # extract values
        edata0 <- extract(x, tmp)

        # sumarize data (extract dominant class)
        x <- sapply(1:length(y), function(p) {
          ind <- which(si==p)
          r0 <- as.vector(edata0[ind[!duplicated(cellFromXY(p, tmp[ind,1:2]))]])
          uc <- unique(r0)
          uc <- uc[!is.na(uc)]
          if (length(uc)>0) {
            count <- sapply(uc, function(p) {sum(r0==p)})
            return(uc[which(count==max(count))[1]])
          } else {return(NA)}})

        rm(tmp, si, edata0)

      }

      # simple query
    } else {x <- extract(x, y@coords)}}

#---------------------------------------------------------------------------------------------------------------------#
# 3. identify segments
#---------------------------------------------------------------------------------------------------------------------#

  # search for segments and return sample indices
  pd <- rle(x)$lengths
  seg.id <- vector('numeric', length(sp))
  for (p in 1:length(pd)) {seg.id[(sum(pd[0:(p-1)])+1):sum(pd[1:p])] <- p}


#---------------------------------------------------------------------------------------------------------------------#
# 4. derive segment statistics
#---------------------------------------------------------------------------------------------------------------------#

  id <- unique(seg.id) # unique ID's

  odf <- do.call(rbind, lapply(id, function(u) {

    ii <- which(seg.id == u) # target samples
    rv <- summary.fun(x[ii]) # segment raster value

    st <- min(z[ii]) # start time
    et <- max(z[ii]) # end time
    tt <- as.numeric(difftime(et, st, units="mins")) # elapsed time

    return(data.frame(sid=u, count=length(ii), start=st, end=et, elapsed=tt, value=rv, stringsAsFactors=FALSE))}

  ))

#---------------------------------------------------------------------------------------------------------------------#
# 5. build plot
#---------------------------------------------------------------------------------------------------------------------#

  if (data.type=='cont') {

    odf$sid <- factor(odf$sid, levels=sort(unique(odf$sid)))

    # plot with time
    if (!is.null(z)) {

      # buid plot object
      p <- ggplot(odf, aes_string(x="sid", y="time", fill="value")) + geom_bar(stat="identity") +
        xlab('') + ylab('Time Spent (minutes)') +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))

      # plot without time
    } else {

      # determine y range scale range
      mv <- max(odf$count)
      nc <- nchar(as.character(mv))
      m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
      mv <- mv / m
      yr <- round(mv)
      if (mv > yr) {yr <- (yr+0.2)*m} else {yr <- yr*m}

      # buid plot object
      p <- ggplot(odf, aes_string(x="sid", y="count", fill="value")) + geom_bar(stat="identity") +
        xlab('\nSegment ID') + ylab('Sample count\n') +
        theme(axis.text.x=element_text(size=10, hjust=1), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}

  if (data.type=='cat') {

    odf$sid <- factor(odf$sid, levels=sort(unique(odf$sid)))
    odf$value <- factor(odf$value, levels=sort(unique(odf$value)))

    # plot with time
    if (!is.null(z)) {

      # buid plot object
      p <- ggplot(odf, aes_string(x="sid", y="elapsed", fill="value")) + geom_bar(stat="identity") +
        xlab('\nSegment ID') + ylab('Time Spent (minutes)\n') + scale_fill_discrete(name="Class") +
        theme(axis.text.x=element_text(size=10, hjust=1), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))

      # plot without time
    } else {

      # buid plot object
      p <- ggplot(odf, aes_string(x="sid", y="count", fill="value")) + geom_bar(stat="identity") +
        xlab('\nSegment ID') + ylab('Sample count\n') + scale_fill_discrete(name="Class") +
        theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size=10, hjust=1), axis.text=element_text(size=12),
              axis.title=element_text(size=14), legend.title=element_text(size=14), legend.text=element_text(size=12))}}

  #---------------------------------------------------------------------------------------------------------------------#
  # 6. return output
  #---------------------------------------------------------------------------------------------------------------------#

  colnames(odf) <- c("segment.id", "nr.samples", "start.time", "end.time", "total.time", "summary.value")
  return(list(segment.id=seg.id, segment.stats=odf, segment.plot=p))

}

