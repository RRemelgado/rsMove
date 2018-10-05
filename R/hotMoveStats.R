#' @title hotMoveStats
#'
#' @description {Segmentation and statistical analysis of the time spent by an animal within a geographical region.}
#' @param x Region unique identifiers. Vector of class \emph{numeric}.
#' @param z Individual identifier. Vector of class \emph{character}.
#' @param y Observation time. Object of class \emph{Date}.
#' @return A list.
#' @importFrom sp spTransform CRS
#' @importFrom ggplot2 ggplot aes_string geom_bar scale_fill_gradientn xlab ylab theme_bw
#' @importFrom grDevices colorRampPalette
#' @details {For each unique region defined by \emph{x}, the function identifies unique temporal segments defined as periods
#' of consecutive days with observations. Then, for each region, the function uses the identified segments to report on the
#' minimum, maximum and mean time spent as well as the total amount of time spent within the region. Moreover, the function
#' provides a detailed report of each segment and informs on the corresponding sample indices. If \emph{z} is specified, the
#' function will in addition count the number of individuals found within each region and within each temporal segment. The
#' final output consists of:
#' \itemize{
#' \item{\emph{region.stats} - \emph{data.frame} with the distribution of samples per region.}
#' \item{\emph{segment.stats} - \emph{data.frame} with all temporal segments assigned to each region.}
#' \item{region.plot} - Plot describing the distribution of samples and recorded time per region.}}
#' @seealso \code{\link{hotMove}}
#' @examples {
#'
#' require(raster)
#'
#' # reference data
#' data(longMove)
#'
#' # extract regions
#' hm <- hotMove(longMove, 0.1, return.shp=TRUE)
#'
#' # plot shapefile (color by region)
#' plot(hm$polygons)
#'
#' # add new information to original shapefile
#' longMove@data <- cbind(longMove@data, hm$region.id)
#'
#' # derive statistics
#' hm.region.stats <- hotMoveStats(hm$region.id, as.Date(longMove@data$timestamp))
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

hotMoveStats <- function(x, y, z) {

  #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  # 1. check input variables
  #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  if (!missing(z)) {
    if (!class(z)%in%c("numeric", "character")) {stop('"z" is not of a valid class')}
    if (length(x)!=length(z)) {stop('"x"/"y" and "z" have different lenghts')}
  } else {z <- replicate(length(x), 1)}
  if (!class(x)%in%c("numeric", "character")) {stop('"x" is not of a valid class')}
  if (class(y)[1]!="Date") {stop('"y" is not of a valid class')}
  if (length(x)!=length(y)) {stop('"x" and "y" have different lenghts')}

  #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
  # 2. evaluate each region
  #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  ur <- sort(unique(x)) # unique regions

  tmp <- lapply(ur, function(r) {

    ri <- which(x == r) # region indices

    # evaluate temporal composition
    y.diff <- c(0, diff(y[ri])) > 1 # running day difference
    seg.ls <- rle(y.diff) # find data sequences
    seg.id <- vector('numeric', length(y.diff)) # label segments (1)
    for (s in 1:length(seg.ls$lengths)) {seg.id[(sum(seg.ls$lengths[0:(s-1)])+1):sum(seg.ls$length[1:s])] <- s} # label segments (1)
    unique.seg <- unique(seg.id) # unique segment ID's

    # single segment stats
    odf1 <- do.call("rbind", lapply(unique.seg, function(s) {
      si <- which(seg.id == s)
      ui <- length(unique(z[ri[si]])) # number of individuals
      ns <- length(which(x == r & y >= y[ri[si[1]]] & y <= y[ri[si[length(si)]]])) # number of samples
      tt <- as.numeric(max(y[ri[si]])-min(y[ri[si]])) + 1
      odf <- data.frame(region.id=r, start.date=min(y[ri[si]]), end.date=max(y[ri[si]]), total.time=tt, nr.individuals=ui, nr.samples=ns, stringsAsFactors=FALSE)
      return(odf)}))

    # region stats
    odf2 <- data.frame(region.id=r, start.date=min(odf1$start.date), end.date=max(odf1$end.date), total.time=sum(odf1$total.time), nr.segments=nrow(odf1),
                       nr.individuals=sum(odf1$nr.individuals), nr.samples=sum(odf1$nr.samples))

    return(list(region=odf2, segment=odf1))

  })

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3. derive output
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  region.df <- do.call("rbind", lapply(tmp, function(i) {i$region})) # region
  segment.df <- do.call("rbind", lapply(tmp, function(i) {i$segment}))# segment report

  rm(tmp)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 4. build plot
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))
  region.df$region.id <- factor(region.df$region.id, levels=sort(unique(region.df$region.id)))
  p <- ggplot(region.df, aes_string(x="region.id", y="total.time", fill="nr.samples")) +
    theme_bw() + geom_bar(stat="identity") + xlab("\nRegion ID") +
    ylab("Number of Days") + scale_fill_gradientn(name="Nr. Samples", colours=cr(10))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# 6. derive output
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

  # return output
  return(list(region.stats=region.df, segment.tats=segment.df, region.plot=p))

}

