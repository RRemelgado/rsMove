#' @title rsComposite
#'
#' @description {Phenological and date driven Pixel Based Compositing (PBC)
#' of remote sensing data supported by GPS tracking date information.}
#' @param img Object of class \emph{RasterSpack} or \emph{RasterBrick}.
#' @param r.dates Object of class \emph{Date} with \emph{img} observation dates.
#' @param obs.time Object of class \emph{Date} with reference dates.
#' @param n.dev Number of deviations from the target date. Default is 1.
#' @param type One of "norm" or "pheno" or "pheno2".
#' @param d.buffer Search buffer (expressed in days).
#' @import raster rgdal
#' @importFrom stats lm
#' @seealso \code{\link{imgInt}} \code{\link{dataQuery}} \code{\link{proSat}}
#' @return A \emph{list}.
#' @details {The function uses a multi-layer raster object to build a composite.
#' It looks at a ginve set of dates (e.g. GPS tracking dates) and estimates a
#' reference date to build the composite for defined by the median of \emph{obs.time}.
#' The median is then used to estimate Median Absolute Deviation (MAD) which
#' specifies the size of the buffer set aroung the target date within which
#' bands will be considered. Here, \emph{n.dev} is used as a multiplier to enlarge
#' the temporal buffer. Alternatively, a user define temporal buffer is allowed
#' by using the keyword \emph{d.buffer}. If \emph{obs.time} countains only one element,
#' the function will use it as a reference date. In this case, if \emph{d.buffer}
#' is NULL the function will set it to 30 by default. The way how the function handles
#' temporal information depends on the \emph{type} keyword. If set to \emph{norm},
#' the function will search for the nearest possible dates within the temporal
#' buffer. However, if \emph{pheno} is set, then the day of the year will be given
#' priority. Thus, if multi-year raster data is provided, older data with a DOY
#' closer to the target that will be used when possible. The output provides:
#' #' \itemize{
#'  \item{\emph{value} - composite of target images}
#'  \item{\emph{dates} - per pixel date code}
#'  \item{\emph{count} - pixel count of \emph{dates}}
#'  \item{\emph{na.count} - count of NA values}
#'  \item{\emph{target} - target date}
#'  \item{\emph{mad} - temporal buffer}}
#'  If \emph{pheno2} is used, for each pixel, the function wilhin estimate a weighted
#'  mean of the clear pixels within the temporal buffer. The weights represent the
#'  inverse time difference between the target and the available dates giver higher
#'  weights to small differences.}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'  rsStk <- stack(rsStk, rsStk, rsStk) # dummy files for the example
#'
#'  # raster dates
#'  r.dates = seq.Date(as.Date("2013-01-01"), as.Date("2013-12-31"), 45)
#'
#'  # target date
#'  obs.time = as.Date("2013-06-01")
#'
#'  # build composite
#'  r.comp <- rsComposite(rsStk, r.dates, obs.time, d.buffer=90)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

rsComposite <- function(img, r.dates, obs.time, n.dev=1, type='norm', d.buffer=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  # raster
  if (!exists('img')) {stop('error: "img" is missing')}
  if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('error: "img" is not of a valid class')}

  # raster dates
  if (!exists('r.dates')) {stop('error: "r.dates" is missing')}
  if (!class(r.dates)[1]%in%c('Date')) {stop('error: "r.dates" is nof of a valid class')}
  if (length(r.dates)!=nlayers(img)) {stop('errorr: "img" and "r.dates" have different lengths')}

  # reference dates
  if (!exists('obs.time')) {stop('"obs.time" is missing')}
  if (!class(obs.time)[1]%in%c('Date')) {stop('error: "obs.time" is nof of a valid class')}

  # auxiliary variables
  if (!is.numeric(n.dev)) {stop('"n.dev" is not numeric')}
  if (!is.null(d.buffer)) {if (!is.numeric(d.buffer)) {stop('"d.buffer" is not numeric')}}
  if (!type%in%c('norm', 'pheno')) {stop('"type" is not a valid keyword')}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. handle time information
#-------------------------------------------------------------------------------------------------------------------------------#

  # determine date format
  if (type=='norm') {
    ot0 <- obs.time
    rd0 <- r.dates}
  if (type=='pheno' | type=='pheno2') {
    bd <- as.Date(paste0(as.character(format(obs.time,'%Y')), '-01-01'))
    ot0 <- as.numeric((obs.time-bd) + 1)
    bd <- as.Date(paste0(as.character(format(r.dates,'%Y')), '-01-01'))
    rd0 <- as.numeric((r.dates-bd) + 1)}

  # determine date range
  if (length(obs.time)>1) {
    t.date <- median(ot0) # target date
    if (is.null(d.buffer)) {d.buffer <- median(abs(ot0-t.date))*n.dev} # search buffer
  } else {
    t.date <- ot0
    if (is.null(d.buffer)) {d.buffer <- 30}
  }


  if (type=="norm" | type=="pheno") {

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. define functions
#-------------------------------------------------------------------------------------------------------------------------------#

    # compositing
    f1 <- function(x) {
      ind <- which(!is.na(x))
      if (length(ind)>0) {
        v <- x[ind]
        d <- rd0[ind]
        diff <- abs(d-t.date)
        ind <- which(diff <= d.buffer)
        if (length(ind)>0) {v[ind[which(diff[ind]==min(diff[ind]))]]
        } else {return(NA)}} else {return(NA)}}

    # retrieve date
    f2 <- function(x) {
      ind <- which(!is.na(x))
      if (length(ind)>0) {
        v <- x[ind]
        d <- rd0[ind]
        diff <- abs(d-t.date)
        ind <- which(diff <= d.buffer)
        if (length(ind)>0) {as.numeric(d[ind[which(diff[ind]==min(diff[ind]))]])
        } else {return(NA)}} else {return(NA)}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. build composites
#-------------------------------------------------------------------------------------------------------------------------------#

    r.value <- calc(img, f1)
    r.date <- calc(img, f2)

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. build composites
#-------------------------------------------------------------------------------------------------------------------------------#

    # build table from date info
    ud <- unique(r.date)
    used <- rd0[as.numeric(rd0)%in%ud]
    ud <- as.numeric(rd0)
    count <- sapply(ud, function(x) {cellStats(r.date==x, sum)})
    df <- data.frame(date=used, value=ud, count=count)

    # check for missing values
    mv <- cellStats(is.na(r.date), sum)

    # return data
    return(list(value=r.value, date=r.date, count=df, na.count=mv, target=t.date, buffer=d.buffer))

  } else {

    # averaging function
    f1 <- function(x) {
      ind <- which(!is.na(x))
      if (length(ind)>0) {
        v <- x[ind]
        d <- rd0[ind]
        diff <- abs(d-t.date)
        w <- 1 - (diff/max(diff))
        ind <- which(diff <= d.buffer)
        if (length(ind)>0) {sum(v[ind]*w[ind]) / sum(w[ind])} else {
          return(NA)}} else {return(NA)}}

    # apply function to raster
    return(calc(img, f1))

  }

}
