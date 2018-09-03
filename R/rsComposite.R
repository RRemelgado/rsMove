#' @title rsComposite
#'
#' @description {Phenological and date driven Pixel Based Compositing (PBC).}
#' @param x Object of class \emph{RasterStack} or \emph{RasterBrick}.
#' @param x.dates Object of class \emph{Date} with \emph{x} observation dates.
#' @param obs.dates Object of class \emph{Date} with reference dates.
#' @param comp.method One of "closest" or "phenological". The default is "closest".
#' @param temporal.buffer Search buffer (expressed in days). The default is NULL.
#' @importFrom raster nlayers calc cellStats
#' @importFrom stats lm
#' @seealso \code{\link{imgInt}} \code{\link{dataQuery}}
#' @return A \emph{list}.
#' @details {The function uses a multi-layer raster object to build a composite for a reference date which
#' corresponds to the median of \emph{obs.dates}. Moreover,the function estimates the Median Absolute Deviation
#' (MAD) of \emph{obs.dates} which determines the temporal buffer that is used to search for usable data. As an
#' alternative, \emph{temporal.buffer} can be specified manually and will be required if \emph{obs.dates} consists
#' of a single value. The user can also specify how the compositing should be done. \emph{comp.method} can be set to:
#' #' \itemize{
#'  \item{\emph{closest} - Selects pixels from layers with the closest possible date.}
#'  \item{\emph{phenological} - Selects pixels from layers with the closest Day of the Year (DoY).}}
#' The final output of \emph{rsComposite} is a list consisting of:
#' #' \itemize{
#'  \item{\emph{composite} - Final image composite}
#'  \item{\emph{dates} - Temporal composition of the composite reporting on the julian day (origin is "1970-01-01").}
#'  \item{\emph{pixel.count} - pixel count of unique values in \emph{dates}. Additionally, it reports on NA values.}
#'  \item{\emph{pixel.count.plot} - Plot with relative frequency of pixels per date extracted from \emph{pixel.count}.}
#'  \item{\emph{target.date} - Reference date used during compositing.}
#'  \item{\emph{temporal.buffer} - Temporal buffer used during compositing.}}
#'  If \emph{pheno2} is used, for each pixel, the function will estimate a weighted
#'  mean of the clear pixels within the temporal buffer. The weights represent the
#'  inverse time difference between the target and the available dates giver higher
#'  weights to small differences.}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # read raster data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)
#'  r.stk <- stack(file)
#'  r.stk <- stack(r.stk, r.stk, r.stk) # dummy files for the example
#'
#'  # raster dates
#'  file.name <- names(r.stk)
#'  x.dates <- as.Date(paste0(substr(file.name, 2, 5), '-',
#'  substr(file.name, 7, 8), '-', substr(file.name, 10, 11)))
#'
#'  # target date
#'  obs.dates = as.Date("2013-06-01")
#'
#'  # build composite
#'  r.comp <- rsComposite(r.stk, x.dates, obs.dates, comp.method="closest", temporal.buffer=90)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

rsComposite <- function(x, x.dates, obs.dates, comp.method='closest', temporal.buffer=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  # raster
  if (!class(x)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"x" is not of a valid class')}

  # raster dates
  if (!class(x.dates)[1]%in%c('Date')) {stop('"x.dates" is nof of a valid class')}
  if (length(x.dates)!=nlayers(x)) {stop('"x" and "x.dates" have different lengths')}

  # reference dates
  if (!class(obs.dates)[1]%in%c('Date')) {stop('"obs.dates" is nof of a valid class')}

  # auxiliary variables
  if (!is.null(temporal.buffer)) {if (!is.numeric(temporal.buffer)) {stop('"temporal.buffer" is not numeric')}}
  if (length(obs.dates) == 1 & is.null(temporal.buffer)) {stop('"obs.dates" has one element (provide specify "temporal.buffer"')}
  if (!comp.method%in%c('closest', 'phenological')) {stop('"comp.method" is not a valid keyword')}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. handle time information
#-------------------------------------------------------------------------------------------------------------------------------#

  # determine date range
  if (length(obs.dates)>1) {
    t.date <- median(obs.dates) # target date
    if (is.null(temporal.buffer)) {temporal.buffer <- median(abs(obs.dates-t.date))} # search buffer
  } else {t.date <- obs.dates}

  # determine date format
  if (comp.method=='closest') {
    obs.dates <- obs.dates
    x.dates0 <- x.dates}
  if (comp.method=='phenological' | comp.method=='pheno2') {
    bd <- as.Date(paste0(as.character(format(obs.dates,'%Y')), '-01-01'))
    obs.dates <- as.numeric((obs.dates-bd) + 1)
    bd <- as.Date(paste0(as.character(format(x.dates,'%Y')), '-01-01'))
    x.dates0 <- as.numeric((x.dates-bd) + 1)}

  t.date1 <- as.numeric(julian.Date(t.date, origin=as.Date("1970-01-01")))
  x.dates1 <- as.numeric(julian.Date(x.dates, origin=as.Date("1970-01-01")))

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. define functions
#-------------------------------------------------------------------------------------------------------------------------------#

  # compositing function
  f1 <- function(i) {
    ind <- which(!is.na(i) & abs(x.dates0-t.date) < temporal.buffer)
    if (length(ind) > 0) {return(i[ind[which(i[ind] == min(i[ind]))[1]]])} else {return(NA)}}

  # function to report on temporal composition
  f2 <- function(i) {
    ind <- which(!is.na(i) & abs(x.dates0-t.date) < temporal.buffer)
    if (length(ind) > 0) {return(x.dates1[which(i[ind] == min(i[ind]))[1]])} else {return(NA)}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. build composites
#-------------------------------------------------------------------------------------------------------------------------------#

  r.value <- calc(x, f1)
  r.date <- calc(x, f2)

#-------------------------------------------------------------------------------------------------------------------------------#
# 5. build output
#-------------------------------------------------------------------------------------------------------------------------------#

  # build table and plot from date info
  ud <- sort(unique(x.dates1))
  count <- sapply(ud, function(i) {cellStats(r.date==i, sum)})
  ud <- do.call("c", lapply(ud, function(i) {x.dates[which(x.dates1==i)[1]]}))
  odf <- data.frame(date=ud, absolute.count=count, relative.count=count/sum(count))
  p <- ggplot(odf, aes_string(x="date", "relative.count")) + geom_bar(stat="identity") +
    xlab("Date") + ylab("Pixel Frequency (%)")

  # return data
  return(list(composite=r.value, dates=r.date, pixel.count=odf, pixel.count.plot=p,
              target.date=t.date, temporal.buffer=temporal.buffer))

}
