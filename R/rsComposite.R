#' @title rsComposite
#'
#' @description {Phenological and date driven Pixel Based Compositing (PBC)
#' of remote sensing data supported by GPS tracking date information.}
#' @param img Object of class \emph{RasterSpack} or \emph{RasterBrick}.
#' @param img.dates Object of class \emph{Date} with \emph{img} observation dates.
#' @param obs.dates Object of class \emph{Date} with reference dates.
#' @param comp.method One of "closest" or "phenological".
#' @param temporal.buffer Search buffer (expressed in days). The default is 30.
#' @importFrom raster nlayers calc cellStats
#' @importFrom stats lm
#' @seealso \code{\link{imgInt}} \code{\link{dataQuery}}
#' @return A \emph{list}.
#' @details {The function uses a multi-layer raster object to build a composite for
#' a reference date which corresponds to the median of \emph{obs.dates}. Moreover,
#' the function determines the Median Absolute Deviation (MAD) of \emph{obs.dates}
#' which determines the temporal buffer that is used to search for usable images.
#' As an alternative, \emph{temporal.buffer} can be specified manually and will
#' be required if \emph{obs.dates} consists of a single value. The user can also
#' specify how the compositing should be done. \emph{comp.method} can be set to:
#' #' \itemize{
#'  \item{\emph{closest} - Uses layer with the closest possible date in relaton to the reference date.}
#'  \item{\emph{phenological} - Uses the layer with the Day of the Year (DoY) in relation to the reference date.}}
#' The final output of \emph{rsComposite} is a list consisting of:
#' #' \itemize{
#'  \item{\emph{composite} - Final image composite}
#'  \item{\emph{dates} - Temporal composition of the composite reporting on the julian day}
#'  \item{\emph{pixel.count} - pixel count of unique values in \emph{dates}. Additionally, it reports on NA values.}
#'  \item{\emph{target.date} - Reference date used during compositing.}
#'  \item{\emph{temporal.buffer} - Temporal buffer used during compositing.}}
#'  If \emph{pheno2} is used, for each pixel, the function wilhin estimate a weighted
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
#'  img.dates <- as.Date(paste0(substr(file.name, 2, 5), '-',
#'  substr(file.name, 7, 8), '-', substr(file.name, 10, 11)))
#'
#'  # target date
#'  obs.dates = as.Date("2013-06-01")
#'
#'  # build composite
#'  r.comp <- rsComposite(r.stk, img.dates, obs.dates, comp.method="closest", temporal.buffer=90)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

rsComposite <- function(img, img.dates, obs.dates, comp.method='closest', temporal.buffer=NULL) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  # raster
  if (!class(img)[1]%in%c('RasterStack', 'RasterBrick')) {stop('"img" is not of a valid class')}

  # raster dates
  if (!class(img.dates)[1]%in%c('Date')) {stop('"img.dates" is nof of a valid class')}
  if (length(img.dates)!=nlayers(img)) {stop('"img" and "img.dates" have different lengths')}

  # reference dates
  if (!class(obs.dates)[1]%in%c('Date')) {stop('"obs.dates" is nof of a valid class')}

  # auxiliary variables
  if (!is.null(temporal.buffer)) {if (!is.numeric(temporal.buffer)) {stop('"temporal.buffer" is not numeric')}}
  if (!comp.method%in%c('closest', 'phenological')) {stop('"comp.method" is not a valid keyword')}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. handle time information
#-------------------------------------------------------------------------------------------------------------------------------#

  # determine date format
  if (comp.method=='closest') {
    obs.dates <- obs.dates
    img.dates <- img.dates}
  if (comp.method=='phenological' | comp.method=='pheno2') {
    bd <- as.Date(paste0(as.character(format(obs.dates,'%Y')), '-01-01'))
    obs.dates <- as.numeric((obs.dates-bd) + 1)
    bd <- as.Date(paste0(as.character(format(img.dates,'%Y')), '-01-01'))
    img.dates <- as.numeric((img.dates-bd) + 1)}

  # determine date range
  if (length(obs.dates)>1) {
    t.date <- median(obs.dates) # target date
    if (is.null(temporal.buffer)) {temporal.buffer <- median(abs(obs.dates-t.date))} # search buffer
  } else {
    t.date <- obs.dates
    if (is.null(temporal.buffer)) {temporal.buffer <- 30}
  }


#-------------------------------------------------------------------------------------------------------------------------------#
# 3. define functions
#-------------------------------------------------------------------------------------------------------------------------------#

  # compositing function
  f1 <- function(x) {
    ind <- which(!is.na(x))
    if (length(ind)>0) {
      v <- x[ind]
      d <- img.dates[ind]
      diff <- abs(d-t.date)
      ind <- which(diff <= temporal.buffer)
      if (length(ind)>0) {v[ind[which(diff[ind]==min(diff[ind]))]]
      } else {return(NA)}} else {return(NA)}}

  # function to report on temporal composition
  f2 <- function(x) {
    ind <- which(!is.na(x))
    if (length(ind)>0) {
      v <- x[ind]
      d <- img.dates[ind]
      diff <- abs(d-t.date)
      ind <- which(diff <= temporal.buffer)
      if (length(ind)>0) {as.numeric(d[ind[which(diff[ind]==min(diff[ind]))]])
      } else {return(NA)}} else {return(NA)}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. build composites
#-------------------------------------------------------------------------------------------------------------------------------#

  r.value <- calc(img, f1)
  r.date <- calc(img, f2)

#-------------------------------------------------------------------------------------------------------------------------------#
# 5. build output
#-------------------------------------------------------------------------------------------------------------------------------#

  # build table from date info
  ud <- unique(r.date)
  used <- img.dates[as.numeric(img.dates)%in%ud]
  ud <- as.numeric(img.dates)
  count <- sapply(ud, function(x) {cellStats(r.date==x, sum)})
  df <- data.frame(date=used, code=ud, count=count)

  # check for missing values
  df <- rbind(df, c(NA, NA, cellStats(is.na(r.date), sum)))

  # return data
  return(list(composite=r.value, dates=r.date, pixel.count=df, target.date=t.date, temporal.buffer=temporal.buffer))

}
