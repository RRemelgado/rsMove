#' @title moveCloud
#'
#' @description {Provides historical information on cloud cover for a set of coordinate
#' pairs. The temporal information is adjusted to the sample observation date}.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.dates Object of class \emph{Date} with observation dates of \emph{xy}.
#' @param data.path Output data path for downloaded data.
#' @param buffer.size Two element vector with temporal buffer size (expressed in days).
#' @param remove.file Logical. Should the files be deleted after usage?
#' @importFrom raster raster extract
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot xlab ylab theme geom_bar
#' @importFrom utils download.file
#' @importFrom RCurl url.exists
#' @return A \emph{list} object reporting on the variability of cloud cover within and around each observation dates.
#' @details {This function uses daily cloud fraction data from NASA's NEO service.
#' For each observation date in \emph{obs.dates}, the function downloads the correspondent image
#' and extracts the percent of cloud cover for the corresponding samples in \emph{xy}. If
#' \emph{data.path} is specified, the function will look within the provided directory for the
#' cloud cover images. If they exist, they won't be downloaded reducing the amount of time required
#' by the function. Moreover, if \emph{buffer.size} is specified, for each date, the function will only
#' consider images that are within the specified temporal buffer. \emph{buffer.size} requires a two
#' element vector which specifies the buffer size before and after the target dates. These additional
#' images will be used to report on the closest time step with the lowest possible cloud cover. The
#' final output provides a \emph{data.frame} ($report) with information on:
#' \itemize{
#'  \item{\emph{cloud cover \% (day)}: cloud cover for the observation dates.}
#'  \item{\emph{best date (after)}: dates before the obsertation dates with the lowest cloud cover.}
#'  \item{\emph{best date cloud cover \% (before)}: cloud cover for best before dates.}
#'  \item{\emph{best date (after)}: dates after the obsertation dates with the lowest cloud cover.}
#'  \item{\emph{best date cloud cover \% (after)}: cloud cover best after dates.}}
#'  Finally, the function generates a plot ($plot) reporting on the variability of cloud cover
#'  within the dates provided by \emph{obs.dates} and the number of samples registered within them.}
#' @references \url{https://cneos.jpl.nasa.gov/}
#' @seealso \code{\link{sMoveRes}} \code{\link{tMoveRes}}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(r))
#'
#'  # test function for 30 day buffer
#'  od <- as.Date(moveData@data$date)
#'  c.cover <- moveCloud(xy=moveData, obs.dates=od, data.path=".", buffer.size=c(30,30))
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

moveCloud <- function(xy=xy, obs.dates=obs.dates, data.path=NULL, buffer.size=NULL, remove.file=FALSE) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  # input keywords
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.na(crs(xy))) {stop('"xy" does not have a valid projection')}
  if (is.null(data.path)) {
    data.path <- tempdir()
    remove.file <- TRUE
  } else {
    if (!dir.exists(data.path)) {stop('"dpath" not found in file system')}
    data.path <- paste0(file.path(data.path), .Platform$file.sep)}
  if (!is.null(buffer.size)) {apply.buffer<-TRUE} else {apply.buffer<-FALSE}

  # ftp servers
  myd <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MYDAL2_D_CLD_FR/" # aqua
  mod <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MODAL2_D_CLD_FR/" # terra

#---------------------------------------------------------------------------------------------------------------------#
# 2. download data and derive statistics
#---------------------------------------------------------------------------------------------------------------------#

  # target dates
  obs.dates <- as.Date(obs.dates)
  ud <- unique(obs.dates)

  # output variables
  d.cc <- vector('numeric', length(xy))
  p.cc.b <- d.cc
  p.dt.b <- d.cc
  class(p.dt.b) <- "Date"
  d.df.b <- d.cc
  p.dt.df <- d.cc
  p.cc.a <- d.cc
  p.dt.a <- p.dt.b
  d.df.a <- d.cc

  for (d in 1:length(ud)) {

    # target observations
    loc <- which(obs.dates==ud[d])

    # set file name
    ifile1 <- paste0(mod, "MODAL2_D_CLD_FR_", ud[d], ".FLOAT.TIFF")
    ofile1 <- paste0(data.path, basename(ifile1))
    ifile2 <- paste0(myd, "MYDAL2_D_CLD_FR_", ud[d], ".FLOAT.TIFF")
    ofile2 <- paste0(data.path, basename(ifile2))

    # check if file exists
    if (!file.exists(ofile1)) {
      if (url.exists(ifile1)) {download.file(ifile1, ofile1, quiet=TRUE, mode="wb")
        mod.r <- TRUE} else {mod.r <- FALSE}} else {mod.r <- TRUE}
    if (!file.exists(ofile2)) {
      if (url.exists(ifile2)) {download.file(ifile2, ofile2, quiet=TRUE, mode="wb")
        myd.r <- TRUE} else {myd.r <- FALSE}} else {myd.r <- TRUE}

    # read data and crop to xy extent
    if (mod.r & myd.r) {d.cc[loc] <- (extract(raster(ofile1), xy[loc,]) + extract(raster(ofile2), xy[loc,])) / 2}
    if (mod.r & !myd.r) {d.cc[loc] <- extract(raster(ofile1), xy[loc,])}
    if (!mod.r & myd.r) {d.cc[loc] <- extract(raster(ofile1), xy[loc,])}

    # search for nearby images
    if(apply.buffer) {

      # determine dates within the buffer
      day.ls <- seq(ud[d]-buffer.size[1], ud[d]+buffer.size[2], 1)

      df <- lapply(day.ls, function(x) {
        ifile1 <- paste0(mod, "MODAL2_D_CLD_FR_", x, ".FLOAT.TIFF")
        ofile1 <- paste0(data.path, basename(ifile1))
        ifile2 <- paste0(myd, "MYDAL2_D_CLD_FR_", x, ".FLOAT.TIFF")
        ofile2 <- paste0(data.path, basename(ifile2))
        if (!file.exists(ofile1)) {
          if (url.exists(ifile1)) {download.file(ifile1, ofile1, quiet=TRUE, mode="wb")
            mod.r <- TRUE} else {mod.r <- FALSE}}
        if (!file.exists(ofile2)) {
          if (url.exists(ifile2)) {download.file(ifile2, ofile2, quiet=TRUE, mode="wb")
            myd.r <- TRUE} else {myd.r <- FALSE}}
        if (mod.r & myd.r) {return((extract(raster(ofile1), xy[loc,]) +
                                     extract(raster(ofile2), xy[loc,])) / 2)}
        if (mod.r & !myd.r) {return(extract(raster(ofile1), xy[loc,]))}
        if (!mod.r & myd.r) {return(extract(raster(ofile2), xy[loc,]))}})

      # extract values
      f.cc <- do.call(cbind, df)

      # find closest minimum
      dq <- lapply(1:length(loc), function(x) {
        diff0 <- day.ls - ud[d]
        ind <- which(diff0 < 0)
        if (length(ind)>0) {
          bv <- min(f.cc[x,ind])
          diff1 <- abs(diff0[ind])
          ind <- ind[which(f.cc[x,ind]==bv)]
          ind <- ind[which(diff0[ind]==min(diff0[ind]))]
          db <- diff0[ind]
          bd <- day.ls[ind]
        } else {
          bd <- NA
          bv <- NA
          db <- NA}
        ind <- which(diff0 > 0)
        if (length(ind)>0) {
          av <- min(f.cc[x,ind])
          diff1 <- abs(diff0[ind])
          ind <- ind[which(f.cc[x,ind]==av)]
          ind <- ind[which(diff0[ind]==min(diff0[ind]))]
          da <- abs(day.ls[ind]-ud[d])
          ad <- day.ls[ind]
        } else {
          ad <- NA
          av <- NA
          da <- NA}
        return(list(bd=bd, bv=bv, db=db, ad=ad, av=av, da=da))})

      # update target variables
      p.dt.b[loc] <- do.call('c', lapply(dq, function(x) {x$bd}))
      p.cc.b[loc] <- unlist(lapply(dq, function(x) {x$bv}))
      d.df.b[loc] <- unlist(lapply(dq, function(x) {x$db}))
      p.dt.a[loc] <- do.call('c', lapply(dq, function(x) {x$ad}))
      p.cc.a[loc] <- unlist(lapply(dq, function(x) {x$av}))
      d.df.a[loc] <- unlist(lapply(dq, function(x) {x$da}))

      rm(f.cc)

    } else {
      p.cc.b <- NA
      p.cc.a <- NA
      p.dt.b <- NA
      p.dt.a <- NA}

    # remove files if required
    if (remove.file) {file.remove(list.files(data.path, '_D_CLD_FR_'))}

  }

  # add column names output
  df <- data.frame(date=obs.dates, day.cover=d.cc, p.day.before=p.dt.b, p.cover.before=p.cc.b,
                   p.cover.after=p.cc.a, p.day.after=p.dt.a, stringsAsFactors=F)

  colnames(df) <- c("date (original)", "cloud cover % (day)", "best date (before)",
                    "best date cloud cover % (before)", "best date (after)",
                    "best date cloud cover % (after)")

# #-------------------------------------------------------------------------------------------#
# # 3. build plot
# #-------------------------------------------------------------------------------------------#

  # plot table
  ud <- sort(ud)
  df0 <- do.call(rbind, lapply(ud, function(d) {
    ind <- which(obs.dates==d)
    cc <- mean(df[ind,2], na.rm=TRUE)
    data.frame(date=d, cover=cc, count=length(ind), stringsAsFactors=FALSE)}))

  # color ramp of fill
  cr <- colorRampPalette(c("khaki2", "forestgreen"))

  # plot
  p <- ggplot(df0, aes(x=date, y=cover, fill=count)) + theme_bw() +
    scale_fill_gradientn(colors=cr(10), name="Nr. Samples") +
    xlab("Observation dates") + ylab("Cloud cover (%)") +
    geom_bar(width=0.7, stat = "identity") +
    theme(axis.text.x=element_text(size=12),
          axis.title.x =element_text(size=14),
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14)) + ylim(0,yr)

#   p <- ggplot(df, aes_string(x="date (original)", y="cloud cover % (day)")) +
#     theme_bw() + geom_bar(stat="identity", colour="black", fill="grey80") +
#     theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
#     xlab("Date") + ylab("Cloud Cover (%)")

  return(list(stats=df, plot=p))

}
