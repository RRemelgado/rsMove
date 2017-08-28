#' @title moveCloud
#'
#' @description Provides historical information on cloud cover.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param o.time Object of class \emph{Date}.
#' @param d.path Output data path for downloaded data.
#' @param b.size Two element vector with temporal buffer size (expressed in days).
#' @param remove.file Logical. Should the files be deleted after usage?
#' @import ggplot2 sp rgdal grDevices
#' @importFrom utils download.file
#' @importFrom RCurl url.exists
#' @return A \emph{list}.
#' @details {This function makes uses daily cloud fraction data from NASA's NEO service.
#' For each observation date (\emph{o.time}), the function downloads the correspondent image
#' and extracts the percent of cloud cover for the samples acquired at the target date. If
#' \emph{d.path} is specified, the function will look within the provided directory for the
#' required files. If so, they won't be downloaded. If \emph{d.buffer} is specified, for each
#' date, the function will consider images before and after within a temporal buffer. \emph{d.buffer}
#' requires two elements which specify the buffer size before and after the target date. These
#' new images will be used to report on the closest time step with the lowest possible cloud
#' cover. The final report provides information on:
#' \itemize{
#'  \item{\emph{day.cover}: cloud cover for the observation date}
#'  \item{\emph{p.day.before}: date before the obsertation date with the lowest cloud cover}
#'  \item{\emph{p.cover.before}: cloud cover for \emph{p.day.before}}
#'  \item{\emph{p.day.after}: date after the obsertation date with the lowest cloud cover}
#'  \item{\emph{p.cover.after}: cloud cover for \emph{p.day.after}}}
#'  The output will also contain a two plots with information on the distance between the
#'  observation dates and the closest date with the lowest cloud cover. The plots show the
#'  amount of samples that are covered in each of the target dates for the best dates before
#'  (\emph{$plot.before}) and after (\emph{$plot.after}) the observation dates.}
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
#'  c.cover <- moveCloud(xy=moveData, o.time=od, d.path=".", b.size=c(30,30))
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

moveCloud <- function(xy=xy, o.time=o.time, d.path=NULL, b.size=NULL, remove.file=FALSE) {

#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#

  # input keywords
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.na(crs(xy))) {stop('"xy" does not have a valid projection')}
  if (is.null(d.path)) {
    d.path <- tempdir()
    remove.file <- TRUE
  } else {
    if (!dir.exists(d.path)) {stop('"dpath" not found in file system')}
    d.path <- paste0(file.path(d.path), .Platform$file.sep)}
  if (!is.null(b.size)) {apply.buffer<-TRUE} else {apply.buffer<-FALSE}

  # ftp servers
  myd <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MYDAL2_D_CLD_FR/" # aqua
  mod <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MODAL2_D_CLD_FR/" # terra

#---------------------------------------------------------------------------------------------------------------------#
# 2. download data and derive statistics
#---------------------------------------------------------------------------------------------------------------------#

  # target dates
  o.time <- as.Date(o.time)
  ud <- unique(o.time)

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

  # table that will contain full data frame of cloud cover
  if (apply.buffer) {o.cc <- vector('list', length(ud))}

  for (d in 1:length(ud)) {

    # target observations
    loc <- which(o.time==ud[d])

    # set file name
    ifile1 <- paste0(mod, "MODAL2_D_CLD_FR_", ud[d], ".FLOAT.TIFF")
    ofile1 <- paste0(d.path, basename(ifile1))
    ifile2 <- paste0(myd, "MYDAL2_D_CLD_FR_", ud[d], ".FLOAT.TIFF")
    ofile2 <- paste0(d.path, basename(ifile2))

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
      day.ls <- seq(ud[d]-b.size[1], ud[d]+b.size[2], 1)

      df <- lapply(day.ls, function(x) {
        ifile1 <- paste0(mod, "MODAL2_D_CLD_FR_", x, ".FLOAT.TIFF")
        ofile1 <- paste0(d.path, basename(ifile1))
        ifile2 <- paste0(myd, "MYDAL2_D_CLD_FR_", x, ".FLOAT.TIFF")
        ofile2 <- paste0(d.path, basename(ifile2))
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

      # raw cloud cover information
      f.cc[f.cc > 1] <- NA
      o.cc[[d]] <- data.frame(value=apply(f.cc, 2, mean, na.rm=T), date=day.ls)

      rm(f.cc)

    } else {
      p.cc.b <- NA
      p.cc.a <- NA
      p.dt.b <- NA
      p.dt.a <- NA}

    # remove files if required
    if (remove.file) {file.remove(list.files(d.path, '_D_CLD_FR_'))}

  }

  # add column names output
  df <- data.frame(day.cover=d.cc, p.cover.before=p.cc.b, p.cover.after=p.cc.a,
                   p.day.before=p.dt.b, p.day.after=p.dt.a, stringsAsFactors=F)

#-------------------------------------------------------------------------------------------#
# 3. build plot
#-------------------------------------------------------------------------------------------#

  # build plot with extended cloud cover information
  if (apply.buffer) {

    df0 <- do.call(rbind, o.cc)
    ud <- unique(df0$date)
    df0 <- data.frame(date=ud, cover=sapply(ud, function(x) {df0$value[which(ud==x)]}))
    p <- ggplot(df0, aes_string(x="date", y="cover")) + theme_bw() +
      geom_bar(stat="identity", colour="black", fill="grey80") +
      theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
      xlab("Date") + ylab("Cloud Cover (%)")

    return(list(stats=df, daily.cover=df0, daily.cover.plot=p))

  } else {return(stats=df)}

}
