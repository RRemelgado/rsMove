#' @title moveCloud
#'
#' @description Provides historical information on cloud cover.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param d.path Output data path for downloaded data.
#' @param d.buffer
#' @param remove.files
#' @param p.res Should the output be ploted on screen? Default is TRUE.
#' @import ggplot2 sp rgdal grDevices
#' @importFrom utils download.file
#' @return A \emph{list}.
#' @details {This function makes uses daily cloud fraction data from NASA's NEO service. 
#' For each observation date, the function downloads the correspondent image and extracts 
#' the percent of cloud cover for the samples acquired at the target date. If \emph{d.path} 
#' is specified, the function will look within the provided directory if the required files 
#' are already. If so, they won't be downloaded. If \emph{d.buffer} is specified, for each 
#' date, the function will consider images before and after within a temporal buffer. These 
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
#'  # test function for 5, 10 20 and 30 m
#'  a.res <- tMoveRes(xy=moveData, dpath='.')
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

moveCloud <- function(xy=xy, d.path=NULL, b.size=NULL, remove.file=TRUE, p.res=T) {
  
#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#
  
  # input keywords
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  if (is.na(rr@projargs)) {stop('"xy" does not have a valid projection')}
  if (is.null(d.path)) {
    d.path <- tempdir()
    remove.file <- TRUE
  } else {
    if (!dir.exists(d.path)) {stop('"dpath" not found in file system')}
    d.path <- paste0(file.path(d.path), .Platform$file.sep)}
  if (!is.null(b.size)) {apply.buffer=TRUE}
  if (!is.logical(p.res)) {stop('"p.res" is not a logical argument')}
  
  # ftp servers
  myd <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MYDAL2_D_CLD_FR/" # aqua
  mod <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MODAL2_D_CLD_FR/" # terra
  
#---------------------------------------------------------------------------------------------------------------------#
# 2. download data and derive statistics
#---------------------------------------------------------------------------------------------------------------------#
  
  # target dates
  ot <- as.Date(ot)
  ud <- unique(ot)
  
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
    loc <- which(ot==ud[d])
    
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
      day.ls <- seq(ud[d]-b.size, ud[d]+b.size, 1)
      
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
        bv <- min(f.cc[x,ind])
        diff1 <- abs(diff0[ind])
        ind <- ind[which(f.cc[x,ind]==av)]
        ind <- ind[which(diff0[ind]==min(diff0[ind]))]
        db <- diff0[ind]
        bd <- day.ls[ind]
        ind <- which(diff0 > 0)
        bv <- min(f.cc[x,ind])
        diff1 <- abs(diff0[ind])
        ind <- ind[which(f.cc[x,ind]==av)]
        ind <- ind[which(diff0[ind]==min(diff0[ind]))]
        da <- abs(day.ls[ind]-ud[d])
        ad <- day.ls[ind]
        return(list(bd=bd, bv=bv, db=db, ad=ad, av=av, da=da))})
      
      # update target variables
      p.dt.b[loc] <- do.call('c', lapply(dq, function(x) {x$bd}))
      p.cc.b[loc] <- unlist(lapply(dq, function(x) {x$bv}))
      d.df.b[loc] <- unlist(lapply(dq, function(x) {x$db}))
      p.dt.a[loc] <- do.call('c', lapply(dq, function(x) {x$ad}))
      p.cc.a[loc] <- unlist(lapply(dq, function(x) {x$av}))
      d.df.a[loc] <- unlist(lapply(dq, function(x) {x$da}))
      
      rm(f.cc)
      
    }
    
    # remove files if required
    if (remove.file) {file.remove(list.files(d.path, '_D_CLD_FR_'))}
    
  }
  
  # add column names output
  df <- data.frame(day.cover=d.cc, p.cover.before=p.cc.b, p.cover.after=p.cc.a, 
                   p.day.before=p.dt.b, p.day.after=p.dt.a, stringsAsFactors=F)
  
#-------------------------------------------------------------------------------------------#
# 3. build plot
#-------------------------------------------------------------------------------------------#
  
  # make color ramp
  cr = colorRampPalette(c("forestgreen", "khaki2", "darkred"))
  
  # build plot
  odf <- data.frame(dd=d.df.b, cc=p.cc.b)
  p1 <- ggplot(odf, aes(x=dd, fill=cc)) + geom_histogram(binwidth=1) + 
    scale_x_continuous(limits=c(-b.size, 0)) + 
    scale_fill_gradientn(name='Cloud cover %', colors=cr(10), limits=c(0,100), breaks=c(0, 50, 100)) +
    xlab("Day Difference") + ylab("Number of Samples") + 
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12))
  
  # determine optimal y range
  mv <- max(ggplot_build(p1)$data[[1]]$count)
  nc <- nchar(as.character(mv))
  m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
  mv <- mv / m
  yr <- round(mv)
  if (mv > yr) {yr <- (yr+0.05)*m} else {yr <- yr*m}
  
  # update plot
  p1 <- p1 + scale_y_continuous(limits=c(0, yr))
  
  # build plot
  odf <- data.frame(dd=d.df.a, cc=p.cc.a)
  p2 <- ggplot(odf, aes(x=dd, fill=cc)) + geom_histogram(binwidth=1) + 
    scale_x_continuous(limits=c(0, b.size)) + 
    scale_fill_gradientn(name='Cloud cover %', colors=cr(10), limits=c(0,100), breaks=c(0, 50, 100)) +
    xlab("Day Difference") + ylab("Number of Samples") + 
    theme(axis.text=element_text(size=10), 
          axis.title=element_text(size=12))
  
  # determine optimal y range
  mv <- max(ggplot_build(p2)$data[[1]]$count)
  nc <- nchar(as.character(mv))
  m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
  mv <- mv / m
  yr <- round(mv)
  if (mv > yr) {yr <- (yr+0.05)*m} else {yr <- yr*m}
  
  # update plot
  p2 <- p2 + scale_y_continuous(limits=c(0, yr))
  
  if (p.res) {p} # plot raster on screen
  
  #---------------------------------------------------------------------------------------------------------------------#
  #  7. derive output
  #---------------------------------------------------------------------------------------------------------------------#
  
  # return data frame and plot
  return(list(stats=df, plot.data=odf, plot.before=p1, plot.after=p2))
  
}