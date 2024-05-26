#' @title moveCloud
#'
#' @description Extract Cloud Cover Fraction (CCF) data for a set of coordinate pairs.
#' @param x Object of class \emph{SpatVector}.
#' @param y Object of class \emph{Date} with observation dates of \emph{y}.
#' @param start First data from when to download CCF data.
#' @param end Last data from when to download CCF data.
#' @param interval Daily interval of data download.
#' @param data.path Output path for downloading data on cloud cover.
#' @importFrom terra rast extract app crds
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes labs theme geom_bar coord_cartesian scale_y_continuous scale_fill_gradientn
#' @importFrom utils download.file
#' @importFrom RCurl url.exists
#' @importFrom stats  weighted.mean
#' @details {The function extracts data on daily Cloud Cover Fractions (CCF)
#' from NASA's Earth Observation (NEO). For a sequence of dates defined by
#' \emph{start}, \emph{end}, and \emph{interval}, the function downloads the
#' correspondent CCF data and stores them in  \emph{data.path}. These data,
#' which have a global coverage and a spatial resolution of 0.1 degrees, will
#' only be downloaded if they do not already exist in \emph{data.path}. When
#' checking for existing data, the function will follow a standard naming
#' convention, so data acquired independently will likely be missed. After
#' downloading the needed data, the function will extract the CCF fraction
#' values for all dates at the coordinates in \emph{x}. This will be used to
#' calculate the mean and standard deviation of the CCF at each date. In
#' addition, the function will extract CCF data for each date in \emph{y} that
#' falls between \emph{start} and \emph{end}.}
#' @return {A \emph{list} containing: \itemize{
#'  \item{\emph{x.stats} - CCF statistics at each entry in \emph{x}}
#'  \item{\emph{r.stats} - CCF statistics at each unique date between\emph{start} and \emph{end} summarized across the elements of \emph{x}}
#'  \item{\emph{x.plot} - Plot of the data in \emph{x.stats}}.
#'  \item{\emph{r.plot} - Plot of the data in \emph{r.stats}}.}}
#' @references \url{https://cneos.jpl.nasa.gov/}
#' @seealso \code{\link{sMoveRes}} \code{\link{tMoveRes}}
#' @examples \dontrun{
#'
#'  require(terra)
#'
#'  # read movement data
#'  longMove <- read.csv(system.file('extdata',
#'  'longMove.csv', package="rsMove"))
#'
#'  # convert observations to vector
#'  longMove = vect(longMove, geom=c("long","lat"), crs="EPSG:4326")
#'
#'  # test function for 30 day buffer
#'  obs.dates <- as.Date(longMove$timestamp)
#'  c.cover <- moveCloud(shortMove, obs.dates)
#'
#' }
#' @export

#-----------------------------------------------------------------------------------------------------------------------#

moveCloud <- function(x, y, start, end, interval, data.path) {

  #---------------------------------------------------------------------------------------------------------------------#
  #  1. check inpur variables
  #---------------------------------------------------------------------------------------------------------------------#

  if (!class(x)[1]%in%c('SpatVector')) stop('"x" is not of a valid class')
  x = as.data.frame(crds(project(x, "EPSG:4326"))) # reproject data if needed
  if (class(y)[1] != 'Date') stop('"start" is not of a valid class')

  if (class(start)[1] != 'Date') stop('"start" is not of a valid class')
  if (class(end)[1] != 'Date') stop('"end" is not of a valid class')
  if (!is.numeric(interval)) stop('"interval" is not of a valid class')

  if (!dir.exists(data.path)) {
    test = tryCatch(dir.create(data.path), error=function(e) return(TRUE))
    if (is.logical(test)) stop('provided "data.path" not valid')
  }

  # data sources
  myd <- "https://neo.gsfc.nasa.gov/archive/geotiff.float/MYDAL2_D_CLD_FR/" # aqua
  mod <- "https://neo.gsfc.nasa.gov/archive/geotiff.float/MODAL2_D_CLD_FR//" # terra

  #---------------------------------------------------------------------------#
  # 2. download data on cloud cover
  #---------------------------------------------------------------------------#

  # target dates
  target_dates = seq(start, end, interval)

  # warn user of the volume of data required
  data.volume = 5.6*length(target_dates)
  message(paste0("warning: about to store ", data.volume, " Mb"))

  # download data
  dates = list()
  files = list()

  for (d in 1:length(target_dates)) {

    date = target_dates[d]

    # file to store
    ofile = file.path(data.path, paste0("MYDAL2-ccf_",
                                        paste0(strsplit(as.character(date),
                                                        "[-]")[[1]],
                                               collapse=""), "_10km.tif"))

    if (!file.exists(ofile)) {

      # list images to be downloaded
      s1 = file.path(myd, paste0("MYDAL2_D_CLD_FR_",
                                 date, ".FLOAT.TIFF"))
      o1 = file.path(data.path, basename(s1))
      if (!file.exists(s1)) c1 = tryCatch(
        download.file(s1, o1, mode="wb", quiet=T),
        error=function(e) return(TRUE))

      s2 = file.path(mod, paste0("MODAL2_D_CLD_FR_",
                                 date, ".FLOAT.TIFF"))
      o2 = file.path(data.path, basename(s2))
      if (!file.exists(s2)) c2 = tryCatch(
        download.file(s2, o2, mode="wb", quiet=T),
        error=function(e) return(TRUE))

      if (!is.logical(c1) | !is.logical(c2)) {

        tmp = c(o1, o2)
        tmp = tmp[file.exists(tmp)]
        cloud.cover.img = rast(tmp)
        cloud.cover.img[cloud.cover.img > 1] = NA

        # average cloud cover measured with AQUA & TERRA
        cloud.cover.img = app(cloud.cover.img, mean, na.rm=T)
        files[[d]] = ofile
        writeRaster(cloud.cover.img, files[[d]])
        dates[[d]] = date

        file.remove(tmp)

      }

    } else {
      files[[d]] = ofile
      dates[[d]] = date
    }

  }

  files = do.call("c", files)
  dates = do.call("c", dates)

  #---------------------------------------------------------------------------#
  # 3. extract cloud cover values for each GPS coordinate
  #---------------------------------------------------------------------------#

  # find GPS coordinates falling within the specified temporal window
  ind = which(y %in% dates)

  r.cloud.cover = extract(rast(files), x[ind,], ID=F)
  r.cloud.cover.data = data.frame(date=dates,
                                  mean=apply(r.cloud.cover, 2, mean, na.rm=T),
                                  sd=apply(r.cloud.cover, 2, sd, na.rm=T))

  x.cloud.cover = rep(0,length(ind))
  unique_dates = unique(y[ind])
  for (d in 1:length(unique_dates)) {
    di = which(y[ind] == unique_dates[d])
    ii = which(dates == unique_dates[d])
    e = extract(rast(files[ii]), x[di,], ID=F)[[1]]
    x.cloud.cover[di] = e
  }

  x.cloud.cover.data = data.frame(index=ind, date=y[ind],
                                  cloud.cover=x.cloud.cover)

  #---------------------------------------------------------------------------#
  # 4. build plot
  #---------------------------------------------------------------------------#

  # table used to plot
  gdf = ddply(x.cloud.cover.data, .(date), summarise,
              cover=mean(cloud.cover),
              mean=mean(cloud.cover),
              sd=sd(cloud.cover))

  p1 <- ggplot(gdf, aes(x=date, y=cover)) +
    theme_bw(base_size=6) +
    geom_bar(width=0.7, stat="identity", fill="grey60") +
    coord_cartesian(xlim=c(start,end), ylim=c(0,1)) +
    labs(x="Observation date", y="Cloud Cover Fraction") +
    scale_y_continuous(expand=c(0,0)) +
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(linewidth=0.2, colour="grey5"))

  # plot cloud cover across the reference period
  p2 <- ggplot(r.cloud.cover.data, aes(x=date, y=mean)) +
    theme_bw(base_size=6) +
    geom_line(colour="grey10") +
    geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha=0.3) +
    coord_cartesian(xlim=c(start,end), ylim=c(0,1)) +
    labs(x="Observation date", y="Cloud Cover Fraction") +
    scale_y_continuous(expand=c(0,0)) +
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(linewidth=0.2, colour="grey5"))

  return(list(x.stats=x.cloud.cover.data,
              r.stats=r.cloud.cover.data,
              x.plot=p1, r.plot=p2))

}
