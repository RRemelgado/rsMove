#' @title tMoveRes
#'
#' @description Analysis of GPS data losses with the choice temporal resolutions.
#' @param x \emph{spatVector}.
#' @param y \emph{Date} vector with observation dates or each entry in \emph{x}.
#' @param time.res Vector of temporal resolutions (expressed in days).
#' @param pixel.res Spatial resolution (unit depends on spatial projection).
#' @importFrom ggplot2 ggplot aes theme geom_bar coord_cartesian element_blank element_line scale_y_continuous scale_fill_gradientn labs
#' @importFrom terra rast ext extend cellFromXY crs
#' @importFrom plyr ddply . summarise
#' @importFrom dbscan dbscan
#' @importFrom grDevices colorRampPalette
#' @details {For each spatial resolution given by  \emph{pixel.res}, and for
#' each temporal resolutions given by \emph{time.res}, the function simulates
#' the number of unique pixels and pixel regions preserved after accounting
#' for pseudo-replication at each hypothetical time step, assuming that the
#' GPS data until the next date with available environmental data is
#' aggregated into unique pixels. Finally, For each temporal aggregation
#' window, the function reports the 'pixel ratio' and the 'region ratio',
#' i.e., the number of pixels and regions divided by the number of GPS records.}
#' @return {A \emph{list} containing:
#' \itemize{
#'  \item{\emph{stats} - Summary statistics reporting on
#'  the number of temporal widows, unique samples and unique
#'  sample regions per temporal resolution.}
#'  \item{\emph{summary.plot} - Plot representing the change in number
#'  of unique pixels and pixel regions per temporal resolution.}
#'  \item{\emph{temporal.plot} - Plot representing the change in number
#'  of unique pixels and pixel regions per temporal resolution and time step.}}}
#' @seealso \code{\link{sMoveRes}}
#' @examples {
#'
#'  require(terra)
#'
#'  #'  # read movement data
#'  longMove = read.csv(system.file('extdata',
#'  'longMove.csv', package="rsMove"))
#'
#'  # convert observations to vector
#'  longMove = vect(longMove, geom=c("long","lat"), crs="EPSG:4326")
#'
#'  # test function for intervals of 1, 8 and 16 days (e.g. of MODIS data)
#'  obs.date = as.Date(longMove$timestamp)
#'  a.res = tMoveRes(longMove, obs.date, c(1,8,16), 0.1)
#'
#' }
#' @export

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

tMoveRes = function(x, y, time.res, pixel.res) {

  #---------------------------------------------------------------------------#
  #  1. check input variables
  #---------------------------------------------------------------------------#

  if (!class(x)[1]%in%c('SpatVector')) stop('"x" is not of a valid class')
  if (!class(y)[1]%in%c('Date')) stop('"y" is not of a valid class')
  if (length(pixel.res)>1) {stop('"pixel.res" has more than one element')}
  if (!is.numeric(time.res)) {stop('"time.res" is not numeric')}

  #---------------------------------------------------------------------------#
  # 2. determine pixel aggregations
  #---------------------------------------------------------------------------#

  st = min(y) # start time
  et = max(y) # end time

  # reference raster (extend to avoid missing samples along the borders)
  e = extend(rast(ext(x), res=pixel.res, crs=crs(x)), pixel.res)

  out = do.call(rbind, lapply(time.res, function(r) {

    nw = round(as.numeric((et - st)) / r + 1) # number of temporal windows

    tmp = do.call(rbind, lapply(1:nw, function(w) {

      # last date
      date = ((st+r)+(r*w))-1

      # GPS observations recorded within the target window
      loc = which(y >= (st+r*(w-1)) & y <= (((st+r)+(r*w))-1))

      if (length(loc) > 0) {

        # cell positions of elements x and region clustering
        odf = data.frame(date=date, resolution=r, nr.observations=length(loc),
                         nr.pixels=length(unique(cellFromXY(e, crds(x[loc,])))),
                         nr.regions=length(unique(dbscan(crds(x[loc,]),
                                                         eps=r,borderPoints=T,
                                                         minPts=1)$cluster)))

        return(odf)
      } else {
        return(data.frame(date=date, resolution=r,
                          nr.observations=0, nr.pixels=0,
                          nr.regions=0))
      }

    }))

    #  estimate final count of pixels/regions
    return(tmp)

  }))

  # calculate the ratio of pixes/regions per number of entries in x
  i = which(out$nr.observations > 0)
  out$pixel.ratio = 0
  out$region.ratio = 0
  out$pixel.ratio[i] = out$nr.pixels[i] / out$nr.observations[i]
  out$region.ratio[i] = out$nr.regions[i] / out$nr.observations[i]

  #---------------------------------------------------------------------------#
  # 3. plot output
  #---------------------------------------------------------------------------#

  # make color palette
  cr = colorRampPalette(c('#8c510a','#bf812d','#dfc27d',
                        '#f6e8c3','#f5f5f5','#c7eae5',
                        '#80cdc1','#35978f','#01665e'))

  # build plot object
  out$resolution = factor(out$resolution, levels=sort(time.res))

  gdf = ddply(out, .(resolution), summarise,
              pixel.ratio=sum(nr.pixels)/sum(nr.observations),
              region.ratio=sum(nr.regions)/sum(nr.observations))

  p1 = ggplot(gdf, aes(x=resolution, y=pixel.ratio, fill=region.ratio)) +
    theme_bw(base_size=6) + geom_bar(stat="identity", width=1) +
    scale_y_continuous(expand=c(0,0)) + coord_cartesian(ylim=c(0,1.0)) +
    scale_fill_gradientn(colors=cr(9), breaks=c(0.0, 0.5, 1.0),
                         limits=c(0,1.0), name="Region ratio\n") +
    labs(x="\nResolution (m)", y="Pixel ratio\n") +
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(linewidth=0.2, colour="grey5"))

  p2 = ggplot(out, aes(x=date, y=pixel.ratio)) +
    theme_bw(base_size=6) +
    geom_bar(stat="identity", width=1, fill="grey60") +
    scale_y_continuous(expand=c(0,0)) + coord_cartesian(ylim=c(0,1.0)) +
    labs(x="\nDate", y="Pixel ratio\n") +
    facet_wrap(~resolution, nrow=1) + # , strip.position="right") +
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          strip.background=element_blank(),
          axis.line=element_line(linewidth=0.2, colour="grey5"))

  # return data frame and plot
  return(list(stats=out, summary.plot=p1, temporal.plot=p2))

}
