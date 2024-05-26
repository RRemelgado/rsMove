#' @title sMoveRes
#'
#' @description Analysis of GPS data losses with the choice spatial resolutions.
#' @param x Object of class \emph{spatVector}.
#' @param y Numeric vector with target  spatial resolutions.
#' @importFrom terra ext crds extend
#' @importFrom dbscan dbscan
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes theme geom_bar coord_cartesian element_blank element_line scale_y_continuous scale_fill_gradientn labs
#' @return A \emph{list}.
#' @details {The function simulates how many unique combinations of animal
#' observations and environmental data would be preserved when choosing an
#' environmental dataset with the resolution given in \emph{y}. It also
#' estimates simulates how many unique pixel regions (i.e., groups of
#' spatially connected pixels) would be preserved after accounting for
#' pseudo-replication. Finally, For each spatial resolution, the function
#' reports the 'pixel ratio' and the 'region ratio', i.e., the number of
#' pixels and regions divided by the number of GPS records.}
#' @return {The function returns a list with:
#' \itemize{
#'  \item{\emph{stats} - Summary statistics reporting on the number
#'  of unique samples and sample regions per spatial resolution.}
#'  \item{\emph{plot} - Plot representing the change in number of
#'  samples and sample regions per spatial resolution.}}}
#' @seealso \code{\link{tMoveRes}}
#' @examples {
#'
#'  require(terra)
#'
#'  # read movement data
#'  shortMove = read.csv(system.file('extdata',
#'  'shortMove.csv', package="rsMove"))
#'
#'  # convert observations to vector
#'  shortMove = vect(shortMove, geom=c("x","y"), crs="EPSG:32632")
#'
#'  # test function for 5, 10 20 and 30 m
#'  a.res = sMoveRes(shortMove, c(5, 10, 20, 30))
#'
#' }
#' @export

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

sMoveRes = function(x, y) {

  #---------------------------------------------------------------------------#
  #  1. check inpur variables
  #---------------------------------------------------------------------------#

  if (!class(x)[1]%in%c('SpatVector')) {stop('"x" is not of a valid class')}
  if (!is.numeric(y)) {stop('"y" is not numeric')}
  if (!is.vector(y)) {stop('"y" is not a vector')}

  #---------------------------------------------------------------------------#
  # 2. find unique sample regions
  #---------------------------------------------------------------------------#

  out = do.call(rbind, lapply(y, function(r) {

    # reference raster (extend to avoid missing samples along the borders)
    e = extend(rast(ext(x), res=r, crs=crs(x)), r)

    # cell positions of elements x and region clustering
    odf = data.frame(resolution=r,
                     nr.pixels=length(unique(cellFromXY(e, crds(x)))),
                     nr.regions=length(unique(dbscan(crds(x),
                                                 eps=r,borderPoints=T,
                                                 minPts=1)$cluster)))

    # output data frame with statistics
    return(odf)

  }))

  # calculate the ratio of pixes/regions per number of entries in x
  out$nr.observations = nrow(x)
  out$pixel.ratio = out$nr.pixels / out$nr.observations
  out$region.ratio = out$nr.regions / out$nr.observations

  #---------------------------------------------------------------------------#
  # 3. plot output
  #---------------------------------------------------------------------------#

  # make color palette
  # make color palette
  cr = colorRampPalette(c('#8c510a','#bf812d','#dfc27d',
                          '#f6e8c3','#f5f5f5','#c7eae5',
                          '#80cdc1','#35978f','#01665e'))

  # build plot object
  out$resolution = factor(out$resolution, levels=y)
  p = ggplot(out, aes(x=resolution, y=pixel.ratio, fill=region.ratio)) +
    theme_bw(base_size=6) + geom_bar(stat="identity", width=1) +
    scale_y_continuous(expand=c(0,0)) + coord_cartesian(ylim=c(0,1.0)) +
    scale_fill_gradientn(colors=cr(9), breaks=c(0.0, 0.5, 1.0),
    limits=c(0,1.0), name="Region ratio\n") +
    labs(x="\nResolution", y="Pixel ratio\n") +
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(linewidth=0.2, colour="grey5"))

  # return data frame and plot
  return(list(stats=out, plot=p))

}

