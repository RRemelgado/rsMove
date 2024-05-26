#' @title map_sites
#'
#' @description Detection of geographic regions of samples using a pixel based approach.
#' @param x Object of class \emph{SpatVector}.
#' @param y Unique individual identifier for each entry in \emph{x}.
#' @param z Timestamps of each element in \emph{x} given as a \emph{POSIXct} object.
#' @param resolution Maximum distance between data points to identify clusters.
#' @param min_size Minimum number of GPS data points per cluster.
#' @importFrom sp Polygon Polygons SpatialPolygons SpatialPolygonsDataFrame
#' @importFrom terra ext extend crs convHull crds
#' @importFrom plyr ddply summarise .
#' @importFrom dbscan dbscan
#' @details {The function provides three outputs. First, it labels each
#' entry in \emph{x} based on their spatial connectivity. Connections are
#' based on \emph{resolution}, which defines the maximum distance allowed
#' among members of a group of data points, below which the elements of
#' that group are treated as a region and labeled with the same unique,
#' non-zero, numeric identifier. Regions where the number of data points
#' is below \emph{min_size} are labelled with a 0. Second, for each region,
#' the function will return a polygon defined by the convex hull of the
#' composing data points. 0-labeled regions are excluded, as well as those
#' where the number of data points is less than 2. Third, for each region,
#' and for each sequence of days with data points within the target region,
#' we provide the following information and statistics:
#' \itemize{
#'  \item{\emph{region_id} - Region unique identifier.}
#'  \item{\emph{start.date} - Segment unique identifier.}
#'  \item{\emph{start.date} - First data date.}
#'  \item{\emph{end.date} - Last data date.}
#'  \item{\emph{day_cover_mean} - Mean percentage of the days in the segment with movement data}
#'  \item{\emph{day_cover_sd} - Standard deviation of the days in the segment with movement data.}
#'  \item{\emph{day_overlap_mean} - Mean percentage of the day across which movement data were recorded.}
#'  \item{\emph{day_overlap_sd} - Standard deviation of the percentage of the day across which movement data were recorded.}
#'  \item{\emph{nr_days} - Number of days with movement data.}
#'  \item{\emph{nr_individuals} - Number of individuals tracked.}
#'  \item{\emph{nr_records} - Number of movement data records.}
#'  \item{\emph{data_frequency_mean} - Mean frequency of GPS entries (in minutes).}
#'  \item{\emph{data_frequency_sd} - Standard deviation of the frequency of GPS entries (in minutes)}
#'  }
#'  Day segments are treated separately because they indicate that a given
#'  area was persistently occupied by one or more tracked individuals.}
#' @return {A \emph{list} containing:
#'  \itemize{
#'  \item{\emph{region.id} - Vector reporting on the region each element in \emph{x} belongs to.}
#'  \item{\emph{region.polygons} - Polygons for each temporal segment in each \emph{region.id.}
#'  \item{\emph{region.stats} - Statistics for each temporal segment.}
#'  }
#' }}
#' @examples {
#'
#' require(terra)
#'
#' # load samples
#' multiMove <- read.csv(system.file('extdata', 'multiMove.csv', package="rsMove"))
#'
#' # convert samples to vector
#' multiMove = vect(multiMove, geom=c("x","y"), crs="EPSG:4326")[1:10000,]
#'
#' # species identifier (only one individual, so set to 1)
#' species_id = multiMove$id
#'
#' # data-time of species observations
#' date = strptime(multiMove$timestamp,format="%Y-%m-%d %H:%M:%S", tz="GMT")
#'
#' # extract regions
#' hm <- map_sites(multiMove, species_id, date, 0.1)
#'
#' }
#' @export

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

map_sites <- function(x, y, z, resolution, min_size=1) {

  #---------------------------------------------------------------------------#
  #  1. check inpur variables
  #---------------------------------------------------------------------------#

  if (is.null(crs(x)))
    stop('"x" is missing a valid projection')
  if (!class(y) %in% c("numeric", "character"))
    stop('"y" must be a "Date" vector')
  if (!sum(class(z) %in% c("POSIXct","POSIXlt","POSIXt")) > 0)
    stop('"z" is not of a valid class')
  if (length(unique(c(nrow(x),length(y),length(z)))) > 1)
    stop('lengths of "x", "y", "z" do not match')
  if (!is.numeric(resolution))
    stop('"resolution" must be a numeric element')

  #---------------------------------------------------------------------------#
  # 2. label clusters of observations
  #---------------------------------------------------------------------------#

  regions <- dbscan(crds(x),
                    eps=resolution,
                    borderPoints=T,
                    minPts=min_size)$cluster

  #---------------------------------------------------------------------------#
  # 3. derive statistics for each regions
  #---------------------------------------------------------------------------#

  tmp = list()
  region_shp = list()

  dates = as.Date(format(z, "%Y-%m-%d"))

  for (r in unique(regions[which(regions > 0)])) {

    # evaluate temporal composition
    ri <- which(regions == r) # region indices

    # find unique dates with data at region r
    region_dates = sort(unique(dates[ri]))

    # between-day interval
    day_diff = c(1,diff(region_dates))

    # search for sequences of days
    segments = rle(day_diff)$lengths
    date_segments = vector("numeric", length(region_dates))
    for (p in 1:length(segments)) {
      date_segments[(sum(segments[0:(p-1)])+1):sum(segments[1:p])] <- p}

    # characterize segments
    unique_segments = unique(date_segments)
    for (s in 1:length(unique_segments)) {

      si = ri[which(dates[ri] %in% region_dates[which(date_segments == s)])]

      # estimate running time difference
      day_freq = ddply(data.frame(date=dates[si], time=z[si]),
                       .(date), summarise,
                       start=min(time), end=max(time),
                       day_cover=day_overlap(time)$day.cover,
                       day_overlap=day_overlap(time)$day.overlap,
                       data_frequency=as.numeric(mean(diff(time))/60))

      # extract segment statistics
      tmp[[(length(tmp)+1)]] = data.frame(
        region_id=r, segment_id=s,
        start.date=min(day_freq$start),
        end.date=max(day_freq$end),
        day_cover_mean=mean(day_freq$day_cover),
        day_cover_sd=sd(day_freq$day_cover),
        day_overlap_mean=mean(day_freq$day_overlap),
        day_overlap_sd=sd(day_freq$day_overlap),
        nr_days=length(unique(z[si])),
        nr_individuals=length(unique(y[si])),
        nr_records=length(si),
        data_frequency_mean=mean(day_freq$data_frequency),
        data_frequency_sd=mean(day_freq$data_frequency)
      )

      # build polygon for segment if sufficient data points exist
      if (length(si) > 2) {
        h = convHull(x[si,])
        h$region_id = r
        h$segment_id = s
        region_shp[[(length(region_shp)+1)]] = h
        rm(h)
      }

    }

  }

  tmp = do.call(rbind, tmp)
  region_shp = do.call(rbind, region_shp)

  #---------------------------------------------------------------------------#
  # 4. compile outputs and report
  #---------------------------------------------------------------------------#

  return(list(region.id=regions, region.polygons=region_shp, region.stats=tmp))



}
