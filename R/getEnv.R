#' @title getEnv
#'
#' @description Interface to download ecologically relevant data.
#' @param d.path Output data path for downloaded data.
#' @param d.source Data source. One of "EarthEnv", "GFC", "GSW", "CCI" or "HSM".
#' @param t.var Target variable.
#' @param ref Projected spatial object from which an extent can be derived.
#' @param p.raster Logical. Should the output be reprojected?
#' @param p.res Pixel resolution (used if p.raster is TRUE).
#' @param pad Mi,ner pf cells used to pad the output raster when resampling.
#' @import grDevices sp rgdal raster
#' @importFrom utils download.file read.csv unzip untar
#' @return One or multiple raster objects.
#' @details {Downloads data from earthenv.org. To check which variables can be downloaded,
#' run the function without specifying \emph{t.var} and specifying \emph{d.source}. This will
#' return a data frame listing the existing variables for a given data source. Here, the user
#' can refer to the column \emph{"code"} to retrieve the keywords that can be passed to
#' the function through \emph{t.var}. The keywords recognized by \emph{d.source} are:
#' \itemize{
#' \item{\emph{"EarthEnv"} - EarthEnv project.}
#' \item{\emph{"GFC"} - Maryland University Global Forest Change.}
#' \item{\emph{"GSW"} - JRC Global Surface Water.}
#' \item{\emph{"CCI"} - ESA CCI Global land cover.}
#' \item{\emph{"HSM"} - JRC Human Settlement map.}
#' \item{\emph{"NEO"} - NASA Earth Observations.}}
#' If \emph{t.var} contains \emph{"DEM90"} from {"EarthEv"} or any variable from \emph{"GFC"}
#' or \emph{"GSW"}, The function will require a reference spatial object (\emph{ref}) to
#' determine which tiles to download. The user can also set \emph{p.raster} to TRUE and define
#' a reference spatial resolution (\emph{s.res}) to re-project the data. The use of \emph{ref})
#' with the remaining variables is optional. If prompted, These files can also be cropped and reprojected.}
#' @references {\url{http://www.earthenv.org/}
#' \url{https://earthenginepartners.appspot.com/science-2013-global-forest}
#' \url{https://global-surface-water.appspot.com/} \url{http://ghsl.jrc.ec.europa.eu/}
#' \url{http://maps.elie.ucl.ac.be/CCI/viewer/} \url{https://neo.sci.gsfc.nasa.gov/}}
#' @seealso \code{\link{dataQuery}} \code{\link{proSat}}
#' @examples \dontrun{
#'
#'  # return list of variables
#'  ee.var <- getEnv(d.source="EarthEnv")
#'  gfc.var <- getEnv(d.source="GFC")
#'  gsw.var <- getEnv(d.source="GSW")
#'  cci.var <- getEnv(d.source="CCI")
#'  hsm.var <- getEnv(d.source="HSM")
#'  neo.var <- getEnv(d.source="NEO")
#'
#'  # download bathimetry data
#'  getEnv(d.source="NEO", t.var="bat", d.path=".")
#'  getEnv(d.source="CCI", t.var="2015", d.path=".")
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

getEnv <- function(d.path=NULL, d.source=NULL, t.var=NULL, ref=NULL, p.raster=FALSE, p.res=NULL, pad=10) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. load variable list
#-------------------------------------------------------------------------------------------------------------------------------#

  # check data source keyword
  if (is.null(d.source)) {stop('please assign a valid keyword to "d.source"')}
  if (!d.source%in%c('EarthEnv', 'GFC', 'GSW', 'CCI', 'HSM', 'NEO')) {stop('"d.source" is not a valid keyword')}

  # read variable list
  var.ls <- system.file('extdata', paste0(d.source, '_variables.csv'), package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=FALSE)

  # return variable list if none is specified
  if (is.null(t.var)) {return(var.ls[,2:3])}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. check parameters
#-------------------------------------------------------------------------------------------------------------------------------#

  # target variables
  ind <- which(var.ls$code == t.var)
  if (length(ind)==0) {stop('"t.var" is not a valid variable')}

  # check reference fil
  if (!is.null(ref)) {
    ext <- extent(projectExtent(ref, crs(var.ls$crs[ind])))} #  project extent

  # check output path
  if (is.null(d.path)) {stop('"d.path" is missing')}
  if (!dir.exists(d.path)) {stop('could not find "d.path" in the file system')}
  d.path <- file.path(d.path)

  # check re-projection parameters
  if (p.raster) {
    if (is.null(p.res)) {stop('"p.raster" set to TRUE. Please define "p.res"')}
    if (!is.numeric(p.res)) {stop('"p.res" is not numeric')}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. download each variable
#-------------------------------------------------------------------------------------------------------------------------------#

  # simple download
  if (var.ls$code[ind] != 'DEM90' & !d.source%in%c('GFC', 'GSW')) {
    ofile <- file.path(d.path, basename(var.ls$link[ind]))
    download.file(var.ls$link[ind], ofile, quiet=TRUE, mode="wb")

    if (d.source=="HSM") {
      files = unzip(ofile, list=TRUE)
      files <- files$Name[grep("tif", files$Name)]
      unzip(ofile, files=files, exdir=d.path, junkpaths=TRUE, overwrite=TRUE)
      file.remove(ofile)
      ofile <- paste0(d.path, .Platform$file.sep, basename(files[grep("tif$", files)]))}

    # crop/re-rpoeject raster
    if (!is.null(ref)) {
      r.data <- raster(ofile)
      if (p.raster) {
        pxr <- res(r.data)[1]
        ext <- extend(ext, pxr*pad)
        r.data <- crop(r.data, ext)
        r.data <- crop(projectRaster(r.data, crs=crs(ref), res=p.res), extent(ref))}
      if (!p.raster) {r.data <- crop(r.data, extent(ref))}
      writeRaster(r.data, ofile, overwrite=TRUE)}
    }

  if (var.ls$code[ind] == 'DEM90' | d.source%in%c('GFC', 'GSW')) {

    # check if res exist
    if (is.null(ref)) {stop('Requested a tiled variable. Please provide "ref"')}

    # read in shapefile with tiles
    file <- system.file('extdata', paste0(d.source, '-tiles.shp'), package="rsMove")
    if (file=='') {
      unzip(system.file('extdata', paste0(d.source, '.zip'), package="rsMove"), exdir=system.file('extdata', package="rsMove"))
      file <- system.file('extdata', paste0(d.source, '-tiles.shp'), package="rsMove")}
    shp <- shapefile(file)

    # determine which tiles are required
    tiles <- as.character(crop(shp, ext)@data$tile)

    # retrieve data
    ofile <- vector('character', length(tiles))
    for (x in 1:length((tiles))) {
      if (d.source=='EarthEnv') {ifile <- paste0(var.ls$link[ind], tiles[x], '.tar.gz')}
      if (d.source=='GFC') {ifile <- paste0(var.ls$link[ind], '_', var.ls$code[ind], '_', tiles[x], '.tif')}
      if (d.source=='GSW') {ifile <- paste0(var.ls$link[ind], var.ls$code[ind], '_', tiles[x], '.tif')}
      ofile[x] <- tempfile(fileext=".tar.gz") # output
      download.file(ifile, ofile[x], quiet=TRUE, mode="wb")} # download file

    # mosaic data only (if GFC or GSW)
    if (d.source=='GFC' | d.source=='GSW') {
      fls <- lapply(ofile, function(x) {raster(x)})
      fls$fun <- mean
      if (length(tiles) > 1 ) {r.data <- do.call(mosaic, fls)}
      if (length(tiles) == 1) {r.data <- fls[[1]]}}

    # untar before mosaic (if EarthEnv)
    if (d.source=="EarthEnv") {
      for (x in 1:length(tiles)) {
        #untar(ofile[x], exdir=d.path, tar="internal")
        untar(ofile[x], exdir=d.path, tar="internal")
        file.remove(ofile)}
      files <- list.files(d.path, ".bil$", full.names=TRUE)
      fls <- lapply(files, function(x) {raster(files)})
      fls$fun <- mean
      if (length(tiles) > 1 ) {r.data <- do.call(mosaic, fls)}
      if (length(tiles) == 1) {r.data <- fls[[1]]}
      files <- list.files(d.path, "EarthEnv", full.names=TRUE)}

      # export raster
    if (!is.null(ref)) {
      if (p.raster) {
        pxr <- res(r.data)[1]
        ext <- extend(ext, pxr*pad)
        r.data <- crop(r.data, ext)
        r.data <- crop(projectRaster(r.data, crs=crs(ref), res=p.res), extent(ref))}
      if (!p.raster) {r.data <- crop(r.data, extent(ref))}}

    # write raster
    writeRaster(r.data, file.path(d.path, paste0(var.ls$code[ind], '_', paste0(tiles, collapse="-"), '.tif')), overwrite=TRUE)

  }

}
