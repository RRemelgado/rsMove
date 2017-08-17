#' @title getEnv
#'
#' @description Interface to download ecologically relevant data.
#' @param d.path Output data path for downloaded data.
#' @param d.source Data source. One of "EarthEnv", "GFC", "GSW", "CCI" or "HSM".
#' @param var Target variables.
#' @param ref Object from which an extent can be derived.
#' @import grDevices sp rgdal
#' @importFrom utils download.file read.csv unzip
#' @return One or multiple raster objects.
#' @details {Downloads data from earthenv.org. To check which variables can be downloaded, 
#' run the function without specifying \emph{var} and specifying \emph{d.source}. This will 
#' return a data frame list the existing variables for a given data source. Here, the user 
#' can refer to the column \emph{"code"} to retrieve the keywords that can be passed to 
#' the function through \emph{var}. The keywords recognized by \emph{d.source} are:
#' \itemize{
#' \item{\emph{"EarthEnv"} - EarthEnv project.}
#' \item{\emph{"GFC"} - Maryland University Global Forest Change.}
#' \item{\emph{"GSW"} - JRC Global Surface Water.}
#' \item{\emph{"CCI"} - ESA CCI Global land cover.}
#' \item{\emph{"HSM"} - JRC Human Settlement map.}
#' \item{\emph{"NEO"} - NASA Earth Observations.}}
#' If \emph{var} contains \emph{"DEM90"} from {"EarthEv"} or any variable from \emph{"GFC"} 
#' or \emph{"GSW"}, \emph{ref} is required. This will be used to determine which tiles should 
#' be downloaded. For details on the specifications of the provided datasets please consult 
#' the referenced websites.}
#' @references {\url{http://www.earthenv.org/} 
#' \url{https://earthenginepartners.appspot.com/science-2013-global-forest/} 
#' \url{https://global-surface-water.appspot.com/} \url{http://ghsl.jrc.ec.europa.eu/}
#' \url{http://maps.elie.ucl.ac.be/CCI/viewer/} \url{https://neo.sci.gsfc.nasa.gov/}}
#' @seealso \code{\link{dataQuery}} \code{\link{spaceDir}}
#' @examples {
#'  
#'  # return list of variables
#'  ee.var <- getEnv(d.source="EarthEnv")
#'  gfc.var <- getEnv(d.source="GFC")
#'  gsw.var <- getEnv(d.source="GSW")
#'  cci.var <- getEnv(d.source="CCI")
#'  hsm.var <- getEnv(d.source="HSM")
#'  neo.var <- getEnv(d.source="NEO")
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

getEnv <- function(d.path=NULL, d.source=NULL, var=NULL, ref=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 1. load variable list  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # check data source keyword
  if (is.null(d.source)) {stop('please assign a valid keyword to "d.source"')}
  if (!d.source%in%c('EarthEnv', 'GFC', 'GSW', 'CCI', 'HSM', 'NEO')) {stop('"d.source" is not a valid keyword')}
  
  # read variable list
  var.ls <- system.file('extdata', paste0(d.source, '_variables.csv'), package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
  
  # return variable list if none is specified
  if (is.null(var)) {return(var.ls[,2:3])}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 2. check parameters
#-------------------------------------------------------------------------------------------------------------------------------#
  
  if (is.null(d.path)) {stop('"d.path" is missing')}
  if (!dir.exists(d.path)) {stop('could not find "d.path" in the file system')}
  if (min(var%in%var.ls$code)==0) {stop('one or more elements in "var" not found')}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 3. download each variable
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # target variables
  ind <- which(var.ls$code %in% var)
  
  # loop through each target variable
  for (i in 1:length(ind)) {
    
    # simple download
    if (var.ls$code[ind[i]] != 'DEM90' & !d.source%in%c('GFC', 'GSW')) {
      ofile <- file.path(d.path, basename(var.ls$link[ind[i]]))
      download.file(var.ls$link[ind[i]], ofile, quiet=TRUE, mode="wb")}
    
    # DEM download
    if (var.ls$code[ind[i]] == 'DEM90' | d.source%in%c('GFC', 'GSW')) {
      
      # control if ref was provided
      if (is.null(ref)) {stop('Selected a tiled product. Provide "ref"')}
      
      # read in shapefile with tiles
      file <- system.file('extdata', paste0(d.source, '-tiles.shp'), package="rsMove")
      if (file=='') {
        unzip(system.file('extdata', paste0(d.source, '.zip'), package="rsMove"), exdir=system.file('extdata', package="rsMove"))
        file <- system.file('extdata', paste0(d.source, '-tiles.shp'), package="rsMove")}
      shp <- shapefile(file)
      
      # project extent
      ext <- projectExtent(ref, crs(shp))
      
      # determine which tiles are required 
      tiles <- as.character(crop(shp, ext)@data$tile)
  
      # retrieve data
      for (x in 1:length((tiles))) {
        if (d.source=='EarthEnv') {ifile <- paste0(var.ls$link[ind[i]], tiles[x], '.tar.gz')}
        if (d.source=='GFC') {ifile <- paste0(var.ls$link[ind[i]], '_', var.ls$code[ind[i]], '_', tiles[x], '.tif')}
        if (d.source=='GSW') {ifile <- paste0(var.ls$link[ind[i]], var.ls$code[ind[i]], '_', tiles[x], '.tif')}
        ofile <- file.path(d.path, basename(ifile)) # output
        download.file(ifile, ofile, quiet=TRUE, mode="wb")} # download file
      
    }
  }
}