#' @title getEnv
#'
#' @description Interface to download data from earthenv.org.
#' @param dpath Output data path for downloaded data.
#' @param d.source Data source. One of "EarthEnv", "GFC" or "GSW".
#' @param var Target variables.
#' @param ref Object from which an extent can be derived.
#' @import grDevices rgeos sp rgdal
#' @importFrom utils download.file
#' @return One or multiple raster objects.
#' @details {Downloads data from earthenv.org. To check which variables can be downloaded, 
#' run the function without specifying \emph{var} and specifying \emph{d.source}. This will 
#' return a data frame list the existing variables for a given data source. Here, the user 
#' can refer to the column \emph{"code"} to retrieve the keywords that can be passed to 
#' the function through \emph{var}. The keywords recognized by \emph{d.source} are:
#' \itemize{
#' \item{\emph{"EarthEnv"} - EarthEnv project.}
#' \item{\emph{"GFC"} - Maryland University Global Forest Change.}
#' \item{\emph{"GSW"} - JRC Global Surface Water.}}
#' If \emph{var} contains \emph{"DEM90"} from {"EarthEv"} or any variable from \emph{"GFC"} 
#' or \emph{"GSW"}, \emph{ref} is required. This will be used to determine which tiles should 
#' be downloaded. For details on the specifications of the provided datasets please consult 
#' the referenced websites.}
#' @references {\url{http://www.earthenv.org/} 
#' \url{https://earthenginepartners.appspot.com/science-2013-global-forest/} 
#' \url{https://global-surface-water.appspot.com//}}
#' @seealso \code{\link{sMoveRes}}
#' @examples {
#'  
#'  # return list of variables
#'  ee.var <- getEnv(source="earthEnv")
#'  gfc.var <- getEnv(source="GFC")
#'  gsw.var <- getEnv(source="GSW")
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

getEnv <- function(dpath=NULL, d.source=NULL, var=NULL, ref=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 1. load variable list  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # check data source keyword
  if (is.null(d.source)) {stop('please assign a valid keyword to "d.source"')}
  if (!d.source%in%c('EartEnv', 'GFC', 'GSW')) {stop('"d.source" is not a valid keyword')}
  
  # read variable list
  var.ls <- system.file('extdata', paste0(d.source, '_variables.csv'), package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
  
  # return variable list if none is specified
  if (is.null(var)) {return(var.ls[,2:3])}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 2. check parameters
#-------------------------------------------------------------------------------------------------------------------------------#
  
  if (is.null(dpath)) {stop('"dpath" is missing')}
  if (dir.exists(dpath)) {stop('could not find "dpath" in the file system')}
  if (min(var%in%var.ls$code)==0) {stop('one or more elements in "var" not found')}
  if ('DEM90'%in%var) {if (is.null(ref)) {
    stop('"DEM90" was set in "var". "ref" is required')}
    ext <- extent(ref)}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 3. download each variable
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # target variables
  ind <- which(var.ls$code %in% var)
  
  # loop through each target variable
  for (i in 1:length(ind)) {
    
    # simple download
    if (var[ind[i]] != 'DEM90' & !d.source%in%c('GFC', 'GSW')) {
      ofile <- file.path(dpath, basename(var.ls$link[ind[i]]))
      download.file(var.ls$link[ind[i]], ofile, mode="wb")}
    
    # DEM download
    if (var.ls$code[ind[i]] == 'DEM90' | d.source%in%c('GFC', 'GSW')) {
      
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
        ofile <- file.path(dpath, basename(ifile)) # output
        download.file(ifile, ofile, mode="wb")} # download file
      
    }
  }
}