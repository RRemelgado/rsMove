#' @title getEnv
#'
#' @description Interface to download data from earthenv.org.
#' @param dpath Output data path for downloaded data.
#' @param var Target variables.
#' @param ref Object from which an extent can be derived.
#' @import grDevices
#' @importFrom utils download.file
#' @return One or multiple raster objects.
#' @details {Downloads data from earthenv.org. To check which variables can be downloaded, 
#' run the function without specifying \emph{var}. This will return a table containing the 
#' codes for each existing variable. If \emph{var} contains \emph{"DEM90"}, \emph{ref} is 
#' required. This will be used to determine which tiles should be downloaded. For details 
#' on the specifications of these datasets please refer to the website of origin.}
#' @references \url{http://www.earthenv.org/}
#' @seealso \code{\link{sMoveRes}}
#' @examples {
#'  
#'  # return table with variables
#'  a.res <- getEnv()
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

getEnv <- function(dpath=NULL, var=NULL, ext=NULL) {
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 1. load variable list  
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # read variable list
  var.ls <- system.file('extdata', 'earth-env_variables.csv', package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
  
  # return variable list if none is specified
  if (is.null(var)) {return(var.ls[,2:3])}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 2. check parameters
#-------------------------------------------------------------------------------------------------------------------------------#
  
  if (is.null(dpath)) {stop('"dpath" is missing')}
  if (dir.exists(dpath)) {stop('could not find "dpath" in the file system')}
  if (min(var%in%var.ls$code)=0) {stop('one or more elements in "var" not found')}
  if ('DEM90'%in%var) {if (is.null(ext)) {
    stop('"DEM90" was set in "var". "ref" is required')}
    ext <- extent(ref)}
  
#-------------------------------------------------------------------------------------------------------------------------------#
# 3. download each variable
#-------------------------------------------------------------------------------------------------------------------------------#
  
  # target variables
  ind <- which(var.ls$code %in% var)
  
  # loop through each target variable
  for (i in 1:length(ind)) {
    
    if (var[ind[i]] != 'DEM90') {
      
      # specify output path
      ofile <- file.path(dpath, basename(var.ls$link[ind[i]]))
      
      # download file
      download.file(var.ls$link[ind[i]], ofile, mode="wb")
      
    } else {
      
      # read in shapefile with tiles
      file <- system.file('extdata', 'DEM90-tiles.shp', package="rsMove")
      if (file=='') {
        unzip(system.file('extdata', 'DEM90.zip', package="rsMove"), exdir=system.file('extdata', package="rsMove"))
        file <- system.file('extdata', 'DEM90-tiles.shp', package="rsMove")}
      shp <- shapefile(file)
      
      # project extent
      ext <- projectExtent(ref, crs(shp))
      
      # determine which tiles are required 
      tiles <- as.character(intersect(shp, spTransform(moveData, crs(shp)))@data$tile)
      
      # retrieve data
      for (x in 1:length((tiles))) {
        ifile <- paste0(var.ls$link[ind[i]], tiles[x], '.tar.gz') # origin
        ofile <- file.path(dpath, basename(var.ls$link[ind[i]])) # output
        download.file(ifile, ofile, mode="wb")} # download file
      
    }
  }
  
}