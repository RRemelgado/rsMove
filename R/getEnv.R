#' @title getEnv
#'
#' @description Interface to download data from earthenv.org.
#' @param dpath Output data path for downloaded data.
#' @param var Target variables.
#' @import grDevices
#' @importFrom utils download.file
#' @return One or multiple raster objects.
#' @details {Downloads data from earthenv.org. To check which variables can 
#' be downloaded, run the function without specifying \emph{var}. This will 
#' return a table containing the codes for each existing variable. For details 
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

getEnv <- function(dpath=NULL, var=NULL) {
  
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
      shp <- shapefile(system.file())
      
      # determine which tiles are required 
      tiles <- crop(shp, extent(xy))@data$tile
      
      # retrieve data
      for (x in 1:length((tiles))) {
        ifile <- paste0(var.ls$link[ind[i]], tiles[x], '.tar.gz') # origin
        ofile <- file.path(dpath, basename(var.ls$link[ind[i]])) # output
        download.file(ifile, ofile, mode="wb")} # download file
      
    }
  }
  
}