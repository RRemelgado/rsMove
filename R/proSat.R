#' @title proSat
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
#' \item{\emph{"HSM"} - JRC Human Settlement map.}}
#' If \emph{var} contains \emph{"DEM90"} from {"EarthEv"} or any variable from \emph{"GFC"} 
#' or \emph{"GSW"}, \emph{ref} is required. This will be used to determine which tiles should 
#' be downloaded. For details on the specifications of the provided datasets please consult 
#' the referenced websites.}
#' @references {\url{http://www.earthenv.org/} 
#' \url{https://earthenginepartners.appspot.com/science-2013-global-forest/} 
#' \url{https://global-surface-water.appspot.com/} \url{http://ghsl.jrc.ec.europa.eu/} 
#' \url{http://maps.elie.ucl.ac.be/CCI/viewer/}}
#' @seealso \code{\link{getEnv}}
#' @examples {
#'  
#'  # return list of variables
#'  modis.var <- getEnv(d.source="EarthEnv")
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

proSat <- function(sensor=sensor, var=NULL, dpath=NULL, ref=NULL) {
  
#-----------------------------------------------------------------------------------------------------------------#
# 1. check input variables  
#-----------------------------------------------------------------------------------------------------------------#
  
  # input sensor exists?
  if (is.null(sensor)) {stop('please assign a valid keyword to "sensor"')}
  if (!sensor%in%c('MODIS')) {stop('"sensor" does not contain a valid keyword')}
  
  # read variable list
  var.ls <- system.file('extdata', paste0(sensor, '_sat-variables.csv'), package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
    
  # if no variable is selected return list of variables
  if (is.null(var)) {return(var.ls[,4:6])}
  
  # check output data path
  if (!dir.exists(dpath)) {stop('"dpath" not found in file system')}
    
  # temporary data path
  if (is.null(rd.path)) {
    rv <- TRUE
    dpath <- tempdir()
  } else {rv <- FALSE}
    
#-----------------------------------------------------------------------------------------------------------------#
#   
#-----------------------------------------------------------------------------------------------------------------#
  
  if (sensor=="MODIS") {
    
    # 16-bit quality conversion variables
    a<-2^(0:15)
    b<-2*a
    
    # read tile shapefile
    file <- system.file('extdata', paste0(sensor, '-tiles.shp'), package="rsMove")
    if (file=='') {
      unzip(system.file('extdata', paste0(d.source, '.zip'), package="rsMove"), exdir=system.file('extdata', package="rsMove"))
      file <- system.file('extdata', paste0(d.source, '-tiles.shp'), package="rsMove")}
    shp <- shapefile(file)
    
    # lst processing
    if (var=="lst") {
      
      
      
      
    }
    
    # vi
    if (var=="ndvi") {
      
      r <- brick("")
      qc <- (((qinfo %% b[1])>=a[1])*2 + ((qinfo %% b[2])>=a[2])*2) == 0
      r
      
      
    }
    
    # ndwi
    if (var=="ndwi") {}
    
    # chlorophile
    if (var=="chlor") {}
    
    # sea surface temperature
    if (var=="sst") {
      
      img <-brick("C:/Users/rus14jh.UNI-WUERZBURG/Downloads/A2003001.L3m_DAY_NSST_sst_4km.nc", varname="qual_sst")
      
      
    }
    
    
    
  }
  
}