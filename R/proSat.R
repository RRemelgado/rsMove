#' @title proSat
#'
#' @description Interface to download and process satelite data.
#' @param d.path Output data path for downloaded data.
#' @param d.source Data source. One of "MODIS", ...
#' @param var Target variables.
#' @param xy Object from which an extent can be derived.
#' @param p.raster Logical. Should the output be re-projected?
#' @param p.res Target pixel resolution (if p.raster is TRUE).
#' @import grDevices sp rgdal curl ncdf4
#' @importFrom utils download.file read.csv unzip
#' @return One or multiple raster objects.
#' @details {Downloads and pre-processes satellite data from selected sources. To check which 
#' variables can be downloaded, run the function without specifying \emph{var} and specifying 
#' \emph{sensor}. This will return a data frame list the existing variables for a given sensor. 
#' The user can selected from the following sensors:
#' \itemize{
#' \item{\emph{"MODIS"} - Moderate Resolution Spectroradiometer.}
#' \item{\emph{...} - }}
#' Independently of the selected sensor, \emph{xy} is required to define a reference extent for 
#' which the data will be downloaded and croped. If \emph{p.raster} is TRUE, the function will 
#' also re-project the output to the projection of \emph{xy}. Note that \emph{xy} does not 
#' required reprojection a priori. The function will find the approapriate projection to use 
#' when deriving a reference extent. If \emph{p.raster} is TRUE \emph{p.res} is also required 
#' to determine the output pixel resolution. The outputs are presented in GTiff format and are 
#' croped and masked by default. In the case of MODIS data, this function uses both 
#' TERRA and AQUA data to build the output image.}
#' @seealso \code{\link{getEnv}} \code{\link{imgInt}} \code{\link{dataQuery}}
#' @examples {
#'  
#'  # return list of variables
#'  modis.var <- proSat(sensor="MODIS")
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

proSat <- function(sensor=sensor, var=NULL, xy=NULL, o.time=NULL, dpath=NULL, p.raster=FALSE, p.res=NULL) {
  
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
  if (length(var)>1) {stop('"var" has more than 1 element')}
  loc <- which(var.ls$variable==var)
  
  # check output data path
  if (!dir.exists(dpath)) {stop('"dpath" not found in file system')}
  
  # check time
  if (is.null(o.time)) {stop('please provide "o.time"')}
  if (class(o.time)[1]!='Date') {stop('"o.time" is not of class "Date"')}
  
  # check reference file
  if (is.null(xy)) {stop('please provide "xy"')}
  ref <- projectExtent(xy, crs(var.ls$crs[loc])) # croping extent
  if (p.raster) {
    if (!is.null(p.res)) {stop('re-projection requested. Please specify "p.res"')}
    if (length(p.res > 2)) {stop('lenght of "p.res" is greater than 2')}
    if (!is.numeric(p.res)) {stop('"p.res is not numeric')}
    rr <- raster(extent(xy), res=p.res, crs=crs(xy))} # reference raster
  
#-----------------------------------------------------------------------------------------------------------------#
# 2. determine year/day  
#-----------------------------------------------------------------------------------------------------------------#
  
  ud <- unique(o.time) # unique dates
  yrs <- sapply(as.character(ud), function(x) {strsplit(x, '-')[[1]][1]})
  doa <- sprintf("%03d", ud - as.Date(paste0(yrs, '-01-01')))
  
#-----------------------------------------------------------------------------------------------------------------#
# 3. download functions
#-----------------------------------------------------------------------------------------------------------------#
  
  # 16-bit quality conversion variables
  a<-2^(0:15)
  b<-2*a
  
  # Normalized Difference Vegetation Index
  dwnVI <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile1)) {
        
        # download image
        download.file(ifile[f], ofile[f], quiet=TRUE, mode="wb")
        
        # dread raster (1)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=2)
        r1 <- raster(tmp)
        
        # read raster (2)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=3)
        r2 <- raster(tmp)
        
        # read raster (3)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=4)
        r3 <- raster(tmp)
        
        # remove temporary download file
        file.remove(ofile[f])
        
        # add raster data to raster list
        ndvi <- (r2-r1) / (r2+r1)
        r3 <- (((r3 %% b[1])>=a[1])*2 + ((r3 %% b[2])>=a[2])*2) == 0
        ndvi[r3!=1] <- NA
        rfiles[[f]] <- ndvi
        
        rm(r1, r2, r3, ndvi)
        
      }
      
      # mosaic raster list
      return(do.call(mosaic, rfiles))
      
    } else {
      
      # download image
      download.file(ifile, ofile, quiet=TRUE, mode="wb")
      
      # dread raster (1)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=2)
      r1 <- raster(tmp)
      
      # read raster (2)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=3)
      r2 <- raster(tmp)
      
      # read raster (3)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=4)
      r3 <- raster(tmp)
      
      # remove temporary download file
      file.remove(ofile)
      
      # add raster data to raster list
      ndvi <- (r2-r1) / (r2+r1)
      r3 <- (((r3 %% b[1])>=a[1])*2 + ((r3 %% b[2])>=a[2])*2) == 0
      ndvi[r3!=1] <- NA
      return(ndvi)
      
      rm(r1, r2, r3)
      
    }
  }
  
  # land surface temperature
  dwnLST <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile1)) {
        
        # download image
        download.file(ifile[f], ofile[f], quiet=TRUE, mode="wb")
        
        # dread raster(1)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=1)
        r1 <- raster(tmp)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=2)
        qc <- raster(tmp)
        qc <- (((qc %% b[1])>=a[1])*2 + ((qc %% b[2])>=a[2])*2) == 0
        r1[qc!=0] <- NA
        
        # dread raster (2)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=7)
        r2 <- raster(tmp)
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=8)
        qc <- raster(tmp)
        qc <- (((qc %% b[1])>=a[1])*2 + ((qc %% b[2])>=a[2])*2) == 0
        r2[qc!=0] <- NA
        
        # update raster list
        rfiles[[f]] <- stack(r1, r2)
        rm(r1, r2,qc)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      download.file(ifile, ofile, quiet=TRUE, mode="wb")
      
      # dread raster(1)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=1)
      r1 <- raster(tmp)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=2)
      qc <- raster(tmp)
      qc <- (((qc %% b[1])>=a[1])*2 + ((qc %% b[2])>=a[2])*2) == 0
      r1[qc!=0] <- NA
      
      # dread raster (2)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=7)
      r2 <- raster(tmp)
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=8)
      qc <- raster(tmp)
      qc <- (((qc %% b[1])>=a[1])*2 + ((qc %% b[2])>=a[2])*2) == 0
      r2[qc!=0] <- NA
      
      # update raster list
      r0 <- stack(r1, r2)
      
      rm(r1, r2,qc)
      
    }
    
    # name and return files
    names(r0) <- c("day", "night")
    return(r0)
    
  }
  
#-----------------------------------------------------------------------------------------------------------------#
# 4. download/mosaic/crop/write raster data
#-----------------------------------------------------------------------------------------------------------------#
    
  if (sensor=="MODIS") {
    
    # read tile shapefile
    file <- system.file('extdata', paste0(sensor, '-tiles.shp'), package="rsMove")
    if (file=='') {
      unzip(system.file('extdata', paste0(d.source, '.zip'), package="rsMove"), exdir=system.file('extdata', package="rsMove"))
      file <- system.file('extdata', paste0(d.source, '-tiles.shp'), package="rsMove")}
    shp <- shapefile(file)
    
    # check which tiles to use
    tiles <- crop(shp, projectExtent(ref, crs(shp)))@data$tile
    
    # download and mosaic each file
    for (d in 1:length(ud)) {
      
      # data to mosaic
      r.data <- vector('list', 2)
      
      # download aqua data
      if (var=="ndvi" |  var=="lst") {
        server <- paste0(var.ls[loc,1], yrs[d], '/', doa[d], '/')
        h = new_handle(dirlistonly=TRUE)
        con = curl(server, "r", h)
        tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
        close(con)
        file <- as.character(sapply(tiles, function(x) {paste0(server, tbl[grep(x, tbl$V1),1])}))
      } else {
        server <- paste0(var.ls[loc,1], yrs[d], '/')
        h = new_handle(dirlistonly=TRUE)
        con = curl(server, "r", h)
        tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
        close(con)
        file <- as.character(sapply(paste0(yrs[d], doa[d]), function(x) {paste0(server, tbl[grep(x, tbl$V1),1])}))}
      if (var=="ndvi") {r.data[[1]] <- crop(dwnVI(file), ref)}
      if (var=="lsp") {r.data[[1]] <- crop(dwnLSP(file), ref)}
      if (var=="chlorofile") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".nc")
        download.file(file, tmp, quiet=TRUE, mode="wb")
        r.data[[1]] <- crop(brick(tmp, var="chlor_a"), ref)}
      if (var=="sst") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".nc")
        download.file(file, tmp, quiet=TRUE, mode="wb")
        r.data[[1]] <- crop(brick(tmp, var="sst"), ref)}
      
      # download terra data
      if (var=="ndvi" |  var=="lst") {
        server <- paste0(var.ls[loc,2], yrs[d], '/', doa[d], '/')
        h = new_handle(dirlistonly=TRUE)
        con = curl(server, "r", h)
        tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
        close(con)
        file <- as.character(sapply(tiles, function(x) {paste0(server, tbl[grep(x, tbl$V1),1])}))
      } else {
        server <- paste0(var.ls[loc,2], yrs[d], '/')
        h = new_handle(dirlistonly=TRUE)
        con = curl(server, "r", h)
        tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
        close(con)
        file <- as.character(sapply(paste0(yrs[d], doa[d]), function(x) {paste0(server, tbl[grep(x, tbl$V1),1])}))}
      if (var=="ndvi") {r.data[[2]] <- crop(dwnVI(file), ref)}
      if (var=="lsp") {r.data[[2]] <- crop(dwnLSP(file), ref)}
      if (var=="chlorofile") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".nc")
        download.file(file, tmp, quiet=TRUE, mode="wb")
        r.data[[2]] <- crop(brick(tmp, var="chlor_a"), ref)}
      if (var=="sst") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".nc")
        download.file(file, tmp, quiet=TRUE, mode="wb")
        r.data[[2]] <- crop(brick(tmp, var="sst"), ref)}
      
      # mosaic, crop and write raster
      r.data <- calc(do.call(stack, r.data), mean, na.rm=T)
      if (p.raster) {projectRaster(r.data, rr)}
      ofile <- paste0(file.path(dpath, fsep=.Platform$file.sep), .Platform$file.sep, as.character(ud[d]), '_', var, '.tif')
      writeRaster(r.data, ofile, driver="GTiff", datatype=dataType(r.data), overwrite=T)
      
      rm(r.data)
      
    }
    
  }
    
}