#' @title proSat
#'
#' @description Interface to download and process satelite data.
#' @param t.var Target variable.
#' @param xy Object from which an extent can be derived.
#' @param o.time Object of class \emph{Date}.
#' @param d.path Output data path for downloaded data.
#' @param p.raster Logical. Should the output be re-projected?
#' @param p.res Target pixel resolution (if p.raster is TRUE).
#' @param user.cred Two element character vector containing username and password.
#' @import grDevices sp rgdal ncdf4
#' @importFrom XML htmlParse readHTMLTable
#' @importFrom httr GET write_disk
#' @importFrom RCurl getURL
#' @importFrom gdalUtils gdal_translate
#' @return One or multiple raster objects.
#' @details {Downloads and pre-processes pre-selected satellite datasets That provide ecologically 
#' meaningful variables. \emph{xy} is required to define a reference extent for which the data will 
#' be downloaded and croped. If \emph{p.raster} is TRUE, the function will also re-project the output 
#' to the projection of \emph{xy}. Note that \emph{xy} does not required reprojection a priori. The 
#' function will find the approapriate projection to use when deriving a reference extent. If \emph{p.raster} 
#' is TRUE \emph{p.res} is also required to determine the output pixel resolution. The outputs are presented 
#' in GTiff format and are croped and masked by default. Note that the download files might not correspond to 
#' the dates supplied through \emph{o.time}. The function will consider the temporal resolution of the target 
#' dataset and compare the possible and the requested dates. If a requested date is not found, the 
#' function will instead provide the closest image time. As a consequence, the number of downloaded 
#' files might be lesser than the number of dates. The function will inform the user on this by providing 
#' a list object which contains the downloaded dates (\emph{$date}) as well as the path for the correspondent 
#' image (\emph{$path}).}
#' @note {To consult the list of provided variables, run the function without specifying any input.
#' This will provide information on the variables, variable codes (used in \emph{t.var}), temporal and 
#' spatial resolution and the sensor of origin. Some variables might require login credentials. Check 
#' the table to know which credentials to use and assign them through \emph{user.cred}.}
#' @seealso \code{\link{getEnv}} \code{\link{imgInt}} \code{\link{dataQuery}}
#' @examples {
#'  
#'  # return list of variables
#'  modis.var <- proSat(sensor="MODIS")
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

proSat <- function(t.var=NULL, xy=NULL, o.time=NULL, d.path=NULL, p.raster=FALSE, p.res=NULL, user.cred=NULL) {
  
  #-----------------------------------------------------------------------------------------------------------------#
  # 1. check input variables  
  #-----------------------------------------------------------------------------------------------------------------#
  
  # read variable list
  var.ls <- system.file('extdata', 'sat-variables.csv', package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
  
  # if no variable is selected return list of variables
  if (is.null(t.var)) {return(var.ls)}
  if (length(t.var)>1) {stop('"var" has more than 1 element')}
  loc <- which(var.ls$code==t.var)
  if (length(loc)==0) {stop('"t.var" is not a recognized variable')}
  sensor <- var.ls$sensor[loc]
  t.res <- var.ls$temporal.resolution..days.[loc]
  
  # read var list for target sensor
  var.ls <- system.file('extdata', paste0(sensor, '_sat-variables.csv'), package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
  
  # check if credentials are needed
  if (!is.na(var.ls$login[loc]) & is.null(user.cred)) {stop(paste0('target variable requires ', var.ls$Login[loc], ' credentials'))}
  
  # check output data path
  if (!dir.exists(d.path)) {stop('"d.path" not found in file system')}
  
  # check time
  if (is.null(o.time)) {stop('please provide "o.time"')}
  if (class(o.time)[1]!='Date') {stop('"o.time" is not of class "Date"')}
  
  # check reference file
  if (is.null(xy)) {stop('please provide "xy"')}
  ref <- projectExtent(xy, crs(var.ls$crs[loc])) # croping extent
  if (p.raster) {
    if (is.null(p.res)) {stop('re-projection requested. Please specify "p.res"')}
    if (length(p.res) > 2) {stop('lenght of "p.res" is greater than 2')}
    if (!is.numeric(p.res)) {stop('"p.res is not numeric')}
    rr <- raster(extent(xy), res=p.res, crs=crs(xy))} # reference raster
  
#-----------------------------------------------------------------------------------------------------------------#
# 2. determine year/day  
#-----------------------------------------------------------------------------------------------------------------#
  
  ud <- unique(o.time) # unique dates
  yrs <- sapply(as.character(ud), function(x) {strsplit(x, '-')[[1]][1]})
  doa <- (ud-as.Date(paste0(yrs, '-01-01')))+1
  
  # check which dates can be downloaded
  potential.doa <- seq(1, 361, t.res)
  
  
  # update temporal information
  ud <- do.call('c', sapply(tmp, function(x) {x$date}))
  ind <- !duplicated(ud)
  ud <- ud[ind]
  doa <- unlist(sapply(tmp, function(x) {x$doa}))[ind]
  yrs <- unlist(sapply(tmp, function(x) {x$year}))[ind]
  
  rm(tmp)
  
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
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], write_disk(ofile[f], overwrite=T))
        
        # dread raster (1)
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp1, sd_index=1)
        r1 <- raster(tmp1)
        
        # read raster (2)
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp2, sd_index=2)
        r2 <- raster(tmp2)
        
        # read raster (3)
        tmp3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp3, sd_index=3)
        r3 <- raster(tmp3)
        
        # add raster data to raster list
        ndvi <- (r2-r1) / (r2+r1)
        r3 <- ((r3 %% b[1])>=a[1])^2 + ((r3 %% b[2])>=a[2])^2 + ((r3 %% b[3])>=a[3])^2
        ndvi[r3>0 | ndvi < -1 | ndvi > 1] <- NA
        rfiles[[f]] <- ndvi
        
        rm(r1, r2, r3, ndvi)
        file.remove(ofile[f], tmp1, tmp2, tmp3)
        
      }
      
      # mosaic raster list
      return(do.call(mosaic, rfiles))
      
    } else {
      
      # download image
      GET(ifile, write_disk(ofile, overwrite=T))
      
      # dread raster (1)
      tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp1, sd_index=1)
      r1 <- raster(tmp1)
      
      # read raster (2)
      tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp2, sd_index=2)
      r2 <- raster(tmp2)
      
      # read raster (3)
      tmp3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp3, sd_index=3)
      r3 <- raster(tmp3)
      
      # derive/return ndvi
      ndvi <- (r2-r1) / (r2+r1)
      r3 <- ((r3 %% b[1])>=a[1])^2 + ((r3 %% b[2])>=a[2])^2 + ((r3 %% b[3])>=a[3])^2
      ndvi[r3>0 | ndvi < -1 | ndvi > 1] <- NA
      return(ndvi)
      
      rm(r1, r2, r3)
      file.remove(ofile, tmp1, tmp2, tmp3)
      
    }
  }
  
  # land surface temperature
  dwnLST <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], write_disk(ofile[f], overwrite=T))
        
        # read raster(1)
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp1, sd_index=1)
        r1 <- raster(tmp1)
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp2, sd_index=2)
        qc <- raster(tmp2)
        qc <- ((qc %% b[1])>=a[1])^2 + ((qc %% b[2])>=a[2])^2
        r1[qc>0] <- NA
        
        # read raster (2)
        tmp3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp3, sd_index=5)
        r2 <- raster(tmp3)
        tmp4 <- tempfile(pattern="tmp4", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp4, sd_index=6)
        qc <- raster(tmp4)
        qc <- ((qc %% b[1])>=a[1])^2 + ((qc %% b[2])>=a[2])^2
        r2[qc>0] <- NA
        
        # update raster list
        rfiles[[f]] <- stack(r1, r2)
        rm(r1, r2,qc)
        file.remove(ofile[f], tmp1, tmp2, tmp3, tmp4)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      GET(ifile, write_disk(ofile, overwrite=T))
      
      # read raster(1)
      tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp1, sd_index=1)
      r1 <- raster(tmp1)
      tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp2, sd_index=2)
      qc <- raster(tmp2)
      qc <- ((qc %% b[1])>=a[1])^2 + ((qc %% b[2])>=a[2])^2
      r1[qc>0] <- NA
      
      # read raster (2)
      tmp3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp3, sd_index=5)
      r2 <- raster(tmp3)
      tmp4 <- tempfile(pattern="tmp4", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp4, sd_index=6)
      qc <- raster(tmp4)
      qc <- ((qc %% b[1])>=a[1])^2 + ((qc %% b[2])>=a[2])^2
      r2[qc>0] <- NA
      
      # update raster list
      r0 <- stack(r1, r2)
      rm(r1, r2, qc)
      file.remove(ofile, tmp1, tmp2, tmp3, tmp4)
      
    }
    
    # name and return files
    names(r0) <- c("day", "night")
    return(r0)
    
  }
  
  # leaf area index
  dwnLAI <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], write_disk(ofile[f], overwrite=T))
        
        # read raster(1)
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp1, sd_index=2)
        r1 <- raster(tmp1)
        
        # read raster(2)
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp2, sd_index=6)
        r2 <- raster(tmp2)
        
        # quality information
        tm3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp3, sd_index=3)
        qc <- raster(tmp3)
        qc <- ((qc %% b[4])>=a[4])^2 + ((qc %% b[5])>=a[5])^2
        r1[qc>0] <- NA
        r2[qc>0] <- NA
        
        # update raster list
        rfiles[[f]] <- stack(r1, r2)
        rm(r1, r2, qc)
        file.remove(ofile[f], tmp1, tmp2, tmp3)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      GET(ifile, write_disk(ofile, overwrite=T))
      
      # read raster(1)
      tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp1, sd_index=2)
      r1 <- raster(tmp1)
      
      # read raster(2)
      tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp2, sd_index=6)
      r2 <- raster(tmp2)
      
      # quality information
      tm3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp3, sd_index=3)
      qc <- raster(tmp3)
      qc <- ((qc %% b[4])>=a[4])^2 + ((qc %% b[5])>=a[5])^2
      r1[qc>0] <- NA
      r2[qc>0] <- NA
      
      # update raster list
      r0 <- stack(r1, r2)
      rm(r1, r2, qc)
      file.remove(ofile, tmp1, tmp2, tmp3)
      
    }
    
    # name and return files
    names(r0) <- c("mean", "sd")
    return(r0)
    
  }
  
  # fraction of photosinthetically active radiation
  dwnFPAR <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], write_disk(ofile[f], overwrite=T))
        
        # read raster(1)
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp1, sd_index=1)
        r1 <- raster(tmp1)
        
        # read raster(2)
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp2, sd_index=5)
        r2 <- raster(tmp2)
        
        # quality information
        tmp3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp3, sd_index=3)
        qc <- raster(tmp3)
        qc <- ((qc %% b[4])>=a[4])^2 + ((qc %% b[5])>=a[5])^2
        r1[qc>0] <- NA
        r2[qc>0] <- NA
        
        # update raster list
        rfiles[[f]] <- stack(r1, r2)
        rm(r1, r2, qc)
        file.remove(ofile[f], tmp1, tmp2, tmp3)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      GET(ifile, write_disk(ofile, overwrite=T))
      
      # read raster(1)
      tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp1, sd_index=1)
      r1 <- raster(tmp1)
      
      # read raster(2)
      tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp2, sd_index=5)
      r2 <- raster(tmp2)
      
      # quality information
      tmp3 <- tempfile(pattern="tmp3", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp3, sd_index=3)
      qc <- raster(tmp3)
      qc <- ((qc %% b[4])>=a[4])^2 + ((qc %% b[5])>=a[5])^2
      r1[qc>0] <- NA
      r2[qc>0] <- NA
      
      # update raster list
      r0 <- stack(r1, r2)
      rm(r1, r2, qc)
      file.remove(ofile, tmp1, tmp2, tmp3)
      
    }
    
    # name and return files
    names(r0) <- c("mean", "sd")
    return(r0)
    
  }
  
  # fire mask
  dwnFIRE <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], write_disk(ofile[f], overwrite=T))
        
        # read raster
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=1)
        rfiles[[f]] <- raster(tmp) == 9
        
        file.remove(ofile[f], tmp)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      GET(ifile, write_disk(ofile, overwrite=T))
      
      # read raster
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp, sd_index=1)
      r0 <- raster(tmp) == 9
      
      file.remove(ofile, tmp)
      
    }
    
    # return file
    return(r0)
    
  }
  
  # snow cover extent
  dwnSNW <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], authenticate(user.cred[1], user.cred[2]), write_disk(ofile[f], overwrite=T))
        
        # read raster
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp1, sd_index=1)
        r0 <- raster(tmp1)
        
        # quality control
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp2, sd_index=2)
        r2 <- raster(tmp2)
        
        # update image
        r0[r0>100 | r2>2] <- NA
        
        file.remove(ofile, tmp1, tmp2)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      GET(ifile, authenticate(user.cred[1], user.cred[2]), write_disk(ofile, overwrite=T))
      
      # read raster
      tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp1, sd_index=1)
      r0 <- raster(tmp1)
      
      # quality control
      tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile, tmp2, sd_index=2)
      r2 <- raster(tmp2)
      
      # update image
      r0[r0>100 | r2>2] <- NA
      
      file.remove(ofile, tmp1, tmp2)
      
    }
    
    # return file
    return(r0)
    
  }
  
  # sea ice extent
  dwnICE <- function(ifile) {
    
    ofile <- tempfile(pattern=basename(ifile), tmpdir=tempdir(), fileext=".hdf")
    
    if (length(ifile)>1) {
      
      rfiles <- vector('list', length(ifile))
      
      for (f in 1:length(ofile)) {
        
        # download image
        GET(ifile[f], authenticate(user.cred[1], user.cred[2]), write_disk(ofile[f], overwrite=T))
        
        # read raster
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        gdal_translate(ofile[f], tmp, sd_index=1)
        rfiles[[f]] <- raster(tmp) == 200
        
        file.remove(ofile[f], tmp)
        
      }
      
      # mosaic raster list
      r0 <- do.call(mosaic, rfiles)
      
      rm(rfiles)
      
    } else {
      
      # download image
      GET(ifile, authenticate(user.cred[1], user.cred[2]), write_disk(ofile, overwrite=T))
      
      # read raster
      tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
      gdal_translate(ofile[f], tmp, sd_index=1)
      rfiles[[f]] <- raster(tmp) == 200
      
      file.remove(ofile[f], tmp)
      
    }
    
    # return file
    return(r0)
    
  }
  
#-----------------------------------------------------------------------------------------------------------------#
# 4. download/mosaic/crop/write raster data
#-----------------------------------------------------------------------------------------------------------------#
  
  # output file list
  f.ls <- vector("list", length(ud))
  
  if (sensor=="MODIS") {
    
    # format doa
    doa <- sprintf("%03d", doa)
    
    # read tile shapefile
    if (var.ls$type[loc]=="tile") {
      file <- system.file('extdata', paste0(sensor, '-tiles.shp'), package="rsMove")
      if (file=='') {
        unzip(system.file('extdata', paste0(sensor, '.zip'), package="rsMove"), exdir=system.file('extdata', package="rsMove"))
        file <- system.file('extdata', paste0(sensor, '-tiles.shp'), package="rsMove")}
      shp <- shapefile(file)
      tiles <- crop(shp, projectExtent(ref, crs(shp)))@data$tile}
    
    # download and mosaic each file
    for (d in 1:length(ud)) {
      
      # data to mosaic
      r.data <- vector('list', 2)
      
      # build file path (terra)
      if (var.ls$type[loc]=="tile") {
        if (t.var!="snow" & t.var!="ice") {
          server <- paste0(var.ls[loc,1], yrs[d], '/', doa[d], '/')
          tbl = readHTMLTable(xmlRoot(htmlParse(GET(url=server))), skip.rows=1)$V1
          file <- as.character(sapply(tiles, function(x) {paste0(server, tbl[grep(x, tbl)])}))}
        if (t.var=="snow" | t.var=="ice") {
          server <- paste0(var.ls[loc,1], paste0(strsplit(as.character(ud[d]), '-')[[1]], collapse='.'), '/')
          tbl <- as.character(readHTMLTable(xmlRoot(htmlParse(GET(url=server, authenticate(user.cred[1], user.cred[2])))))$V2)
          file <- as.character(sapply(tiles, function(x) {
            paste0(server, tbl[grep(paste0(x, ".*.hdf$"), tbl)])}))}}
      if (t.var=="chlorofile" | t.var=="sst") {
        server <- paste0(var.ls[loc,1], yrs[d], '/')
        tbl <- as.character(readHTMLTable(xmlRoot(htmlParse(GET(url=server))), skip.rows=1)$V1)
        file <- tbl[grep(paste0(yrs[d], doa[d]), tbl)]}
      
      # download data (aqua)
      if (t.var=="ndvi") {r.data[[1]] <- crop(dwnVI(file), ref)}
      if (t.var=="lst") {r.data[[1]] <- crop(dwnLST(file), ref)}
      if (t.var=="lai") {r.data[[1]] <- crop(dwnLAI(file), ref)}
      if (t.var=="fpar") {r.data[[1]] <- crop(dwnFPAR(file), ref)}
      if (t.var=="fire") {r.data[[1]] <- crop(dwnFIRE(file), ref)}
      if (t.var=="snow") {r.data[[1]] <- crop(dwnSNW(file), ref)}
      if (t.var=="ice") {r.data[[1]] <- crop(dwnICE(file), ref)}
      if (t.var=="chlorofile") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".nc")
        GET(paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/", file), write_disk(tmp, overwrite=T))
        r.data[[1]] <- crop(brick(tmp, var="chlor_a")[[1]], ref)}
      if (t.var=="sst") {
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".nc")
        GET(paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/", file[grep("_SST", file)]), write_disk(tmp1, overwrite=T))
        r1 <- crop(brick(tmp1, var="sst")[[1]], ref)
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".nc")
        GET(paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/", file[grep("_NSST", file)]), write_disk(tmp2, overwrite=T))
        r2 <- crop(brick(tmp2, var="sst")[[1]], ref)
        r0 <- stack(r1, r2)
        names(r0) <- c('day', 'night')
        r.data[[1]] <- r0
        rm(r1, r2, r0)}
      if (t.var=="cw") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        GET(paste0(var.ls[loc,1], "MODAL2_D_CLD_WP_", as.character(ud[d]), ".FLOAT.TIFF"), write_disk(tmp, overwrite=T))
        r0 <- crop(brick(tmp), ref)
        r0[r0>1000] <- NA
        r.data[[1]] <- r0
        rm(r0)}
      
      #-----------------------------------------------------------------------------------------------------------------#
      
      # build file path (aqua)
      if (var.ls$type[loc]=="tile") {
        if (t.var!="snow" & t.var!="ice") {
          server <- paste0(var.ls[loc,2], yrs[d], '/', doa[d], '/')
          tbl = readHTMLTable(xmlRoot(htmlParse(GET(url=server))), skip.rows=1)$V1
          file <- as.character(sapply(tiles, function(x) {paste0(server, tbl[grep(x, tbl)])}))}
        if (t.var=="snow" | t.var=="ice") {
          server <- paste0(var.ls[loc,2], paste0(strsplit(as.character(ud[d]), '-')[[1]], collapse='.'), '/')
          tbl <- as.character(readHTMLTable(xmlRoot(htmlParse(GET(url=server, authenticate(user.cred[1], user.cred[2])))))$V2)
          file <- as.character(sapply(tiles, function(x) {
            paste0(server, tbl[grep(paste0(x, ".*.hdf$"), tbl)])}))}}
      if (t.var=="chlorofile" | t.var=="sst") {
        server <- paste0(var.ls[loc,2], yrs[d], '/')
        tbl <- as.character(readHTMLTable(xmlRoot(htmlParse(GET(url=server))), skip.rows=1)$V1)
        file <- tbl[grep(paste0(yrs[d], doa[d]), tbl)]}
      
      # download data (aqua)
      if (t.var=="ndvi") {r.data[[2]] <- crop(dwnVI(file), ref)}
      if (t.var=="lst") {r.data[[2]] <- crop(dwnLST(file), ref)}
      if (t.var=="lai") {r.data[[2]] <- crop(dwnLAI(file), ref)}
      if (t.var=="fpar") {r.data[[2]] <- crop(dwnFPAR(file), ref)}
      if (t.var=="fire") {r.data[[2]] <- crop(dwnFIRE(file), ref)}
      if (t.var=="snow") {r.data[[2]] <- crop(dwnSNW(file), ref)}
      if (t.var=="ice") {r.data[[2]] <- crop(dwnICE(file), ref)}
      if (t.var=="chlorofile") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".nc")
        GET(paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/", file), write_disk(tmp, overwrite=T))
        r.data[[2]] <- crop(brick(tmp, var="chlor_a")[[1]], ref)}
      if (t.var=="sst") {
        tmp1 <- tempfile(pattern="tmp1", tmpdir=tempdir(), fileext=".nc")
        GET(paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/", file[grep("_SST", file)]), write_disk(tmp1, overwrite=T))
        r1 <- crop(brick(tmp1, var="sst")[[1]], ref)
        tmp2 <- tempfile(pattern="tmp2", tmpdir=tempdir(), fileext=".nc")
        GET(paste0("https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/", file[grep("_NSST", file)]), write_disk(tmp2, overwrite=T))
        r2 <- crop(brick(tmp2, var="sst")[[1]], ref)
        r0 <- stack(r1, r2)
        names(r0) <- c('day', 'night')
        r.data[[2]] <- r0
        rm(r1, r2, r0)}
      if (t.var=="cw") {
        tmp <- tempfile(pattern="tmp", tmpdir=tempdir(), fileext=".tif")
        GET(paste0(var.ls[loc,2], "MYDAL2_D_CLD_WP_", as.character(ud[d]), ".FLOAT.TIFF"), write_disk(tmp, overwrite=T))
        r0 <- crop(brick(tmp), ref)
        r0[r0>1000] <- NA
        r.data[[2]] <- r0
        rm(r0)}
      
      # average, crop and write raster
      nb <- nlayers(r.data[[1]]) # number of bands
      bn <- names(r.data[[1]]) # band names
      r.data <- lapply(1:nb, function(x) {calc(do.call(stack, lapply(1:2, function(y) {(r.data[[y]])[[x]]})), mean, na.rm=T)})
      if (nb > 0) {r.data <- do.call(stack, r.data)} else {r.data <- r.data[[1]]}
      names(r.data) <- bn
      if (p.raster) {r.data <- projectRaster(r.data, rr)}
      ofile <- paste0(file.path(d.path, fsep=.Platform$file.sep), .Platform$file.sep, as.character(ud[d]), '_', t.var, '.tif')
      writeRaster(r.data, ofile, driver="GTiff", datatype=dataType(r.data), overwrite=T)
      
      # update file list
      f.ls[[d]] <- ofile
      
      rm(r.data)
      
    }
  }
  
  # return information on downloaded files
  return(list(path=unlist(f.ls), date=ud))
  
}