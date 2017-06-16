#' @title labelSample
#'
#' @description Region labeling of samples based on their spatial connectivity. First, the samples are converted to pixel coordinates and the connectivity between neighboring samples is evaluated. The resulting regions are filtered when required removing small sample groups. Then, the pixels occupied by samples are dilated and a new region labelling is performed. Regions which are within close proximity of each other are thus aggregated. Finnally, for each labeled pixel, the function provides mean sample coordinates and the corresponding pixel index and region ID.
#' @param x Vector of x coordinates.
#' @param y Vector of y coordinates.
#' @param rad Minimum radius within which sample regions are aggregated. Unit depends on the projection of the data.
#' @param npx Minimum number of pixels wihin a region. Before turning the final sample regions, the function removes regions where the sample count is below the specified value. Removing these samples allows to account for sampling anomalies. The default is 2.
#' @param var Optional. Object of class "matrix" or "data.frame". It provides additional numeric variables (e.g. time spent) that wiill be average on a per pixel basis and added to the output.
#' @param pxr Pixel resolution used to convert between geographic and pixel coordinates. Required if "layer" is not provided.
#' @param layer Raster object which can be provided as an alternative to "pxr".
#' @return Matrix containing the coordinates for selected pixels (x, y) accompanied by their matrix index (index), region ID (region). When provided, the output will contain a summary of "var".
#' @import raster, grDevices
#' @seealso \code{\link{sampleMove}}
#' @aliases
#' @examples \dontrun{
#'
#' }

#-------------------------------------------------------------------------------------------------------------------------------#

labelSample <- function(x=x, y=y, rad=rad, npx=2, var=NULL, pxr=NULL, layer=NULL) {
  
  # check input variables
  if (!exists(x)) {stop('error: missing x coordinates')}
  if (!exists(y)) {stop('error: missing y coordinates')}
  if (length(x)!=length(y)) {stop('error: "x" and "y" have different lenghts')}
  if (!exists(rad)) {stop('error: missing "rad"')}
  if (!is.null(var)) {
    if (!(class(var) %in% c('matrix', 'data.frame'))) {stop('error: "var" not a matrix/data.frame')}
    if (dim(var)==NULL) {if (length(var)!=length(x)) {stop('error: "var" has a different length from x/y')}}
    if (dim(var)[1]!=NULL) {if (nrow(var)!=length(x)) {stop('error: "var" has a different length from x/y')}}
  }
  
  #-------------------------------------------------------------------------------------------------------------------------------#
  
  # extract extent of study area
  if (is.null(layer)) {
    if (is.null(pxr)) {stop('error: no layer. pxr is required')}
    ext <- c(min(x), max(x), min(y), max(y))
    nc <- round((ext[2]-ext[1]) / pxr) + 1 # number of columns
    nr <- round((ext[4]-ext[3]) / pxr) + 1 # number of rows
    
  } else {
    if (!grDevices::is.raster(layer)) {stop('error: "layer" is not a valid raster object')}
    ext <- raster::extent(layer)
    pxr <- raster::res(layer)[1]
    nr <- dim(s)[1]
    nc <- dim(s)[2]
  }
  
  # derive pixel coordinates
  sp <- (round((ext[4]-y)/pxr)+1) + nr * round((x-ext[1])/pxr) # convert coordinates to pixel positions
  up <- unique(sp) # unique pixel positions
  if (length(up)==1) {stop('warning: only one pixel with data found. Processing aborted (is pxr correct?)')}
  
  #-------------------------------------------------------------------------------------------------------------------------------#
  
  # evaluate pixel connectivity and filter (if rad > 1)
  if (rad > 1) {
    regions <- matrix(0, nc, nr)
    for (r in 1:length(up)) {
      xp <- up %% nr
      yp <- up / nr
      if (xp > 1) {sc<-xp-1} else {sc<-xp}
      if (xp < nc) {ec<-xp+1} else {ec<-xp}
      if (yp > 1) {sr<-yp-1} else {sr<-yp}
      if (yp < nr) {er<-yr+1} else {er<-yp}
      if (max(regions[sc:ec,sr:er])>0) {
        mv <- min(regions[sc:ec,sr:er])
        uv <- unique((regions[sc:ec,sr:er])[which(regions[sc:ec,sr:er]>0)])
        for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
      } else {regions[xpos,ypos] <- max(regions)+1}
    }
    
    # estimate per region pixel count
    uv <- unique(regions[which(regions>0)])
    count <- sapply(uv, function(x) {length(which(regions==x))})
    
    # remove samples related to regions with a pixel count bellow npx
    uv <- uv[which(count>=npx)]
    for (r in 1:length(uv)) {regions[which(regions==uv[r])]=0}
    up <- up[which(regions[up]>0)]
    
    rm(regions, count, uv)
    
  }
  
  #-------------------------------------------------------------------------------------------------------------------------------#
  
  # define 
  
  # dilate samples
  upd <- sapply(up, function(x) {(((x %% nr)-k):((x %% nr)+k)) + nr * (((x / nr)-k):((x / nr)+k))})
  upd <- unique(upd[which(upd >= 1 & upd <= (nr*nc))])
  
  # evaluate sample connectivity
  regions <- matrix(0, nc, nr)
  for (r in 1:length(upd)) {
    xp <- upd %% nr
    yp <- upd / nr
    if (xp > 1) {sc<-xp-1} else {sc<-xp}
    if (xp < nc) {ec<-xp+1} else {ec<-xp}
    if (yp > 1) {sr<-yp-1} else {sr<-yp}
    if (yp < nr) {er<-yr+1} else {er<-yp}
    if (max(regions[sc:ec,sr:er])>0) {
      mv <- min(regions[sc:ec,sr:er])
      uv <- unique((regions[sc:ec,sr:er])[which(regions[sc:ec,sr:er]>0)])
      for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
    } else {regions[xpos,ypos] <- max(regions)+1}
  }
  
  # retrieve original samples
  regions <- regions[up]
  uv <- unique(regions[which(regions > 0)])
  uRegions <- regions * 0
  for (r in 1:length(uv)) {uregions[which(regions==uv[r])] <- r}
  
  rm(regions)
  
  #-------------------------------------------------------------------------------------------------------------------------------#
  
  # summarize input variables
  xr <- vector('numeric', length(up))
  yr <- vector('numeric', length(up))
  cn <- c('x', 'y', 'index', 'region')
  if (!is.null(var)) {
    ov <- matrix(0, length(up), ncol(var))
    cn <- c(cn, colnames(var))
  }
  for (r in 1:length(up)) {
    ind <- which(sp==x)
    xr[r] <- mean(x[ind])
    yr[r] <- mean(y[ind])
    if (!is.null(var)) {ov[r,] <- apply(var[ind,], 1, mean, na.rm=T)}
  }
  
  rm(x, y, sp, var)
  
  #-------------------------------------------------------------------------------------------------------------------------------#
  
  # build/return output
  if (!is.null(var)) {out <- cbind(xr, yr, up, region, ov)} else {
    out <- cbind(xr, yr, up, region)
  }
  rm(xr, yr, up, region, ov)
  colnames(out) <- cn
  return(out)
  
}