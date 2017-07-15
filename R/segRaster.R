#' @title segRaster
#'
#' @description Connencte-region based raster segmentation.
#' @param prob Object of class \emph{RasterLayer}.
#' @param pt Difference threshold. Default is 0.05.
#' @param mp Minimum value. Default is 0.5.
#' @import raster grDevices
#' @importFrom stats sd
#' @return A list object.
#' @details {The function segments an input layer using a connected 
#' component region labeling approach. For each pixel, the function 
#' estimates the difference between it and its imediate 8 neighbors. 
#' The pixels where the difference is below the defined threshold 
#' (\emph{ct}) are aggregated into a single region. The user can define 
#' a minimum pixel value using \emph{mp} which will limit the range of 
#' pixels under evaluation. The result contains a raster with unique 
#' values for each segment region region (\emph{$segment}) as well as a 
#' data frame (\emph{$stats}) with statistics for each region. The data 
#' frame report on the minimum (\emph{min}), maximum (\emph{max}), mean (\emph{mean}) 
#' and standard deviation (\emph{sd}) of the pixels contained in each region.}
#' @seealso \code{\link{moveModel}} \code{\link{modelApply}}
#' @examples {
#'  
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'  
#'  # load example probability image
#'  file <- system.file('extdata', 'konstanz_probabilities.tif', package="rsMove")
#'  probImg <- raster(file)
#'  
#'  # segment probabilities
#'  rs <- segRaster(probImg)
#'  
#' }
#' @export

#-----------------------------------------------------------------------------------#

segRaster <- function(prob, pt=0.1, mp=0.5) {
  
#-----------------------------------------------------------------------------------#
# 1. check input variables
#-----------------------------------------------------------------------------------#
  
  if (class(prob)[1]!='RasterLayer') {stop('"prob" is not a "RasterLayer"')}

#-----------------------------------------------------------------------------------#
# 2. segment regions
#-----------------------------------------------------------------------------------#
  
  nr <- dim(prob)[1]
  nc <- dim(prob)[2] 

  # identify usable pixels
  data <- as.matrix(prob)
  pos <- which(data > mp)
  
  # evaluate pixel connectivity
  regions <- matrix(0, nr, nc)
  for (r in 1:length(pos)) {
    rp <- ((pos[r]-1) %% nr)+1
    cp <- ((pos[r]-1) %/% nr)+1
    if (cp > 1) {sc<-cp-1} else {sc<-cp}
    if (cp < nc) {ec<-cp+1} else {ec<-cp}
    if (rp > 1) {sr<-rp-1} else {sr<-rp}
    if (rp < nr) {er<-rp+1} else {er<-rp}
    if (max(regions[sr:er,sc:ec])>0) {
      diff <- abs(data[sr:er,sc:ec]-data[rp,cp]) <= pt
      uv <- unique(c(regions[sr:er,sc:ec]*diff))
      uv <- uv[which(uv>0)]
      if (length(uv>0)) {
        mv <- min(uv)
        regions[rp,cp]<-min(uv)
        for (u in 1:length(uv)) {regions[which(regions==uv[u])] <- mv}
      } else {regions[rp,cp]<-max(regions)+1}
    } else {regions[rp,cp] <- max(regions)+1}
  }
  
#-----------------------------------------------------------------------------------#
# 3. derive segment statistics
#-----------------------------------------------------------------------------------#
    
  # update region id
  uv = sort(unique(regions[which(regions>0)]))
  uregions = regions
  nr <- length(uv)
  pmn <- vector('numeric', nr) # min
  pmx <- vector('numeric', nr) # max
  pav <- vector('numeric', nr) # mean
  psd <- vector('numeric', nr) # sd
  npx <- vector('numeric', nr) # count
  for (u in 1:nr) {
    pos <- which(regions==uv[u])
    uregions[pos] <- u
    pmn[u] <- min(data[pos], na.rm=T)
    pmx[u] <- max(data[pos], na.rm=T)
    pav[u] <- mean(data[pos], na.rm=T)
    psd[u] <- sd(data[pos], na.rm=T)
    npx[u] <- sum(!is.na(data[pos]))
  }
  
  rm(regions, data)
  
#-----------------------------------------------------------------------------------#
# 4. return output
#-----------------------------------------------------------------------------------#
  
  # convert data back to raster
  uregions = raster(uregions)
  extent(uregions) <- extent(prob)
  res(uregions) <- res(prob)
  crs(uregions) <- crs(prob)
  
  # build/return data frame
  df <- data.frame(segment=uv, min=pmn, max=pmx, mean=pav, sd=psd, cout=npx)
  
  # return matrix/df
  return(list(segment=uregions, stats=df))

}