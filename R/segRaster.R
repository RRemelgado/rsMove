#' @title segRaster
#'
#' @description {Connencte-region based raster segmentation that 
#' preserves spatial gradients and reports on region statistics.}
#' @param prob Object of class \emph{RasterLayer}.
#' @param break.point Difference threshold. Default is 0.05.
#' @param min.prob Minimum value. Default is 0.5.
#' @import raster rgdal
#' @importFrom stats sd
#' @return A list object.
#' @details {The function segments an input layer using a connected
#' component region labeling approach. For each pixel, the function
#' estimates the difference between it and its imediate 8 neighbors.
#' The pixels where the difference is below the defined threshold
#' (\emph{ct}) are aggregated into a single region. The user can define
#' a minimum pixel value using \emph{min.prob} which will limit the range of
#' pixels under evaluation. The result contains a raster with unique
#' values for each segment region region (\emph{$segment}) as well as a
#' data frame (\emph{$stats}) with statistics for each region. The data
#' frame report on the minimum (\emph{min}), maximum (\emph{max}), mean (\emph{mean})
#' and standard deviation (\emph{sd}) of the pixels contained in each region.}
#' @seealso \code{\link{moveModel}} \code{\link{modelApply}}
#' @examples {
#'
#'  require(raster)
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

segRaster <- function(prob, break.point=0.1, min.prob=NULL) {

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
  if (is.null(min.prob)) {pos <- which.max(!is.na(prob))} else {pos <- which.max(prob >= min.prob)}
  
  # evaluate pixel connectivity
  regions <- raster(extent(prob), vals=0, res=res(lc), crs=crs(prob))
  for (r in 1:length(pos)) {
    rp <- rowFromCell(prob, pos[r])
    if (rp > 1) {sr<-rp-1} else {sr<-rp}
    if (rp < nr) {er<-rp+1} else {er<-rp}
    cp <- colFromCell(prob, pos[r])
    if (cp > 1) {sc<-cp-1} else {sc<-cp}
    if (cp < nc) {ec<-cp+1} else {ec<-cp}
    if (max(regions[sr:er,sc:ec])>0) {
      diff <- abs(prob[sr:er,sc:ec]-prob[rp,cp]) <= break.point
      uv <- unique(regions[sr:er,sc:ec][which(diff)])
      uv <- uv[which(uv>0)]
      if (length(uv)>0) {
        mv <- min(uv)
        regions[rp,cp]<-min(uv)
        for (u in 1:length(uv)) {regions[which.max(regions==uv[u])] <- mv}
      } else {regions[rp,cp]<- cellStats(regions,max, na.rm=T)+1}
    } else {regions[rp,cp] <- cellStats(regions,max, na.rm=T)+1}
  }

#-----------------------------------------------------------------------------------#
# 3. derive segment statistics
#-----------------------------------------------------------------------------------#

  # update region id
  uv <- sort(unique(regions))
  uv <- uv[which(uv > 0)]
  uregions = regions
  nr <- length(uv)
  pmn <- vector('numeric', nr) # min
  pmx <- vector('numeric', nr) # max
  pav <- vector('numeric', nr) # mean
  psd <- vector('numeric', nr) # sd
  npx <- vector('numeric', nr) # count
  for (u in 1:nr) {
    pos <- which.max(regions==uv[u])
    uregions[pos] <- u
    pmn[u] <- min(prob[pos], na.rm=T)
    pmx[u] <- max(prob[pos], na.rm=T)
    pav[u] <- mean(prob[pos], na.rm=T)
    psd[u] <- sd(prob[pos], na.rm=T)
    npx[u] <- sum(!is.na(prob[pos]))
  }

  rm(regions, prob)

#-----------------------------------------------------------------------------------#
# 4. return output
#-----------------------------------------------------------------------------------#

  # convert data back to raster
  uregions <- raster(uregions)
  extent(uregions) <- extent(prob)
  res(uregions) <- res(prob)
  crs(uregions) <- crs(prob)

  # build/return data frame
  df <- data.frame(segment=uv, min=pmn, max=pmx, mean=pav, sd=psd, count=npx)

  # return matrix/df
  return(list(segment=uregions, stats=df))

}

