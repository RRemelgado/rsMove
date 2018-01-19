#' @title segRaster
#'
#' @description {Connected-region based raster segmentation that preserves spatial gradients.}
#' @param img Object of class \emph{RasterLayer}.
#' @param break.point Difference threshold. Default is 0.05.
#' @param min.value Minimum value. Default is 0.5.
#' @importFrom raster raster extent crs res
#' @importFrom stats sd
#' @return A list object.
#' @details {The function segments an input raster layer (\emph{img}) using a
#' connected component region labeling approach. For each pixel, the function
#' estimates the difference between it and its imidiate neighbors. If the
#' difference is below the threshold defined by \emph{break.point} these are
#' aggregated into a single region. Moreover, the user can define a minimum pixel
#' value using \emph{min.value} which will ignore all pixels below that value. The
#' output of this function consists of:
#' \itemize{
#'  \item{\emph{regions} - Region raster image.}
#'  \item{\emph{stats} - Basic statistics for each pixel region.}}}
#' @seealso \code{\link{predictResources}}
#' @examples {
#'
#'  require(raster)
#'
#'  # load example probability image
#'  file <- system.file('extdata', 'konstanz_probabilities.tif', package="rsMove")
#'  r <- raster(file)
#'
#'  # segment probabilities
#'  rs <- segRaster(r)
#'
#' }
#' @export

#-----------------------------------------------------------------------------------#

segRaster <- function(img, break.point=0.1, min.value=0.5) {

  #-----------------------------------------------------------------------------------#
  # 1. check input variables
  #-----------------------------------------------------------------------------------#

  if (class(img)[1]!='RasterLayer') {stop('"img" is not a "RasterLayer"')}
  if (is.null(min.value)) {min.value <- 0}

  #-----------------------------------------------------------------------------------#
  # 2. segment regions
  #-----------------------------------------------------------------------------------#

  # raster dimensions
  nr <- dim(img)[1]
  nc <- dim(img)[2]

  # identify usable pixels
  data <- as.matrix(img)
  pos <- which(data >= min.value)

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
      diff <- abs(data[sr:er,sc:ec]-data[rp,cp]) <= break.point &
        is.finite(data[sr:er,sc:ec])
      uv <- unique(c(regions[sr:er,sc:ec][diff]))
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
  uv <- sort(unique(regions[which(regions>0)]))
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
  uregions <- raster(uregions)
  extent(uregions) <- extent(img)
  res(uregions) <- res(img)
  crs(uregions) <- crs(img)

  # build/return data frame
  df <- data.frame(segment=uv, min=pmn, max=pmx, mean=pav, sd=psd, count=npx)

  # return matrix/df
  return(list(regions=uregions, stats=df))

}
