#' @title segRaster
#'
#' @description {Connected-region based raster segmentation that preserves spatial gradients.}
#' @param x Object of class \emph{RasterLayer}.
#' @param break.point Difference threshold. Default is 0.05.
#' @param min.value Minimum value. Default is 0.5.
#' @importFrom raster raster extent crs res
#' @importFrom stats sd
#' @return A list object.
#' @details {The function segments an input raster layer (\emph{x}) using a
#' connected component region labeling approach. For each pixel, the function
#' estimates the difference between it and its immediate neighbors. If the
#' difference is below the threshold defined by \emph{break.point} these are
#' aggregated into a single region. Moreover, the user can define a minimum pixel
#' value using \emph{min.value} which will ignore all pixels below that value. The
#' output of this function consists of:
#' \itemize{
#'  \item{\emph{regions} - Region raster image.}
#'  \item{\emph{stats} - Basic statistics for each pixel region.}}}
#' @seealso \code{\link{predictResources}}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # load example probability image
#'  file <- system.file('extdata', 'probabilities.tif', package="rsMove")
#'  r <- raster(file)
#'
#'  # segment probabilities
#'  rs <- segRaster(r)
#'
#' }
#' @export

#-----------------------------------------------------------------------------------#

segRaster <- function(x, break.point=0.1, min.value=0.5) {

  #-----------------------------------------------------------------------------------#
  # 1. check input variables
  #-----------------------------------------------------------------------------------#

  if (class(x)[1]!='RasterLayer') {stop('"x" is not a "RasterLayer"')}
  if (is.null(min.value)) {min.value <- 0}

  #-----------------------------------------------------------------------------------#
  # 2. segment regions
  #-----------------------------------------------------------------------------------#

  # raster dimensions
  nr <- dim(x)[1]
  nc <- dim(x)[2]

  # identify usable pixels
  pos <- which.max(x >= min.value)

  # evaluate pixel connectivity
  regions <- x
  regions[] <- 0

  tmp <- sapply(pos, function(p) {

    rp <- rowFromCell(x, p)
    cp <- colFromCell(x, p)

    if (cp > 1) {sc<-cp-1} else {sc<-cp}
    if (cp < nc) {ec<-cp+1} else {ec<-cp}
    if (rp > 1) {sr<-rp-1} else {sr<-rp}
    if (rp < nr) {er<-rp+1} else {er<-rp}

    if (max(regions[sr:er,sc:ec]) > 0) {

      diff <- abs(x[sr:er,sc:ec]-x[rp,cp]) <= break.point & is.finite(x[sr:er,sc:ec])
      uv <- unique(c(regions[sr:er,sc:ec][diff]))
      uv <- uv[which(uv>0)]
      if (length(uv>0)) {
        mv <- min(uv)
        regions[rp,cp]<-min(uv)
      } else {regions[rp,cp] <<- cellStats(regions, max)+1}

    } else {regions[rp,cp] <<- cellStats(regions, max)+1}

  })

  #-----------------------------------------------------------------------------------#
  # 3. derive segment statistics
  #-----------------------------------------------------------------------------------#

  # update region id
  uv <- sort(unique(regions[which.max(regions > 0)]))
  nr <- length(uv)
  pmn <- vector('numeric', nr) # min
  pmx <- vector('numeric', nr) # max
  pav <- vector('numeric', nr) # mean
  psd <- vector('numeric', nr) # sd
  npx <- vector('numeric', nr) # count
  for (u in 1:nr) {
    pos <- which.max(regions==uv[u])
    regions[pos] <- u
    pmn[u] <- min(x[pos], na.rm=T)
    pmx[u] <- max(x[pos], na.rm=T)
    pav[u] <- mean(x[pos], na.rm=T)
    psd[u] <- sd(x[pos], na.rm=T)
    npx[u] <- sum(!is.na(x[pos]))
  }

  #-----------------------------------------------------------------------------------#
  # 4. return output
  #-----------------------------------------------------------------------------------#

  # build/return data frame
  df <- data.frame(segment=uv, min=pmn, max=pmx, mean=pav, sd=psd, count=npx)

  # return matrix/df
  return(list(regions=regions, stats=df))

}
