#' @title segRaster
#'
#' @description Connencte-region based raster segmentation.
#' @param prob Object of class "RasterLayer".
#' @param pt Probability threshold. Default is 0.05.
#' @param mp Minimum value. Default is 0.5.
#' @inport raster
#' @return {A list containing a segmented image ("$segment") 
#' and statistics for each segment ("$stats")}
#' @details {The function segments an input layer using a connected 
#' component region labeling approach. After selecting the target 
#' pixels ("mp"), the function estimates the difference
#' in between each pixel and its imediate 8 neighbors. The pixels where the
#' difference is below the defined threshold ("ct") are aggregated into 
#' a single region. The result contains a raster with unique values for 
#' each segment ("$segment") as well as a table ("$stats) which reports 
#' on statistics for their raster values ("min", "max", "mean", "sd").
#' @export
#' @keywords
#' @seealso
#' @aliases
#' @examples \dontrun{
#'  
#' }
#' 
#-----------------------------------------------------------------------------------#

segRaster <- function(prob, pt=0.05, mp=0.5)

#-----------------------------------------------------------------------------------#
# check input variables
#-----------------------------------------------------------------------------------#
  
  if (!exists("prob")) {return('error: "prob" is missing')}
  if (!class(prob)[1]!='RasterLayer') {return('error: "prob" is not a "RasterLayer"')}

#-----------------------------------------------------------------------------------#
# segment regions
#-----------------------------------------------------------------------------------#
  
  # identify usable pixels
  data <- as.matrix(prob)
  pos <- which(data > pt)
  
  # evaluate pixel connectivity
  regions <- as.matrix(data)
  for (r in 1:length(up)) {
    xp <- pos[r] %% nr
    yp <- pos[r] / nr
    if (xp > 1) {sc<-xp-1} else {sc<-xp}
    if (xp < nc) {ec<-xp+1} else {ec<-xp}
    if (yp > 1) {sr<-yp-1} else {sr<-yp}
    if (yp < nr) {er<-yr+1} else {er<-yp}
    if (data[xp,yp]==0) {
      diff <- abs(data[sc:ec,sr:er]-data[xp,yp]) <= ct
      rv <- regions[sc:ec,sr:er]*ct
      if (max(rv)>0) {
        mv <- min(rv)
        uv = unique(rv)
        regions[xp,yp] <- mv
        for (u in 1:length(uv)) {regions[which(regions==uv[u])]<-mv}
      } else {regions[sc:ec,sr:er] <- diff*(max(regions)+1)
    }
  }
  
#-----------------------------------------------------------------------------------#
# derive segment statistics
#-----------------------------------------------------------------------------------#
    
  # update region id
  uv = unique(regions[which(regions>0)])
  uregions = regions
  nr <- length(uv)
  pmn <- vector('numeric', nr) # min
  pmx <- vector('numeric', nr) # max
  pav <- vector('numeric', nr) # mean
  psd <- vector('numeric', nr) # sd
  for (u in 1:nr) {
    pos <- which(regions==uv[u])
    uregions[pos] <- u
    pmn[u] <- min(data[pos])
    pmx[u] <- max(data[pos])
    pav[u] <- mean(data[pos])
    psd[u] <- sd(data[pos])
  }
  
  rm(regions, data)
  
#-----------------------------------------------------------------------------------#
# return output
#-----------------------------------------------------------------------------------#
  
  # convert data back to raster
  uregions <- raster(uregions, rr)
  
  # build/return data frame
  df <- data.frame(min=pmn, max=pmx, mean=pav, sd=psd, stringsAsFactors=F)
  
  # return matrix/df
  return=list(segment=uregions, stats=df)

}