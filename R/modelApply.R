#' @title modelApply
#'
#' @description Apply a model or an ensemble of models to raster data.
#' @param model List object as provided by \emph{moveModel()}.
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @import raster sp caret rgdal
#' @importFrom stats complete.cases
#' @return A \emph{Raster}.
#' @details {The function uses the output of \emph{moveModel()}. If this contains a 
#' list of models from multiple runs, the function creates a stack of predictions 
#' and summarizes it on a pixel-by-pixel basis using a weighted mean. The weights 
#' are defined by the average performance for \emph{presence} and \emph{background} 
#' samples in each iteration.}
#' @seealso \code{\link{segRaster}} \code{\link{moveModel}}
#' @examples \dontrun{
#'  
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'  
#'  # read remote sensing data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'  
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130805-20130811.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(rsStk))
#'
#'  # extract samples
#'  ot = as.Date(moveData@data$date)
#'  samples <- sampleMove(xy=moveData, ot=ot, error=10, method='m')
#'  
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=samples,img=rsStk, rd=TRUE)
#'  
#'  # identify unique sample regions
#'  label <- labelSample(xy=rsQuery, rad=90, npx=1, pxr=rsStack)
#'  
#'  # select background samples
#'  ind <- which(label>0) # selected samples
#'  bSamples <- backSample(xy=moveData[ind,], rid=label[ind], img=rsStk, method='pca')
#'  
#'  # derive model predictions
#'  p.model <- moveModel(pxy=rsQuery, axy=bSamples, label=label)
#'  
#'  # derive prediction from model ensemble
#'  prob <- modelApply(p.model, rsStack)
#'  
#'  # see output
#'  plot(prob)
#'  
#' }
#' @export

#--------------------------------------------------------------------------------#

modelApply <- function(model, img) {
  
#--------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------#
  
  if (!exists('model')) {stop('error: "model" is missing')}
  if (is.null(model$f1)) {stop('error: "model" is not a valid input')}
  if (!exists('img')) {stop('error: "img" is missing')}
  if (is.null(model$model)) {stop('error: "model" is not a valid input')}
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    stop('error: "img" is not a valid raster layer')}
  
#---------------------------------------------------------------------------------------#
# 2. read/prepare raster data
#---------------------------------------------------------------------------------------#
  
  idata <- getValues(img) # get values
  cc <- complete.cases(idata) # rows with no NA's
  tmp <- idata[,1] # reference vector to assign predicted values
  tmp[!cc] <- NA # set missing values to NA
  rb <- brick(img, nl=length(model$model)) # will contain prob. images
  
#---------------------------------------------------------------------------------------#
# 3. apply models (build ensemble)
#---------------------------------------------------------------------------------------#
  
  for (m in 1:length(model$model)) {
    tmp[cc] <- predict(model$model[m], idata[cc,], type='prob')[[1]]$`1` # presence prob.
    rb[[m]] <- setValues(img[[1]], tmp) # translate prediction into a raster
  }
  
#---------------------------------------------------------------------------------------#
# 4. summarize model ensemble (if required)
#---------------------------------------------------------------------------------------#
  
  if (length(model$model) > 1) {
    w <- apply(model$f1, 1, mean) # derive weights
    sf <- function(x) {sum(x*w) / sum(w)} # weighted mean
    rb <- calc(rb, sf)} # estimate final prediction
  return(rb)
  
}