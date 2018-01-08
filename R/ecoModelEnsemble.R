#' @title ecoModelEnsemble
#'
#' @description Emsemble-based predictive modeling of resource suitability based on presence/absence samples.
#' @param presences
#' @param absences
#' @param sample.label
#' @param img
#' @param img Object of class \emph{RasterLayer}, \emph{RasterStack} or \emph{RasterBrick}.
#' @param fun Function that specifies how the model should be applied.
#' @import raster sp caret rgdal
#' @importFrom stats complete.cases
#' @return A \emph{Raster}.
#' @details {This function is built on the premise that areas which recorded more samples are more ecologically
#' relevant than those that had less samples in comparison. \emph{modelApply} uses the information on validation
#' sample count and provided by \emph{moveModel} to imform on the relative important of the models derived by this
#' same function. First, for each model in this list, \emph{modelApply} derives a spatial prediction and adds it to
#' a \emph{RasterBrick}. Then, the number of samples used for validation
#' of models
#'
#'
#'
#' If this contains a
#' list of models from multiple runs, the function creates a stack of predictions
#' and summarizes it on a pixel-by-pixel basis using a weighted mean. The weights
#' are defined by the average performance for \emph{presence} and \emph{background}
#' samples in each iteration.}
#' @seealso \code{\link{segRaster}} \code{\link{stratModel}}
#' @examples \dontrun{
#'
#'  require(raster)
#'
#'  # read remote sensing data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'tc.*tif', full.names=TRUE)
#'  rsStk <- stack(file)
#'
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130805-20130811.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(rsStk))
#'
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=moveData,img=rsStk, remove.dup=TRUE)
#'
#'  # identify unique sample regions
#'  label <- labelSample(xy=rsQuery, rad=90, npx=1, pxr=rsStk)
#'
#'  # select background samples
#'  ind <- which(label>0) # selected samples
#'  bSamples <- backSample(xy=moveData[ind,], rid=label[ind], img=rsStk, nb=4000, method='pca')
#'
#'  # derive model predictions
#'  fun <- function(x,y) {train(x, y, method="rf", trControl=trainControl(method='oob'))}
#'  p.model <- moveModel(p.data=rsQuery@data, a.data=bSamples@data, label=label, fun=fun, nruns=1)
#'
#'  # derive prediction from model ensemble
#'  fun <- function(x,y) {predict(x, y, type='prob')[[1]]$`1`}
#'  prob <- modelApply(p.model, rsStk, fun=fun)
#'
#'  # see output
#'  plot(prob)
#'
#' }
#' @export

#--------------------------------------------------------------------------------#

modelApply <- function(model, img, fun=NULL) {

#--------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------#

  if (!exists('model')) {stop('error: "model" is missing')}
  if (is.null(model$f1)) {stop('error: "model" is not a valid input')}
  if (!exists('img')) {stop('error: "img" is missing')}
  if (is.null(model$model)) {stop('error: "model" is not a valid input')}
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    stop('error: "img" is not a valid raster layer')}
  if (is.null(fun)) {fun <- function(x,y) {predict(x, y, type='prob')[[1]]$`1`}}
  if (!is.null(fun)) {if (!is.function(fun)) {stop('"fun" is not a function')}}

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
    tmp[cc] <- fun(model$model[m], idata[cc,]) # presence prob.
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
