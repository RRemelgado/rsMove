#' @title modelApply
#'
#' @description Apply a model or an ensemble of models to raster data.
#' @param model Model list as provided by moveModel.
#' @param img Object of class "RasterLayer", "RasterStack" or "RasterBrick".
#' @inport raster sp caret
#' @return Raster of predicted probabilities.
#' @details {The function uses the output of moveModel(). If this contains a 
#' list of models from multiple runs, the function creates a stack of precitions 
#' and summarizes it on a pixel-by-pixel basis using a weighted mean. The weights 
#' are defined by the average performance for "presence" and "background" samples 
#' in each iteration.}
#' @seealso \code{\link{segRaster}} \code{\link{moveModel}}
#' @examples \dontrun{
#'
#' }
#' @export

#--------------------------------------------------------------------------------#

modelApply <- function(model, img) {
  
#--------------------------------------------------------------------------------#
# 1. check input variables
#--------------------------------------------------------------------------------#
  
  if (!exists('model')) {'error: "model" is missing'}
  if (is.null(model$f1)) {return('error: "model" is not a valid input')}
  if (!exists('img')) {'error: "img" is missing'}
  if (is.null(model$model)) {return('error: "model" is not a valid input')}
  if (!class(img)[1]%in%c('RasterLayer', 'RasterStack', 'RasterBrick')) {
    return('error: "img" is not a valid raster layer')}
  
#---------------------------------------------------------------------------------------#
# 2. read/prepare raster data
#---------------------------------------------------------------------------------------#
  
  idata <- getValues(img) # get values
  cc <- complete.cases(idata) # rows with no NA's
  tmp <- idata[,1] # reference vector to assign predicted values
  tmp[!cc] <- NA # set missing values to NA
  rb <- brick(extent(img), crs=crs(img), nl=length(model$model)) # brick of predictions
  res(rb) <- res(img) # adjust brick resolution (automaticaly fixes dimensions)
  
#---------------------------------------------------------------------------------------#
# 3. apply models (build ensemble)
#---------------------------------------------------------------------------------------#
  
  for (m in 1:length(model$model)) {
    tmp[cc] <- predict(model, idata[cc,], type='prob')$'1' # presence prob.
    rb[[m]] <- setValues(img[[1]], tmp) # translate prediction into a raster
  }
  
#---------------------------------------------------------------------------------------#
# 4. summarize model ensemble (if required)
#---------------------------------------------------------------------------------------#
  
  if (length(model$model) > 1) {
    w <- apply(model$f1, 1, mean) # derive weights
    sf <- function(x) {if(sum(is.na(x))==0) {sum(x*w) / sum(w)}} # weighted mean
    rb <- calc(rb, sf)} # estimate final prediction
  return(rb)
  
}