#' @title moveModel
#' 
#' @description Spatially stratified predictive modeling.
#' @param p.data Object of class \emph{data.frame} with environmental variables for presence samples.
#' @param a.data Object of class \emph{data.frame} with environmental variables for background samples.
#' @param label Region labels. If missing, "p.data" is assumed as one region.
#' @param fun A function with the modeling algorithm to use.
#' @param nruns Number of runs. Default is 1.
#' @import caret raster rgdal
#' @importFrom stats complete.cases
#' @return A \emph{list}.
#' @details {For n iterations, where n is the number of unique sample regions, 
#' the function uses one sample region for validation while the remaining ones 
#' are used for training. The background samples are split randomly at each 
#' iteration. The final accuracy, provided as a F1-score for both presence and 
#' background sampels, is derived from the total of true and false positives 
#' (\emph{$f1}). Additionaly, for each run, the function returns a model (\emph{$model}) 
#' which is trained using all the samples. This output can be passed to modelApply(). 
#' By default, the function uses a Random Forest classifier. However, the a user specified 
#' function can be passed using \emph{fun}.}
#' @seealso \code{\link{sampleMove}} \code{\link{labelSample}} \code{\link{backSample}} \code{\link{modelApply}} \code{\link[caret]{train}}
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
#'  # retrieve remote sensing data for samples
#'  rsQuery <- dataQuery(xy=moveData,img=rsStk, remove.dup=TRUE)
#'  
#'  # identify unique sample regions
#'  label <- labelSample(xy=rsQuery, rad=3000, pxr=rsStk)
#'  
#'  # select background samples
#'  ind <- which(label>0) # selected samples
#'  bSamples <- backSample(xy=moveData[ind,], rid=label[ind], img=rsStk, method='pca')
#'  
#'  # derive model predictions
#'  fun <- function(x,y) {train(x, y, method="rf", trControl=trainControl(method='oob'))}
#'  out <- moveModel(p.data=rsQuery@data, a.data=bSamples@data, label=label, fun=fun, nruns=1)
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------#

moveModel <-function(p.data=p.data, a.data=a.data, label=NULL, fun=fun, nruns=1) {

#----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables and define auxiliary functions
#----------------------------------------------------------------------------------------------------------------------------------#
  
  # check variables in data frame
  if (ncol(p.data)!=ncol(a.data)) {stop('different number of variables in "p.data" and "a.data"')}
  if (min(colnames(p.data)==colnames(a.data))==0) {stop('column names of "p.data" and "a.data" differ')}
  
  # check labels
  if (is.null(label)) {label <- vector('numeric', length(p.data))+1}
  if (length(label)!=nrow(p.data)) {stop('"p.data" and "label" have different lengths')}
  
  # remove duplicates
  cc <- complete.cases(p.data)
  p.data <- p.data[cc,]
  label <- label[cc]
  a.data <- a.data[complete.cases(a.data),]
  
  # define modeling algorithm and model control method
  if (is.null(fun)) {fun <- function(x,y) {train(x, y, method="rf", trControl=trainControl(method='oob'))}}
  if (!is.null(fun)) {if (!is.function(fun)) {stop('"fun" is not a function')}}
  
  if (!is.numeric(nruns)) {stop('"nruns" is not numeric')}
  if (length(nruns)>1) {stop('"nrus" should be a single numeric element')}
  
#----------------------------------------------------------------------------------------------------------------------------------#
# 2. define class codes
#----------------------------------------------------------------------------------------------------------------------------------#
  
  i1 <- vector('numeric', nrow(p.data))+1 # presence class code
  i0 <- vector('numeric', nrow(a.data))+2 # absence class code
  uv <- unique(label) # unique region id's
  
#----------------------------------------------------------------------------------------------------------------------------------#
# 3. Build model
#----------------------------------------------------------------------------------------------------------------------------------#
  
  # initiate outputs
  mean.acc <- data.frame(presence=matrix(0,nruns), background=matrix(0,nruns))
  m.ls <- vector('list', nruns) # model list
  
  # perform n training runs (used to evaluate model consistency)
  for (n in 1:nruns) {
    
    pp <- 0
    cp <- 0
    tp <- 0
    pa <- 0
    ca <- 0
    ta <- 0
    
    # use each unique region for validation
    for (v in 1:length(uv)) {
      
      # split presences
      i1.t <- which(label!=uv[v])
      i1.v <- which(label==uv[v])
      
      # split absence indices (1/2 split)
      si <- sample(1:length(i0), length(i0))
      i0.t <- i0[si[seq(from=1, to=length(i0), by=2)]]
      i0.v <- i0[si[seq(from=2, to=length(i0), by=2)]]
      
      # build model for training set
      model <- fun(rbind(p.data[i1.t,], a.data[i0.t,]), as.factor(c(i1[i1.t],i0[i0.t])))
      
      # estimate/store accuracies
      vi <- c(i1[i1.v], i0[i0.v])
      pred <- as.numeric(predict(model, rbind(p.data[i1.v,], a.data[i0.v,])))
      pp <- pp + sum(pred == 1)
      cp <- cp + sum(pred == 1 & vi == 1)
      tp <- tp + sum(vi==1)
      pa <- pa + sum(pred == 2)
      ca <- ca + sum(pred == 2 & vi == 2)
      ta <- ta + sum(vi==2)
      
      # remove temporary variables
      rm(model, i1.t, i1.v, i0.t, i0.v, vi, si)
      
    }
    
    # estimate final accuracy list
    p <- cp / pp
    r <- cp / tp
    mean.acc$presence[n] = 2 * ((p * r) / (p + r))
    p <- ca / pa
    r <- ca / ta
    mean.acc$background[n] = 2 * ((p * r) / (p + r))
    
    # update model list
    m.ls[[n]] <- fun(rbind(p.data, a.data), as.factor(c(i1,i0)))
    
    # remove temporary variables
    rm(p, r, cp, pp, tp, ca, pa, ta)
    
  }
  
  # write output
  return(list(f1=mean.acc, model=m.ls))
  
}