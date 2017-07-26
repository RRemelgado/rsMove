#' @title moveModel
#' 
#' @description Spatially stratified predictive modeling.
#' @param pxy Object of class \emph{SpatialPoinsDataFrame} with presence environmental variables.
#' @param axy Object of class \emph{SpatialPoinsDataFrame} with background environmental variables.
#' @param label Region labels. If missing, "pxy" is assumed as one region.
#' @param method Classification algorithm (see \url{http://topepo.github.io/caret/index.html}. Default is \emph{rf} (Radom Forest).
#' @param control Object derived by \emph{trainControl} (see \code{\link[caret]{trainControl}}). Default used out-of-bag (\emph{oob}) accuracies.
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
#' which is trained using all the samples. This output can be passed to modelApply().}
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
#'  out <- moveModel(pxy=rsQuery, axy=bSamples, label=label)
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------#

moveModel <-function(pxy=pxy, axy=axy, label=NULL, method=NULL, control=NULL, nruns=1) {

#----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables and define auxiliary functions
#----------------------------------------------------------------------------------------------------------------------------------#
  
  # check variables in data frame
  if (crs(pxy)@projargs!=crs(axy)@projargs) {stop('"pxy" and "axy" have different projections')}
  if (ncol(pxy@data)!=ncol(axy@data)) {stop('variables in "pxy" and "axy" have different dimensions')}
  if (min(colnames(pxy@data)==colnames(pxy@data))==0) {stop('one or more variables have different names')}
  
  # check labels
  if (is.null(label)) {label <- vector('numeric', length(pxy))+1}
  if (length(label)!=length(pxy)) {stop('"pxy" and "label" have different lengths')}
  
  # remove duplicates
  cc <- complete.cases(pxy@data)
  pxy <- pxy[cc,]
  label <- label[cc]
  axy <- axy[complete.cases(axy@data),]
  
  # define modeling algorithm and model control method
  if (!is.null(method)) {if (is.null(control)) {stop('"method" defined. Missing "control"')}}
  if (!is.null(control)) {if (is.null(method)) {stop('"control" defined. Missing "method"')}}
  if (is.null(method) & is.null(control)) {
    method <- 'rf'
    control <- trainControl(method='oob')}
  if (!exists('nruns')) {nruns <- 1}

#----------------------------------------------------------------------------------------------------------------------------------#
# 2. define class codes
#----------------------------------------------------------------------------------------------------------------------------------#
  
  i1 <- vector('numeric', length(pxy))+1 # presence class code
  i0 <- vector('numeric', length(axy))+2 # absence class code
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
      model <- train(rbind(pxy@data[i1.t,], axy@data[i0.t,]), as.factor(c(i1[i1.t],i0[i0.t])), method=method, trControl=control)
      
      # estimate/store accuracies
      vi <- c(i1[i1.v], i0[i0.v])
      pred <- as.numeric(predict(model, rbind(pxy@data[i1.v,], axy@data[i0.v,])))
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
    m.ls[[n]] <- train(rbind(pxy@data, axy@data), as.factor(c(i1,i0)), method=method, trControl=control)
    
    # remove temporary variables
    rm(p, r, cp, pp, tp, ca, pa, ta)
    
  }
  
  # write output
  return(list(f1=mean.acc, model=m.ls))
  
}