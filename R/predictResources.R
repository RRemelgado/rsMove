#' @title predictResources
#'
#' @description Spatially stratified predictive modeling of resource suitability based on presence/absence samples.
#' @param x Object of class \emph{data.frame} with environmental variables for presence samples.
#' @param y Object of class \emph{data.frame} with environmental variables for background samples.
#' @param z \emph{Numeric} or \emph{character} vector with sample region labels. If missing, \emph{x} is assumed as being one region.
#' @param env.data Object of class \emph{RasterStack} or \emph{RasterBrick} with environmental variables in \emph{x} and \emph{y}.
#' @importFrom stats complete.cases
#' @importFrom caret train trainControl
#' @importFrom raster calc nlayers
#' @return A \emph{list}.
#' @references \href{10.1002/rse2.70}{Remelgado, R., Leutner, B., Safi, K., Sonnenschein, R., Kuebert, C. and Wegmann, M. (2017), Linking animal movement and remote sensing - mapping resource suitability from a remote sensing perspective. Remote Sens Ecol Conserv.}
#' @details {Modeling of resource suitability using animal movement data following the method of Remelgado et al (2017). Each
#' unique label in \emph{z} is kept for validation while the remaining samples are used for training. Then, the function evaluates
#' the performance of this model reporting (internally) on the number of true positives, false positives and the number of cases for
#' both presences and absences. Once all sample regions are used for validation, the reported values are summed and used to derive a
#' F1-measure. The F1-measure is estimated as \emph{2 * (P * R) / (P + R)} where \emph{P} is the Precision (ratio of true positives
#' within the number of predicted values) and \emph{R} is the Recall (ratio of true positives within the number of validation samples).
#' As a consequence, rather than reporting on an average performance, the final performance assessment reported by \emph{predictResources}
#' depicts an objective picture on how the model performed among the different sets sample regions. This metric is provided for presences
#' (\emph{x}) and absences ({\emph{y}}) separately informing on the stability of the model. This analysis is performed using a Random Forest
#' model as provided within the \code{\link[caret]{train}} function of the caret package. The final predictive model is then derived with all
#' samples. The output of the function is a list object consisting of:
#' \itemize{
#'  \item{\emph{f1} - \emph{data.frame} with final F1-measure for presences and absences.}
#'  \item{validation - \emph{data.frame} with region identifiers and validation sample count at each iteration.}
#'  \item{\emph{iteration.models} - List of models estimated at each iteration.}
#'  \item{\emph{final.model} - Final predictive model based on all samples.}
#'  \item{\emph{probabilities} - Predicted probability image. Given if \emph{env.data} is set.}}}
#' @seealso \code{\link{sampleMove}} \code{\link{labelSample}} \code{\link{backSample}} \code{\link[caret]{train}}
#' @examples \dontrun{
#'
#'  require(rgdal)
#'  require(raster)
#'  require(sp)
#'
#'  # read remote sensing data
#'  file <- list.files(system.file('extdata', '', package="rsMove"), 'ndvi.tif', full.names=TRUE)
#'  r.stk <- stack(file)
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # observation time
#'  obs.time <- strptime(paste0(shortMove@data$date, ' ', shortMove@data$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # remove redundant samples
#'  shortMove <- moveReduce(shortMove, r.stk, obs.time)$points
#'
#'  # retrieve remote sensing data for samples
#'  rsQuery <- extract(r.stk, shortMove)
#'
#'  # identify unique sample regions
#'  label <- labelSample(shortMove, r.stk, agg.radius=30)
#'
#'  # select background samples
#'  bSamples <- backSample(shortMove, r.stk, label, sampling.method='pca')
#'
#'  # derive model predictions
#'  out <- predictResources(rsQuery, bSamples@data, label, env.data=r.stk)
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------#

predictResources <-function(x, y, z, env.data=NULL) {

#----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables and define auxiliary functions
#----------------------------------------------------------------------------------------------------------------------------------#

  # check variables in data frame
  if (ncol(x)!=ncol(y)) {stop('different number of variables in "x" and "y"')}
  if (min(colnames(x)==colnames(y))==0) {stop('column names of "x" and "y" differ')}

  # check labels
  if (length(z)!=nrow(x)) {stop('"x" and "z" have different lengths')}
  uv <- unique(z) # unique region id's
  if (length(uv) == 1) {stop('"z" only has one unique element')}

  # remove duplicates
  cc <- complete.cases(x)
  x <- x[cc,]
  z <- z[cc]
  y <- y[complete.cases(y),]

  # check environmental data
  if (!is.null(env.data)) {
    if (!class(env.data)[1]%in%c("RasterStack", "RasterBrick")) {stop('"env.data" is not a valid raster object')}
    if (nlayers(env.data)!=ncol(x)) {stop('"env.data" has a different amount of variables from the training data')}
    if (length(colnames(x))!=length(names(env.data))) {stop('the variable names of "env.data" are differ from the predictive data')}}

  # training metric
  tc <- trainControl(method='oob')

#----------------------------------------------------------------------------------------------------------------------------------#
# 2. define class codes
#----------------------------------------------------------------------------------------------------------------------------------#

  i1 <- vector('numeric', nrow(x))+1 # presence class code
  i0 <- vector('numeric', nrow(y))+2 # absence class code

#----------------------------------------------------------------------------------------------------------------------------------#
# 3. Build model
#----------------------------------------------------------------------------------------------------------------------------------#

  # initiate outputs
  val.set <- data.frame(region=vector(class(z), length(uv)), count=vector('numeric', length(uv)))
  m.ls <- vector('list', length(uv)) # model list

  pp <- 0
  cp <- 0
  tp <- 0
  pa <- 0
  ca <- 0
  ta <- 0

  # use each unique region for validation
  for (v in 1:length(uv)) {

    # split presences
    i1.t <- which(z!=uv[v])
    i1.v <- which(z==uv[v])

    # split absence indices (1/2 split)
    si <- sample(1:length(i0), length(i0))
    i0.t <- i0[si[seq(from=1, to=length(i0), by=2)]]
    i0.v <- i0[si[seq(from=2, to=length(i0), by=2)]]

    # build model for training set
    m.ls[[v]] <- train(rbind(x[i1.t,], y[i0.t,]), as.factor(c(i1[i1.t],i0[i0.t])), method="rf", trControl=tc)

    # estimate/store accuracies
    vi <- c(i1[i1.v], i0[i0.v])
    pred <- as.numeric(predict(m.ls[[v]], rbind(x[i1.v,], y[i0.v,])))
    pp <- pp + sum(pred == 1) # class predictions (1)
    cp <- cp + sum(pred == 1 & vi == 1) # correct class predictions (1)
    tp <- tp + sum(vi==1) # validation set class samples (1)
    pa <- pa + sum(pred == 2) # class predictions (2)
    ca <- ca + sum(pred == 2 & vi == 2) # correct class predictions (2)
    ta <- ta + sum(vi==2) # validation set class samples (2)

    # information on validation sample set
    val.set$region[v] = uv[v] # unique identifyer
    val.set$count[v] = length(i1.v) # sample count

    # remove temporary variables
    rm(i1.t, i1.v, i0.t, i0.v, vi, si)

  }

  # estimate final accuracy list
  p <- cp / pp
  r <- cp / tp
  ap = 2 * ((p * r) / (p + r))
  p <- ca / pa
  r <- ca / ta
  aa = 2 * ((p * r) / (p + r))
  acc <- data.frame(presence=ap, absence=aa)

  # remove temporary variables
  rm(p, r, cp, pp, tp, ca, pa, ta, ap, aa)

  # write output
  model <- train(rbind(x, y), as.factor(c(i1,i0)), method="rf", trControl=tc)
  if (!is.null(env.data)) {prob <-  calc(env.data, function(x){stats::predict(model, x, type='prob')$'1'})} else {prob <- NULL}
  return(list(f1=acc, validation=val.set, iteration.models=m.ls, final.model=model, probabilities=prob))

}
