#' @title predictResources
#'
#' @description Spatially stratified predictive modeling of resource suitability based on presence/absence samples.
#' @param presence.data Object of class \emph{data.frame} with environmental variables for presence samples.
#' @param absence.data Object of class \emph{data.frame} with environmental variables for background samples.
#' @param sample.label Numeric or character vector with sample region labels. If missing, "presence.data" is assumed as one region.
#' @param env.data Object of class \emph{RasterStack} or \emph{RasterBrick} with environmental variables in \emph{presence.data} and \emph{absence.data}.
#' @importFrom stats complete.cases
#' @importFrom caret train trainControl
#' @importFrom raster predict
#' @return A \emph{list}.
#' @references \href{10.1002/rse2.70}{Remelgado, R., Leutner, B., Safi, K., Sonnenschein, R., Kuebert, C. and Wegmann, M. (2017), Linking animal movement and remote sensing - mapping resource suitability from a remote sensing perspective. Remote Sens Ecol Conserv.}
#' @details {Modeling of resource suitability using animal movement data following a recent paper (Remelgado et al, 2017). For each
#' unique label in \emph{sample.label}, the function keeps it for validation and uses the remaining samples for training. Then, the
#' function evaluates the performance of this model reporting (internally) on the number of true positives, false positives and the
#' number of validation and predicted cases for both presences and absences. Once all sample regions are used for validation, the
#' reported values are summed and used to derive the F1-measure. The F1-measure is estimated as \emph{2 * (P * R) / (P + R)} where
#' \emph{P} is the Precision (ratio of true positives within the number of predicted values) and \emph{R} is the Recall (ratio of
#' true positives within the number of validation samples). As a consequence, rather than reporting on an average performance, the
#' final performance assessment reported by \emph{predictResources} depicts an objective picture on how the model performed among the
#' different sets sample regions. This metric is provided for presences (\emph{presence.data}) and absences ({\emph{absence.data}})
#' separately offering an overview on the stability of the model. This analysis is performed using a Random Forest model as provided
#' within the \code{\link[caret]{train}} function of the caret package. The final predictive model is then derived with all samples.
#' The output of \emph{predictResources} is a list object consisting of:
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
#'  shortMove <- moveReduce(xy=shortMove, obs.time=obs.time, img=rsStk)$points
#'
#'  # retrieve remote sensing data for samples
#'  rsQuery <- extract(rsStk, shortMove)
#'
#'  # identify unique sample regions
#'  label <- labelSample(xy=shortMove, agg.radius=90, pixel.res=rsStk)
#'
#'  # select background samples
#'  ind <- which(!is.na(label)) # selected samples
#'  bSamples <- backSample(xy=shortMove[ind,], region.id=label[ind],
#'  img=rsStk, sampling.method='pca')
#'
#'  # derive model predictions
#'  out <- predictResources(presence.data=rsQuery,
#'  absence.data=bSamples@data, sample.label=label, env.data=rsStk)
#'
#' }
#' @export

#----------------------------------------------------------------------------------------------------------------------------------#

predictResources <-function(presence.data=presence.data, absence.data=absence.data, sample.label=NULL, env.data=NULL) {

#----------------------------------------------------------------------------------------------------------------------------------#
# 1. check input variables and define auxiliary functions
#----------------------------------------------------------------------------------------------------------------------------------#

  # check variables in data frame
  if (ncol(presence.data)!=ncol(absence.data)) {stop('different number of variables in "presence.data" and "absence.data"')}
  if (min(colnames(presence.data)==colnames(absence.data))==0) {stop('column names of "presence.data" and "absence.data" differ')}

  # check labels
  if (is.null(sample.label)) {sample.label <- vector('numeric', length(presence.data))+1}
  if (length(sample.label)!=nrow(presence.data)) {stop('"presence.data" and "sample.label" have different lengths')}

  # remove duplicates
  cc <- complete.cases(presence.data)
  presence.data <- presence.data[cc,]
  sample.label <- sample.label[cc]
  absence.data <- absence.data[complete.cases(absence.data),]

  # check environmental data
  if (!is.null(env.data)) {
    if (!class(env.data)[1]%in%c("RasterStack", "RasterBrick")) {stop('"env.data" is not a valid raster object')}
    if (nlayers(env.data)!=ncol(presence.data)) {stop('"env.data" has a different amount of variables from the training data')}
    if (length(colnames(presence.data))!=length(names(env.data))) {stop('the variable names of "env.data" are differ from the predictive data')}}

  # training metric
  tc <- trainControl(method='oob')

#----------------------------------------------------------------------------------------------------------------------------------#
# 2. define class codes
#----------------------------------------------------------------------------------------------------------------------------------#

  i1 <- vector('numeric', nrow(presence.data))+1 # presence class code
  i0 <- vector('numeric', nrow(absence.data))+2 # absence class code
  uv <- unique(sample.label) # unique region id's

#----------------------------------------------------------------------------------------------------------------------------------#
# 3. Build model
#----------------------------------------------------------------------------------------------------------------------------------#

  # initiate outputs
  val.set <- data.frame(region=vector(class(sample.label), length(uv)), count=vector('numeric', length(uv)))
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
    i1.t <- which(sample.label!=uv[v])
    i1.v <- which(sample.label==uv[v])

    # split absence indices (1/2 split)
    si <- sample(1:length(i0), length(i0))
    i0.t <- i0[si[seq(from=1, to=length(i0), by=2)]]
    i0.v <- i0[si[seq(from=2, to=length(i0), by=2)]]

    # build model for training set
    m.ls[[v]] <- train(rbind(presence.data[i1.t,], absence.data[i0.t,]),
                       as.factor(c(i1[i1.t],i0[i0.t])), method="rf", trControl=tc)

    # estimate/store accuracies
    vi <- c(i1[i1.v], i0[i0.v])
    pred <- as.numeric(predict(m.ls[[v]], rbind(presence.data[i1.v,], absence.data[i0.v,])))
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
  model <- train(rbind(presence.data, absence.data), as.factor(c(i1,i0)), method="rf", trControl=tc)
  if (!is.null(env.data)) {prob <-  calc(env.data, function(x){predict(model, x, type='prob')$'1'})} else {prob <- NULL}
  return(list(f1=acc, validation=val.set, iteration.models=m.ls, final.model=model, probabilities=prob))

}
