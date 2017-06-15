#' @title moveModel
#'
#' @description modeling of
#' @param xyz matrix provided by labelSample()
#' @param var
#' @param runs
#' @param ma
#' @param tc
#' @param rs
#' @inport caret, raster
#' @return
#' @export
#' @usage
#' @keywords
#' @seealso
#' @aliases
#' @examples \dontrun{
#'
#' }

moveModel <-function(samples=samples, features=features, runs=runs, method=method, control=control) {

  #--------------------------------------------------------------------------------------------------#
  # 1. check input variables and define auxiliary functions
  #--------------------------------------------------------------------------------------------------#

  # check variables in data frame
  if (!exists(samples)) {stop('missing sample information')}
  if (!is.null(samples$x)) {stop('missing x coordinates')}
  if (!is.null(samples$y)) {stop('missing y coordinates')}
  if (!is.null(samples$c)) {stop('missing class information')}
  if (!is.null(samples$r)) {stop('missing region information')}

  # define modeling algorithm and model control method
  if (!exists(method)) {method = 'rf'}
  if (!exists(control)) {control = trainControl(method='oob')}

  # identidy incomplete rows
  cc <- complete.cases(var)
  if (length(cc)!=0) {stop('no usable data in "features" (NA values in all positions')}

  # update variabes (remove incomplete rows)
  input <- input[cc,]

  # function used to select number of pc's
  pcf = function(x) {which((x$sdev^2) > 1)}

  #--------------------------------------------------------------------------------------------------#
  # 2. extract background samples
  #--------------------------------------------------------------------------------------------------#

}

#
#
#
# # build output matrixes
# mean.acc = matrix(0, nr, 2)
#
# # read original data
# input = read.table(paste0('./', population, '/', sfile1), header=F, sep=';')
#
#
#
#
# #----------------------------------------------------------------------------------------------------------#
# # perform n training runs (used to evaluate model consistency)
# for (n in 1:nr) {
#
#   pp = 0
#   cp = 0
#   tp = 0
#   pa = 0
#   ca = 0
#   ta = 0
#
#   # use each unique region for validation
#   for (v in 1:length(uv)) {
#
#     # split absence indices (1/2 split)
#     si = sample(1:length(ai), length(ai))
#     ai.t = ai[si[seq(from=1, to=length(ai), by=2)]]
#     ai.v = ai[si[seq(from=2, to=length(ai), by=2)]]
#
#     # build model for training set
#     ti = which(g!=uv[v] & y==1) # training indices (presence)
#     ci = c(ai.t[sample(1:length(ai.t), length(ti), replace=T)], ti)
#     model = train(x[ci,], as.factor(y[ci]), method='rf', trControl=tc)
#
#     # estimate/store accuracies
#     vi = which(g==uv[v])
#     ci = c(ai.v[sample(1:length(ai.v), length(vi), replace=T)], vi)
#     pred = as.numeric(predict(model, x[ci,]))
#     pp = pp + sum(pred == 1)
#     cp = cp + sum(pred == 1 & y[ci] == 1)
#     tp = tp + sum(y[ci]==1)
#     pa = pa + sum(pred == 2)
#     ca = ca + sum(pred == 2 & y[ci] == 2)
#     ta = ta + sum(y[ci]==2)
#
#     # remove temporary variables
#     rm(model, vi, ti, si, ai.t, ai.v)
#
#   }
#
#   # estimate final accuracy list
#   p = cp / pp
#   r = cp / tp
#   mean.acc[n,1] = 2 * ((p * r) / (p + r))
#   p = ca / pa
#   r = ca / ta
#   mean.acc[n,2] = 2 * ((p * r) / (p + r))
#
#   # remove temporary variables
#   rm(p, r, cp, pp, tp, ca, pa, ta)
#
# }
#
# mean.acc = as.data.frame(mean.acc)
# colnames(mean.acc) = c('Presence', 'Absence')
#
# # accuracy reports
# write.csv(mean.acc, paste0('./', population, '/', population, '_f1-mean.txt'))
