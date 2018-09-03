#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#' @title plausibilityTest
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Quantifies and plots the distribution of pixels within a mask over a reference categorical raster object.
#' @param x Object of class \emph{RasterLayer} and \emph{RasterStack}.
#' @param y Object of class \emph{RasterLayer}.
#' @param class.labels Labels of classes in \emph{y} provided as a character vector.
#' @return A \emph{list}.
#' @details {For each lazer in \emph{x}, (e.g. classification mask) the function returns the absolute and relative count of non-NA
#' pixels within each unique value of \emph{y} (e.g. land cover map). Then, the results for each layer are compared in a combined
#' plot. The output of the function is a list consisting of:
#'  \itemize{
#'  \item{\emph{absolute.count} - Absolute pixel count for each layer of \emph{x} overlapping with each value of \emph{y}.}
#'  \item{\emph{relative.count} - Relative pixel count for each layer of \emph{x} overlapping with each value of \emph{y}.}
#'  \item{\emph{relative.plot} - Plot comparing the relative pixel count of the layers in \emph{x} within each value of \emph{y}.}}}
#' @importFrom raster crs extract nlayers
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw ylim theme xlab ylab scale_fill_manual facet_wrap element_blank unit
#' @importFrom stats as.formula
#' @importFrom grDevices rainbow
#' @examples {
#'
#'  require(raster)
#'
#'  # load example probability image
#'  file <- system.file('extdata', 'probabilities.tif', package="rsMove")
#'  p <- raster(file) > 0.5
#'
#'  # land cover map
#'  lc <- raster(system.file('extdata', 'landCover.tif', package="rsMove"))
#'
#'  # segment probabilities
#'  pt <- plausibilityTest(p, lc)
#'
#'  # show plot
#'  pt$relative.plot
#'
#'  # see relative sample count
#'  head(pt$relative.count)
#'
#' }
#' @export
#'
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

plausibilityTest <- function(x, y, class.labels=NULL) {

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

  # check sample layer
  if (!class(x)[1] %in% c("RasterLayer", "RasterStack")) {stop('"x" is not a raster object')}

  # check reference layer
  if (class(y)[1]!="RasterLayer") {stop('"y" is not a valid raster object')}
  if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "y" have different projections')}
  if (sum(dim(x)[1:2]-dim(y)[1:2])!=0) {stop('"x" & "y" have different dimensions')}
  n.runs <- nlayers(x)

  # check auxiliary information
  if (!is.null(class.labels)) {
    if (!is.character(class.labels)) {'"class.labels" is not a "character" vector'}}

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Extract unique cases in "y" and check "class.labels"
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

  # unique raster values
  unique.values <- sort(unique(y))
  unique.values <- unique.values[!is.na(unique.values)]

  # define class labels and check the number of elements
  if (is.null(class.labels)) {class.labels <- as.character(unique.values)} else {
    if (length(class.labels)!=length(unique.values)) {stop('the length of "class.labels" is different from the number of cases in "y"')}}

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3. Check sample distribution per class
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

  # extract layer names
  layer.names <- sapply(1:n.runs, function(j) {paste0('Layer_', as.character(j))})

  # count unique raster values per mask
  sample.count <- lapply(1:n.runs, function(i) {
    erv <-y[which.max(x[[i]])]
    erv <- sapply(unique.values, function(c) {sum(erv==c, na.rm=TRUE)})
    return(list(absolute=erv, relative=erv/sum(erv)))})

  # build final table (absolute count)
  absolute.count <- do.call(cbind, lapply(sample.count, function(i) {i$absolute}))
  colnames(absolute.count) <- layer.names
  absolute.count <- data.frame(code=unique.values, label=class.labels, absolute.count, stringsAsFactors=FALSE)

  # build final table (relative count)
  relative.count <- do.call(cbind, lapply(sample.count, function(i) {data.frame(i$relative)}))
  colnames(relative.count) <- layer.names
  relative.count <- data.frame(code=unique.values, label=class.labels, relative.count, stringsAsFactors=FALSE)

  # data.frame used to a comparative plot
  plot.data <- do.call(rbind, lapply(1:n.runs, function(i) {data.frame(count=sample.count[[i]]$relative, label=class.labels, layer=layer.names[i])}))

  # remove redundant information
  rm(sample.count)

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 4. Build plot
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

  p <- ggplot(plot.data, aes_string(x="layer", y="count", fill="layer")) + geom_bar(stat='identity') +
    ylab('Relative Importance\n') + theme(legend.position='bottom', axis.text.y=element_text(size=8),
                                          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                                          axis.text=element_text(), axis.title=element_text(size=10),
                                          plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm")) +
    facet_wrap(as.formula(paste("~", "label")), nrow=floor(length(unique.values)/4)+1)

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 5. return output
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

  list(absolute.count=absolute.count, relative.count=relative.count, relative.plot=p)

}
