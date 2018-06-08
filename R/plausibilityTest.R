#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#' @title plausibilityTest
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#' @description Quantifies and plots the distribution of pixels within a mask over a reference categorical raster object. 
#' @param sample.mask Object of class \emph{raster} and \emph{RasterStack} with the sample mask(s).
#' @param reference.map Reference \emph{raster} (e.g. Land cover map).
#' @param class.labels Labels of classes in \emph{reference.map} provided as a character vector.
#' @param sample.colors Hex color codes for each layer in \emph{sample.mask} provided as a character vector.
#' @return A \emph{list} with statistical information on the distribution of samples per class and a comparative plot.
#' @details {The function counts the number of non-NA pixels in \emph{sample.mask} within each class of \emph{reference.map}. 
#' Then, the sample count is normalized by its largest value. The output of the function is a list consisting of:
#'  \itemize{
#'  \item{\emph{absolute.count} - Absolute sample count for each layer of \emph{sample.mask} within each class of \emph{reference.map.}
#'  \item{\emph{relative.count} - Relative sample count for each layer of \emph{sample.mask} within each class of \emph{reference.map.}
#'  \item{\emph{relative.plot} - Plot comparing the relative sample count of the layers in \emph{sample.mask} within each class of \emph{reference.map}.}}}
#' @importFrom raster crs extract nlayers
#' @importFrom ggplot2 ggplot aes_string geom_bar theme_bw ylim theme xlabl ylab scale_fill_manual facet_wrap
#' @importFrom grDevices rainbow

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------------------------------------------#

plausibilityTest <- function(sample.mask=sample.mask, reference.map=reference.map, class.labels=NULL, sample.colors=NULL) {
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 1. Check input variables
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  # check sample layer
  if (!class(sample.mask)[1] %in% c("RasterLayer", "RasterStack")) {stop('"sample.mask" is not a raster object')}
  
  # check reference layer
  if (class(reference.map)[1]!="RasterLayer") {stop('"reference.map" is not a valid raster object')}
  if (crs(sample.mask)@projargs!=crs(reference.map)@projargs) {stop('"sample.mask" and "reference.map" have different projections')}
  if (sum(dim(sample.mask)[1:2]-dim(reference.map)[1:2])!=0) {stop('"sample.mask" & "reference.map" have different dimensions')}
  n.runs <- nlayers(sample.mask)
  
  # check auxiliary information
  if (!is.null(class.labels)) {if (!is.character(class.labels)) {'"class.labels" is not a "character" vector'}}
  if (!is.null(sample.colors)) {
    if (!is.character(sample.colors)) {'"sample.colors" is not a "character" vector'}
    if (length(sample.colors)!=n.runs) {'"sample.mask" and "sample.colors" have different lenghts'}
  } else {sample.colors <- rainbow(n.runs, start=0.1, end=0.9)}
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 2. Extract unique cases in "reference.map" and check "class.labels" and "class.colors"
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  # unique raster values
  unique.values <- sort(unique(reference.map))
  unique.values <- unique.values[!is.na(unique.values)]
  
  # define class labels and check the number of elements
  if (is.null(class.labels)) {class.labels <- as.character(unique.values)} else {
    if (length(class.labels)!=length(unique.values)) {stop('the length of "class.labels" is different from the number of cases in "reference.map"')}}

#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 3. Check sample distribution per class
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  # extract layer names
  layer.names <- names(sample.mask)
  
  # count unique raster values per mask
  sample.count <- lapply(1:n.runs, function(i) {
    erv <-reference.map[which.max(sample.mask[[i]])]
    erv <- sapply(unique.values, function(c) {sum(erv==c, na.rm=TRUE)})
    return(list(absolute=erv, relative=erv/max(erv)))})
  
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
    facet_wrap(as.formula(paste("~", "label")), nrow=floor(length(unique.values)/4)+1) + scale_fill_manual(values=sample.colors)
  
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
# 5. return output
#----------------------------------------------------------------------------------------------------------------------------------------------------------#
  
  list(absolute.count=absolute.count, relative.count=relative.count, relative.plot=p)
  
}