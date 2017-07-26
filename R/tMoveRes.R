#' @title tMoveRes
#'
#' @description Provides historical information on cloud cover.
#' @param xy Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param dpath Output data path for downloaded data.
#' @param p.res Should the output be ploted on screen? Default is TRUE.
#' @import ggplot2 sp rgdal grDevices
#' @importFrom utils download.file
#' @return A \emph{list}.
#' @details {This function makes use of the global cloud cover statistics provided 
#' through EarthEnv (\url{http://www.earthenv.org/cloud}) to extract base statistics 
#' for a shapefile containing GPS tracking data. The function starts by downloading 
#' monthly means of fractional cloud cover. This requires 9.63 Gb of space on the disk. 
#' If the data already exists in \emph{dpath} this step will ignored. If the data 
#' already exists, please assure it has not been modified. Then, the values for each 
#' point are extracted (\emph{$point.data}) and used to derive the min, max, mean and 
#' sd for each month (\emph{$stats}). This statistics are then used to build a plot 
#' with ggplot expressing the mean and maximum cloud cover fraction per month. This 
#' plot will be added to the output (\emph{$plot}). If \emph{p.res} is TRUE, it is 
#' also ploted on screen.}
#' @references \url{https://doi.org/10.1371/journal.pbio.1002415}
#' @seealso \code{\link{sMoveRes}}
#' @examples \dontrun{
#'  
#'  require(raster)
#'  
#'  # read raster data
#'  r <- raster(system.file('extdata', 'tcb_1.tif', package="rsMove"))
#'  
#'  # read movement data
#'  moveData <- read.csv(system.file('extdata', 'konstanz_20130804.csv', package="rsMove"))
#'  moveData <- SpatialPointsDataFrame(moveData[,1:2], moveData, proj4string=crs(r))
#'  
#'  # test function for 5, 10 20 and 30 m
#'  a.res <- tMoveRes(xy=moveData, dpath='.')
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

tMoveRes <- function(xy=xy, dpath=dpath, p.res=T) {
  
#---------------------------------------------------------------------------------------------------------------------#
#  1. check inpur variables
#---------------------------------------------------------------------------------------------------------------------#
  
  if (!class(xy)[1]%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"xy" is not of a valid class')}
  rr <- crs(xy) # reference projection
  if (is.na(rr@projargs)) {stop('"xy" does not have a valid projection')}
  if (!dir.exists(dpath)) {stop('"dpath" not found in file system')}
  if (!is.logical(p.res)) {stop('"p.res" is not a logical argument')}
  
#---------------------------------------------------------------------------------------------------------------------#
#  2. list, check for and download files
#---------------------------------------------------------------------------------------------------------------------#
  
  # target files
  i.files <- c(paste0('http://data.earthenv.org/cloud/MODCF_monthlymean_', sprintf('%02d', 1:12), '.tif'))
  o.files <- file.path(dpath, basename(i.files))
  
  # download missing files
  cc <- file.exists(o.files)
  if (sum(!cc)>0) {
    print('downloading files')
    ind <- which(!cc)
    for (f in 1:length(ind)) {download.file(i.files[ind[f]], o.files[ind[f]], mode="wb")}}

#---------------------------------------------------------------------------------------------------------------------#
#  3. read, crop and scale data
#---------------------------------------------------------------------------------------------------------------------#
  # column names for extracted values
  cn <- c('Jan', 'Feb', 'Mar', 'Apr', 
          'May', 'Jun', 'Jul', 'Aug', 
          'Sep', 'Oct', 'Nov', 'Dec')

  # read and crop image stack
  c.data <- lapply(o.files, function(x){extract(raster(x), xy)*0.01})
  c.data <- do.call(cbind, lapply(c.data, data.frame))
  colnames(c.data) <- cn

#---------------------------------------------------------------------------------------------------------------------#
#  4. derive statistics
#---------------------------------------------------------------------------------------------------------------------#
  
  # estimate base area statistics
  mnv <- apply(c.data, 2, min, na.rm=T)
  mxv <- apply(c.data, 2, max, na.rm=T)
  avv <- apply(c.data, 2, mean, na.rm=T)
  sdv <- apply(c.data, 2, sd, na.rm=T)
  
  # aggregate stats in common data frame (used for output)
  odf <- data.frame(month=cn, min=mnv, max=mxv, mean=avv, sd=sdv, stringsAsFactors=F)
  
  rm(mnv, mxv, avv, sdv)

#---------------------------------------------------------------------------------------------------------------------#
#  5. build plot
#---------------------------------------------------------------------------------------------------------------------#
  
  # make color ramp
  cr = colorRampPalette(c("forestgreen", "khaki2", "darkred"))
  
  # plot
  p <- ggplot(odf, aes(x=factor(month, levels=cn), y=mean, fill=max)) + 
    scale_fill_gradientn(colors=cr(10), breaks=c(0, 50, 100), 
                         limits=c(0,100), name="Max. %\n") + 
    xlab("\nMonth") + ylab("Mean Cloud Cover (%)\n") + ylim(0,100) + 
    geom_bar(width=0.7, stat="identity") + 
    theme(axis.text.x=element_text(size=12), 
          axis.title.x =element_text(size=14), 
          axis.text.y=element_text(size=12),
          axis.title.y =element_text(size=14),
          legend.text=element_text(size=12), 
          legend.title=element_text(size=14))
  
  if (p.res) {p} # plot raster on screen
    
#---------------------------------------------------------------------------------------------------------------------#
#  7. derive output
#---------------------------------------------------------------------------------------------------------------------#

  # return data frame and plot
  return(list(stats=odf, point.data=c.data, plot=p))
  
}