#' @title satTime
#'
#' @description Finds available satellite dates for tracking data.
#' @param o.time Object of class \emph{Date}.
#' @param t.var Target variable.
#' @param p.res Logical. Should the output be ploted?
#' @import ggplot2
#' @return One or multiple plots.
#' @details {This function compares a set of input dates (\emph{o.time}) to the possible dates for 
#' which satellite data can be downloaded. This analysis is performed for the set of variables 
#' provided through the function proSat(). The function provides a report with the closest dates 
#' available for each of the samples and a faceted plot for each unique year showing the observed 
#' dates, the dates covered by the satellite data and the dates covered by both.}
#' @seealso \code{\link{moveCloud}} \code{\link{proSat}}
#' @examples {
#'  
#'  # return list of variables
#'  var.ls <- satTime(o.time=as.Date("2013-08-04"), t.var="ndvi")
#'  
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

satTime <- function(o.time=o.time, t.var=t.var, p.res=T) {
  
#-----------------------------------------------------------------------------------------------------------------#
# 1. check input variables  
#-----------------------------------------------------------------------------------------------------------------#
  
  if (class(o.time)[1] != 'Date') {stop('"o.time" is not of class "Date"')}
  if (!is.character(t.var)) {stop('"t.var" is not a "character"')}
  if (length(t.var)>1) {stop('"t.var" has more than 1 element')}
  
  # read variable list
  var.ls <- system.file('extdata', 'sat-variables.csv', package="rsMove")
  var.ls <- var.ls <- read.csv(var.ls, stringsAsFactors=F)
  
#-----------------------------------------------------------------------------------------------------------------#
# 2. determine potential temporal distribution  
#-----------------------------------------------------------------------------------------------------------------#
  
  # temporal resolution of target variable
  t.res <- var.ls$temporal.resolution..days.[var.ls$code==t.var]
  
  # determine day of year for provided dates
  o.yrs <- sapply(as.character(o.time), function(x) {strsplit(x, '-')[[1]][1]})
  o.doa <- (o.time-as.Date(paste0(o.yrs, '-01-01')))+1
  s.doa <- (seq(1, 365, t.res))
  
#-----------------------------------------------------------------------------------------------------------------#
# 3. compare provided and target dates
#-----------------------------------------------------------------------------------------------------------------#
  
  tmp <- lapply(1:length(o.time), function(x) {
    diff <- s.doa - o.doa[x]
    i0 <- which(diff==0)
    if (length(i0>0)) {
      return(list(bd=(as.Date(paste0(o.yrs[x],'-01-01'))+(s.doa[i0]-1)), 
                  ad=(as.Date(paste0(o.yrs[x],'-01-01'))+(s.doa[i0]-1))))
    } else {
      ib <- which(diff < 0)
      if(length(ib)>0) {bd<-s.doa[which(diff[ib]==max(diff[ib]))]} else {bd<-NA}
      ia <- which(diff > 0)
      if(length(ia)>0) {ad<-s.doa[which(diff[ia]==max(diff[ia]))]} else {ad<-NA}
      return(list(bd=(as.Date(paste0(o.yrs[x],'-01-01'))+(bd-1)), 
                  ad=(as.Date(paste0(o.yrs[x],'-01-01'))+(ad-1))))}})
  
  # output data frame
  df <- data.frame(original=o.time, 
                   best.before=do.call('c', lapply(tmp, function(x) {x$bd})), 
             best.after=do.call('c', lapply(tmp, function(x) {x$ad})))
  
#-----------------------------------------------------------------------------------------------------------------#
# 4. build plot per year
#-----------------------------------------------------------------------------------------------------------------#
  
  # build data frame used for plotting
  uy <- sort(unique(yrs))
  gg <- do.call(rbind, lapply(uy, function(x) {
    ind <- which(o.yrs==x)
    od <- unique(sort(c(o.doa[ind], s.doa)))
    yc <- vector('character', length(od))
    yc[which(od%in%o.doa[ind] & od%in%s.doa)] <- 'Satellite/Observed' # obs. dates with sat. data
    yc[which(!od%in%o.doa[ind] & od%in%s.doa)] <- 'Satellite' # sat. dates not covering obs.
    yc[which(od%in%o.doa[ind] & !od%in%s.doa)] <- 'Observed' # obs. dates not covering sat.
    yc[yc%in%o.doa[ind]] <- yc[yd%in%o.doa[ind]]
    return(data.frame(class=yc, doy=od, code=1, year=x, stringsAsFactors=F))}))
  
  # define color/class scheme
  uc <- unique(gg$class)
  cn <- c("Observed", "Satellite", "Satellite/Observed")
  cc <- c("firebrick3", "dodgerblue4", "green3")
  ind <- which(cn%in%uc)
  cn <- cn[ind]
  cc <- cc[ind]
  
  # build plot
  p <- ggplot(gg, aes(x=doy, y=code, fill=factor(class, levels=cn))) + 
    theme_bw() + geom_bar(stat="identity", width=1) + ylab("") + xlab("Day of Year") + 
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
          axis.text.x=element_text(size=12, vjust=1), 
          legend.title=element_blank(), legend.position="bottom", 
          legend.text=element_text(size=12)) + 
    scale_x_continuous(breaks=seq(0, 370, 50), lim=c(1, 366), expand=c(0, 0)) + 
    scale_y_continuous(expand=c(0, 0)) + 
    scale_fill_manual(values=cc, labels=cn) + 
    facet_wrap(~factor(year, levels=uy), nrow=length(uy))
  
  if(p.res) {p}
  
  return(list(report=df, plot=p))
  
  
}