#' @title spaceDir
#'
#' @description Analysis of environmental change in space along a movement track.
#' @param x Object of class \emph{SpatialPoints} or \emph{SpatialPointsDataFrame}.
#' @param obs.time Object of class \emph{Date}, \emph{POSIXlt} or \emph{POSIXct} with \emph{x} observation dates.
#' @param y Object of class \emph{RasterLayer}.
#' @param sample.direction One of \emph{forward}, \emph{backward} or \emph{both}. Default is \emph{both}.
#' @param data.type One of 'cont' or 'cat'. Defines which type of variable is in use.
#' @param distance.method One of 'm' or 'deg' specifying the projection unit. Default is 'm'.
#' @param buffer.size Spatial buffer size expressed in the map units.
#' @param stat.fun Output statistical metric.
#' @param min.count Minimum number of pixels required by \emph{stat.fun}. Default is 2.
#' @importFrom raster extract
#' @importFrom stats lm
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 aes_string geom_line theme scale_color_gradientn xlab ylab fortify geom_path element_text element_blank geom_histogram xlim
#' @seealso \code{\link{timeDir}} \code{\link{dataQuery}} \code{\link{imgInt}}
#' @return A \emph{list} containing shapefiles with information on environmental change and travel distance/time and a plot of the results.
#' @details {This function quantifies environmnetal changes in space along a movement track. For each set of consecutive points, the
#' function will derive coordinates for all unique pixels that overlap with - and are beteween - each point in a spatial moving widow
#' with a shape defined by \emph{sample.direction}. The user can define \emph{buffer.size} to consider all pixels within a predefined
#' spatial distance of the initially selected pixels. Once all pixels are selected, the function will extract the corresponding values
#' in \emph{y} summarize them using a pre-defined statistical function. If \emph{data.type} is \emph{cont}, a statistical function can
#' be provided through \emph{stat.fun}. However, if \emph{data.type} is \emph{cat}, the function will report on the dominant class and
#' on the shannon index for each window. If the number of pixels within the a window is lower than \emph{min.count} the function will
#' return \emph{NA}. On top of this, \emph{spaceDir} will also report on the linear distance traveled between endpoints (in meters) and
#' the travel time (in minutes). The output of the function is a list consisting of:
#' \itemize{
#'  \item{\emph{endpoints} - Point shapefile with endpoints of each spatial segment. Reports on a given statistical metric, traveled distance, travel time and the mean timestamp.}
#'  \item{\emph{segments} - Line shapefile with spatial segments. Reports on the same information as \emph{endpoints}.
#'  \item{\emph{plot} - Plotting of \emph{segments} where each segment is colored according to the extracted raster value.}}}}
#'
#' @examples {
#'
#'  require(raster)
#'
#'  # read raster data
#'  r <- raster(system.file('extdata', '2013-07-16_ndvi.tif', package="rsMove"))
#'
#'  # read movement data
#'  data(shortMove)
#'
#'  # observation time
#'  obs.time <- strptime(paste0(shortMove@data$date, ' ',shortMove@data$time),
#'  format="%Y/%m/%d %H:%M:%S")
#'
#'  # perform directional sampling
#'  of <- function(i) {lm(i~c(1:length(i)))$coefficients[2]}
#'  s.sample <- spaceDir(shortMove, r, "backward", "cont", obs.time=obs.time, stat.fun=of)
#'
#' }
#' @export

#-------------------------------------------------------------------------------------------------------------------------------#

spaceDir <- function(x, y, sample.direction, data.type, obs.time=NULL, distance.method='m', buffer.size=NULL, stat.fun=NULL, min.count=2) {

#-------------------------------------------------------------------------------------------------------------------------------#
# 1. check variables
#-------------------------------------------------------------------------------------------------------------------------------#

  # samples
  if (!exists('x')) {stop('"x" is missing')}
  if (!class(x)%in%c('SpatialPoints', 'SpatialPointsDataFrame')) {stop('"x" is not of a valid class')}

  # sample dates
  if (!is.null(obs.time)) {
    if (!class(obs.time)[1]%in%c('Date', 'POSIXct', 'POSIXlt')) {stop('"obs.time" is nof of a valid class')}
    if (length(obs.time)!=length(x)) {stop('errorr: "x" and "obs.time" have different lengths')}}

  # raster
  if (!exists('y')) {stop('"y" is missing')}
  if (!class(y)[1]=='RasterLayer') {stop('"y" is not of a valid class')}
  if (crs(x)@projargs!=crs(y)@projargs) {stop('"x" and "y" have different projections')}

  # query direction
  if (!is.null(sample.direction)) {
    if (length(sample.direction)>1) {stop('"sample.direction" has too many entries')}
    if (!sample.direction%in%c('backward', 'forward', 'both')) {stop('"sample.direction" is not a valid entry')}
  } else {sample.direction <- 'both'}

  # variable data.type
  if (is.null(data.type)) {stop('"data.type" is missing')} else {
    if (!data.type%in%c('cont', 'cat')) {stop('"data.type" is not a recognized keyword')}}

  # check/define input metrics
  if (is.null(stat.fun)) {stat.fun <- function(i) {lm(i~c(1:length(i)))$coefficients[2]}} else {
    if(!is.function(stat.fun)) {stop('"stat.fun" is not a valid function')}}

#-------------------------------------------------------------------------------------------------------------------------------#
# 2. select pixels between consecutive points
#-------------------------------------------------------------------------------------------------------------------------------#

  # base raster info
  rProj <- crs(y)
  pxr <- res(y)

  # backward sampling
  if (sample.direction=='backward') {
    f1 <- function(i) {
      si <- i-1
      ei <- i
      d0 <- sqrt((x@coords[si,1]-x@coords[ei,1])^2 + (x@coords[si,2]-x@coords[ei,2])^2)
      if (!is.null(obs.time)) {t0 <- difftime(obs.time[ei], obs.time[si], units='mins')} else {t0 <- NA}
      x0 <- x@coords[si:ei,1]
      y0 <- x@coords[si:ei,2]
      if((d0 > (pxr[1]*2))) {
        dx <- x0[1]-x0[2]
        dy <- y0[1]-y0[2]
        if (abs(dx)>abs(dy)) {
          m <- lm(y0~x0)$coefficients
          if (dx>0) {cm<- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dx)/pxr[1]))*pxr[1])}
          x0 <- (x0[1]+cm)
          y0 <- m[1] + x0 * m[2]
        } else {
          m <- lm(x0~y0)$coefficients
          if (dy>0) {cm<- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm<- ((0:round(abs(dy)/pxr[1]))*pxr[1])}
          y0 <- (y0[1]+cm)
          x0 <- m[1] + y0 * m[2]}}
      if (distance.method=='deg') {
        x00 <- x@coords[si:ei,1]*pi/180
        y00 <- x@coords[si:ei,2]*pi/180
        sf <- function(o) {
          xDiff <- abs(x00[(o-1)]-x00[o])
          yDiff <- abs(y00[(o-1)]-y00[o])
          aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y00[o]) * cos(x00[o]) * sin(xDiff/2.) * sin(xDiff/2.)
          cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
          return(6371000 * cCoef)}
        d0 <- sum(sapply(2:length(x00), sf))}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), i)))}
    op <- lapply(2:length(x), f1)}

  # forward sampling
  if (sample.direction=='forward') {
    f1 <- function(i) {
      si <- i
      ei <- i+1
      d0 <- sqrt((x@coords[si,1]-x@coords[ei,1])^2 + (x@coords[si,2]-x@coords[ei,2])^2)
      if (!is.null(obs.time)) {t0 <- difftime(obs.time[ei], obs.time[si], units='mins')} else {t0 <- NA}
      x0 <- x@coords[si:ei,1]
      y0 <- x@coords[si:ei,2]
      if((d0 > (pxr[1]*2))) {
        dx <- x0[1]-x0[2]
        dy <- y0[1]-y0[2]
        if (abs(dx)>abs(dy)) {
          m <- lm(y0~x0)$coefficients
          if (dx>0) {cm<- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm<-((0:round(abs(dx)/pxr[1]))*pxr[1])}
          x0 <- (x0[1]+cm)
          y0 <- m[1] + x0 * m[2]
        } else {
          m <- lm(x0~y0)$coefficients
          if (dy>0) {cm<- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm<-((0:round(abs(dy)/pxr[1]))*pxr[1])}
          y0 <- (y0[1]+cm)
          x0 <- m[1] + y0 * m[2]}}
      if (distance.method=='deg') {
        x00 <- x@coords[si:ei,1]*pi/180
        y00 <- x@coords[si:ei,2]*pi/180
        sf <- function(o) {
          xDiff <- abs(x00[(o-1)]-x00[o])
          yDiff <- abs(y00[(o-1)]-y00[o])
          aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y00[o]) * cos(x00[o]) * sin(xDiff/2.) * sin(xDiff/2.)
          cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
          return(6371000 * cCoef)}
        d0 <- sum(sapply(2:length(x00), sf))}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), i)))}
    op <- lapply(1:(length(x)-1), f1)}

  # backward-forward sampling
  if (sample.direction=='both') {
    f1 <- function(i) {
      si <- i-1
      ei <- i+1
      d0 <- sqrt((x@coords[si,1]-x@coords[ei,1])^2 + (x@coords[si,2]-x@coords[ei,2])^2)
      if (!is.null(obs.time)) {t0 <- difftime(obs.time[ei], obs.time[si], units='mins')} else {t0 <- NA}
      x0 <- x@coords[si:ei,1]
      y0 <- x@coords[si:ei,2]
      if((d0 > (pxr[1]*2))) {
        x00 <- vector('list', 2)
        y00 <- vector('list', 2)
        for (o in 2:3) {
          dx <- x0[(o-1)]-x0[o]
          dy <- y0[(o-1)]-y0[o]
          if (abs(dx)>abs(dy)) {
            m <- lm(y0~x0)$coefficients
            if (dx>0) {cm <- -((0:round(abs(dx)/pxr[1]))*pxr[1])} else {cm <- ((0:round(abs(dx)/pxr[1]))*pxr[1])}
            tx <- (x0[(o-1)]+cm)
            ty <- m[1] + x0[(o-1):o] * m[2]
          } else {
            m <- lm(x0~y0)$coefficients
            if (dy>0) {cm <- -((0:round(abs(dy)/pxr[1]))*pxr[1])} else {cm <- ((0:round(abs(dy)/pxr[1]))*pxr[1])}
            y0 <- (y0[(o-1)]+cm)
            x0 <- m[1] + y0[(o-1):o] * m[2]}
            if (o==3) {
              x00[(o-1)] <- tx[2:length(tx)]
              y00[(o-1)] <- ty[2:length(ty)]
            } else {
              x00[(o-1)] <- tx
              y00[(o-1)] <- ty}}
        x0 <- unlist(x00)
        y0 <- unlist(y00)}
      if (distance.method=='deg') {
        x00 <- x@coords[si:ei,1]*pi/180
        y00 <- x@coords[si:ei,2]*pi/180
        sf <- function(o) {
          xDiff <- abs(x00[(o-1)]-x00[o])
          yDiff <- abs(y00[(o-1)]-y00[o])
          aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y00[o]) * cos(x00[o]) * sin(xDiff/2.) * sin(xDiff/2.)
          cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
          return(6371000 * cCoef)}
        d0 <- sum(sapply(1:length(x00), sf))}
      return(list(x=x0, y=y0, d=d0, t=t0, p=replicate(length(x0), i)))}
    op <- lapply(2:(length(x)-1), f1)}

  # retrieve x, y and p
  xc <- unlist(lapply(op, function(i){i$x}))
  yc <- unlist(lapply(op, function(i){i$y}))
  pd <- unlist(lapply(op, function(i){i$d}))
  td <- unlist(lapply(op, function(i){i$t}))
  us <- unlist(lapply(op, function(i){i$p}))

  rm(op)

#-------------------------------------------------------------------------------------------------------------------------------#
# 3. query samples
#-------------------------------------------------------------------------------------------------------------------------------#

  # I. apply spatial buffer (if prompted)
  # II. retrieve environmental variables
  if (!is.null(buffer.size)) {

    # dilate samples and update sample indices
    tmp <- lapply(unique(us), function(i) {
      ind <- which(us==i)
      ind <- do.call(rbind, lapply(ind, function(y) {
        r0 <- raster(extent((xc[y]-buffer.size), (xc[y]+buffer.size),
                            (yc[y]-buffer.size), (yc[y]+buffer.size)),
                     res=pxr[1], crs=rProj)
        return(xyFromCell(r0, 1:ncell(r0)))}))
      ind <- ind[!duplicated(cellFromXY(y, ind)),]
      return(list(c=ind, s=replicate(nrow(ind), i)))})
    us1 <- unlist(lapply(tmp, function(i) {i$s}))
    tmp <- do.call(rbind, lapply(tmp, function(i) {i$c}))

    # retrieve environmental data
    edata <- extract(y, tmp)

    rm(tmp)

  } else {
    us1 = us
    edata <- extract(y, cbind(xc,yc))
  }

#-------------------------------------------------------------------------------------------------------------------------------#
# 4. analyze samples
#-------------------------------------------------------------------------------------------------------------------------------#

  if (data.type=='cont') {

    # query function
    f2 <- function(i) {
      ind <- which(us1==i)
      u <- which(!is.na(edata[ind]))
      if (sum(u) >= min.count) {return(as.numeric(stat.fun(as.numeric(edata[ind[u]]))))} else {return(NA)}}

    # apply user provided functon
    ov <- data.frame(stat=sapply(unique(us1), f2))

  }

  if (data.type=='cat') {

    # unique classes
    uc <- unique(y)

    # function to sumarize class composition
    f2 <- function(i) {
      ind <- which(us1==i)
      return(sapply(uc, function(y) {sum(edata[ind]==y, na.rm=T)}))}
    ov <- do.call(rbind, lapply(unique(us1), f2))

    # derive additional statistics
    dc <- apply(ov, 1, function(i) {(uc[which(i==max(i))])[1]})
    ind <- !is.na(uc) # avoid NA values are used in calculation
    si <- apply(ov[,ind], 1, function(i) {
      ind <- which(i>0)
      i[ind] <- (i[ind]/sum(i[ind]))*log(i[ind]/sum(i[ind]))
      return(-sum(i))})

    # update output table
    ov <- as.data.frame(cbind(ov, dc, si))
    colnames(ov) <- c(as.character(uc), 'main', 'shannon')

  }

#-------------------------------------------------------------------------------------------------------------------------------#
# 5. build shapefiles
#-------------------------------------------------------------------------------------------------------------------------------#

  # build dara table
  us0 <- unique(us)

  # subset variables (depends on sampling strategy)
  df <- data.frame(x=x@coords[us0,1], y=x@coords[us0,2], timestamp=obs.time[us0], travel.distance=pd, travel.time=td, sid=us0)
  df <- cbind(df, ov)

  # build segment endpoint shapefile
  p.shp <- SpatialPointsDataFrame(df[,1:2], df, proj4string=rProj)

  # build segment path shapefile
  f <- function(i) {
    loc <- which(us==us0[i])
    return(Lines(list(Line(cbind(xc[loc],yc[loc]))), i))}
  l.shp = SpatialLinesDataFrame(SpatialLines(lapply(1:length(us0), f), proj4string=rProj), df)

  #-------------------------------------------------------------------------------------------------------------------------------#
  # 5. build plot
  #-------------------------------------------------------------------------------------------------------------------------------#

  df <- data.frame(df, id=1:nrow(df))
  fl <- fortify(l.shp)
  fl <- merge(fl, df, by.x="id", by.y="id")

  cr <- colorRampPalette(c("dodgerblue3", "khaki2", "forestgreen"))

  if (data.type=="cont") {

    p <- ggplot(fl, aes_string(x="long", y="lat", color="stat", group="group")) + theme_bw() +
      geom_path(size=2) + xlab("X") + ylab("Y") +
      theme(legend.text=element_text(size=10),
            axis.text=element_text(size=10),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank()) +
      scale_color_gradientn(name="Stat\n", colors=cr(10),
                            limits=c(min(df$stat), max((df$stat))))

  }

  if (data.type=="cat") {

    mv <- round(max(df$shannon, na.rm=T))
    nc <- nchar(as.character(mv))
    m <- as.numeric(paste0(1, paste0(replicate((nc-1), '0'), collapse='')))
    mv <- mv / m
    vl <- round(mv)
    if (mv > vl) {vl <- (vl+0.2)*m} else {tb <- vl*m}

    p <- ggplot(fl, aes_string(x="long", y="lat", color="shannon")) + theme_bw() +
      geom_line(size=2) + xlab("X") + ylab("Y") +
      theme(legend.text=element_text(size=10),
            axis.text=element_text(size=10),
            panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank()) +
      scale_color_gradientn(name="Shannon\n", colors=cr(10), breaks=c(0.0, (vl/2), vl), limits=c(0.0, vl))

  }

  #-------------------------------------------------------------------------------------------------------------------------------#
  # 7. derive output
  #-------------------------------------------------------------------------------------------------------------------------------#

  # output
  return(list(endpoints=p.shp, segments=l.shp, plot=p))

}
