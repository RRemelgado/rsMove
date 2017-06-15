#' @title sampleMove
#'
#' @description Sampling strategy to identify areas where an animal showed little or no movement based on GPS tracking data.
#' @param x A Vector of x coordinates
#' @param y A Vector of y coordinates
#' @param time vVlid "POSIXlt" or "POSIXct" type vector with the time of observation for each record
#' @param method Selects the appropriate method to estimate the distance between consecutive points. if 'm' it estimates the ecludian distance. If 'deg' it uses the haversine formula
#' @import raster
#' @return Data frame with sample coordinates ('x' and 'y'), total time spent per sample ('time' expressed in minutes) and the total number of observations per sample ('count')
#' @examples \dontrun{
#'
#' # without reference grid
#' moveData <- data('whiteStork')
#' output <- sampleMove(input$x, input$y, input$t, 10, 'm')
#'
#' # with reference grid
#' rasterData <- data('grid')
#' output <- sampleMove(input$x, input$y, input$t, 10, 'm', rasterData)
#' }

#-------------------------------------------------------------------------------------------------------------------------------#

sampleMove <- function(x, y, time, error, method, layer, tUnit) {

  #-----------------------------------------------------------------------------------------------------------------------------#
  # 1. extract samples
  #-----------------------------------------------------------------------------------------------------------------------------#

  # check input parameters
  if (((length(x)+length(y)+length(time)+length(error))/length(x))!=length(x)) {
    stop('error: variable lengths do not match')
  }
  time <- as.numeric(time) # set time as numeric
  if (method!='m' & method!='deg') {stop(paste0('error: method ', method, ' not recognized'))}

  # conver coordinates to radians (if deg)
  if (method=='deg') {
    x0 <- x * pi/180
    y0 <- y * pi/180
  }

  # Identify time segments
  sc <- list()
  sp <- 0
  xs <- list()
  ys <- list()
  ts <- list()
  ss <- list()
  for (r in 2:length(x)) {

    # Estimate distance (harvesine method)
    if (method=='deg') {
      if (sp==0) {
        yDiff <- abs(y0[r]-y0[(r-1)])
        xDiff <- (x0[r]-x0[(r-1)])
      } else {
        yDiff <- abs(y0[r]-y0[sp])
        xDiff <- (x0[r]-x0[sp])
      }
      aCoef <- sin(yDiff/2) * sin(yDiff/2) + cos(y0[r]) * cos(x0[r]) * sin(xDiff/2.) * sin(xDiff/2.)
      cCoef <- 2 * atan2(sqrt(aCoef), sqrt(1.-aCoef))
      lDist <- 6371000 * cCoef
    }

    # estimate distance (ecludian method)
    if (method=='m') {
      if (sp==0) {lDist <- sqrt((x[r]-x[(r-1)])^2 + (y[r]-y[(r-1)])^2)} else {lDist <- sqrt((x[r]-x[sp])^2 + (y[r]-y[sp])^2)}
    }

    # determine if the sample belongs to a new segment
    if (lDist < de & sp==0) {sp <- r-1}
    if (lDist > de & sp>0) {
      sc[[length(sc)+1]] <- c(sp,(r-1))
      sp <- 0
    }

  }

  if (method=='deg') {rm(x0, y0)}

  # summarize time segments
  ns <- length(sc)
  if (ns > 0) {

    xr <- 1:ns
    yr <- 1:ns
    tr <- 1:ns
    nr <- 1:ns
    for (r in 1:length(sc)) {
      loc <- sc[[r]]
      xr[r] <- median(x[loc[1]:loc[2]])
      yr[r] <- median(y[loc[1]:loc[2]])
      tr[r] <- as.numeric(difftime(time[loc[1]], time[loc[2]], units=tUnit))
      nr[r] <- length(loc)
    }

    # update final samples
    xs[[d]] <- xr
    ys[[d]] <- yr
    ts[[d]] <- tr
    ss[[d]] <- nr

    rm(xr, yr, tr, nr, sc)

  }

  #-----------------------------------------------------------------------------------------------------------------------------#
  # 2. build output
  #-----------------------------------------------------------------------------------------------------------------------------#

  # if requested, return unique samples in grid
  if (exists('layer')) {

    xs <- unlist(xs) # x coordinates
    ys <- unlist(ys) # y coordinates
    ts <- unlist(ts) # time spent
    ss <- unlist(ss) # number of samples

    # check layer
    if (!grDevices::is.raster(layer)) {stop('error: "layer" is not a valid raster layer')}

    # derive pixel coordinates
    e <- raster:extent(layer) # layer extent
    pr <- res(layer)[1] # pixel resolution
    nr <- round((e[4]-e[3]) / pr) + 1 # number of rows
    sp <- (round((e[4]-ys[,2])/pr)+1) + nr * round((xs-e[1])/pr) # convert coordinates to pixel positions
    up <- unique(sp) # unique pixel positions

    rm(pr, nr, e)

    # evaluate time spent/number samples
    xr <- 1:length(up)
    yr <- 1:length(up)
    tt <- 1:length(up) # total time
    st <- 1:length(up) # total n° samples
    for (r in 1:length(up)) {
      ind <- which(sp==up[r])
      xr[r] <- mean(xs[ind])
      yr[r] <- mean(ys[ind])
      tt[r] <- sum(ts[ind])
      ts[r] <- sum(ss[ind])
    }

    rm(xs, ys, ts, ss, sp)

    # build/return output table
    os <- data.frame(xr, yr, up, tt, st, stringsAsFactors=F)
    colnames(os) <- c('x', 'y', 'indices', 'time', 'count')
    return(os)

  } else {

    # if no layer is provided return the original sample set
    os <- data.frame(unlist(xs), unlist(ys), unlist(ts), unlist(ss), stringsAsFactors=F)
    os <- c('x', 'y', 'indices', 'time', 'count')
    return(os)

  }

}
