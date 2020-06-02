#' rsMove.
#'
#' @name rsMove
#' @description Example datasets of the rsMove package
#' @docType package
#' @import raster sp
NULL

#' Example data of animal movements during a migration.
#'
#'Movement data for one White Stork collected during it migration between Germany and Spain.
#'
#' \itemize{
#'   \item{timestamp}{observation timestamp.}
#'   \item{long}{longitude.}
#'   \item{lat}{latitude.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name longMove
#' @usage data(longMove)
#' @format A SpatialPointsDataFrame
NULL

#' Example data of animal movements during the nesting period.
#'
#'Movement data for one White Stork collected within its nesting site.
#'
#' \itemize{
#'   \item{x}{x coordinate.}
#'   \item{y}{y coordinate.}
#'   \item{date}{observation date.}
#'   \item{time}{observation time.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name shortMove
#' @usage data(shortMove)
#' @format A SpatialPointsDataFrame
NULL
