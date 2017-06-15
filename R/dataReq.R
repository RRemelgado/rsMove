#' @title dataReq
#'
#' @description
#' @param x
#' @export
#' @import raster
#' @import rgeos
#' @return
#' @export
#' @usage
#' @keywords
#' @seealso
#' @aliases
#' @examples \dontrun{
#'
#' }

# lshp = 'E:/01_DATA/LANDSAT/infos/WRS_Grid/wrs2_descending.shp'
# mshp = 'E:/01_DATA/MODIS/infos/modis_sinusoidal_grid_world.shp'
# sbbx = 'E:/01_DATA/AM_data/White-Stork/infos/macro-regions_composition.txt'
#
# library(raster)
#
# # read data on tiles and stopovers
# L_tiles = shapefile(lshp)
# M_tiles = shapefile(mshp)
# binfo = read.table(sbbx, header=T, sep=';', stringsAsFactors=F)
#
# binfo = binfo[which(binfo$population=='Germany'),]
#
# # determine tiles that overlap with each stopover
# nr = nrow(binfo)
# L_PR = vector('character', nr)
# M_PR = vector('character', nr)
# for (r in 1:nr) {
#   e = extent(binfo$xMin[r], binfo$xMax[r], binfo$yMin[r], binfo$yMax[r])
#   tmp = crop(L_tiles, e)@data
#   L_PR[r] = paste(tmp$WRSPR, collapse=",")
#   tmp = raster(e)
#   crs(tmp) = crs(L_tiles)
#   e = projectExtent(tmp, crs(M_tiles))
#   tmp = crop(M_tiles, e)@data
#   M_PR[r] = paste(paste0('h', tmp$h, 'v', tmp$v), collapse=",")
# }
#
# lpr = unique(unlist(strsplit(L_PR, ',')))
# mpr = unique(unlist(strsplit(M_PR, ',')))
