#' @title backSampling
#'
#' @description updates a given set of samples derived with labelSample(). It performs a selection of background samples based on a PCA analysis.
#' @param input matrix provided by labelSample()
#' @param features
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

# i1 = which(label!=0)
# i0 = which(label==0)
#
# # select absences
# pca = prcomp(x, scale=T, center=T)
# pca = pca$x[,pcf(pca)]
# ai = vector('list', ncol(pca))
# uv = unique(g[i1])
# for (p in 1:ncol(pca)) {
#   usr = vector('list', length(uv))
#   for (z in 1:length(uv)) {
#     ri = which(g==uv[z])
#     s1 = median(pca[ri,p])
#     s2 = median(abs(pca[ri,p]-s1))
#     usr[[z]] = i0[which(abs(pca[i0,p]-s1) > (s2*sdm))]
#   }
#   usr = unlist(usr)
#   ui = unique(usr)
#   count = vector('numeric', length(ui))
#   for (z in 1:length(ui)) {count[z] = length(which(usr==ui[z]))}
#   ai[[p]] = ui[which(count==length(uv))]
# }
# ai = unique(unlist(ai))
#
# rm(pca, s1, s2, usr, count)
