#' Cluster Identity
#'
#' @name clusterID
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `factor`.
#'
#' @examples
#' clusterID(sce_small)
NULL



#' @rdname clusterID
#' @export
setMethod(
    "clusterID",
    signature("SingleCellExperiment"),
    function(object) {
        colData(object)[["ident"]]
    }
)



#' @rdname clusterID
#' @export
setMethod(
    "clusterID",
    signature("seurat"),
    function(object) {
        slot(object, "ident")
    }
)
