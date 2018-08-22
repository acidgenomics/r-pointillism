#' Cluster Identity
#'
#' @name ident
#' @family Marker Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `factor`.
#'
#' @examples
#' ident(sce_small)
NULL



#' @rdname ident
#' @export
setMethod(
    "ident",
    signature("SingleCellExperiment"),
    function(object) {
        colData(object)[["ident"]]
    }
)



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    "ident",
    signature("seurat"),
    function(object) {
        slot(object, "ident")
    }
)
