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
#' x <- clusterID(sce_small)
#' table(x)
#' head(x)
NULL



#' @rdname clusterID
#' @export
setMethod(
    "clusterID",
    signature("SingleCellExperiment"),
    function(object) {
        object <- as(object, "SingleCellExperiment")
        x <- colData(object)[["ident"]]
        assert_is_factor(x)
        names(x) <- colnames(object)
        x
    }
)



#' @rdname clusterID
#' @export
setMethod(
    "clusterID",
    signature("seurat"),
    getMethod("clusterID", "SingleCellExperiment")
)
