#' Cluster Identity
#'
#' @name clusterID
#'
#' @inheritParams general
#'
#' @return `factor`.
#'
#' @examples
#' data(seurat_small)
#' x <- clusterID(seurat_small)
#' head(x)
#' table(x)
NULL



clusterID.SingleCellExperiment <-  # nolint
    function(object) {
        object <- as(object, "SingleCellExperiment")
        x <- colData(object)[["ident"]]
        assert_is_factor(x)
        names(x) <- colnames(object)
        x
    }



#' @rdname clusterID
#' @export
setMethod(
    f = "clusterID",
    signature = signature("SingleCellExperiment"),
    definition = clusterID.SingleCellExperiment
)



#' @rdname clusterID
#' @export
setMethod(
    f = "clusterID",
    signature = signature("seurat"),
    definition = getMethod(
        f = "clusterID",
        signature = signature("SingleCellExperiment")
    )
)
