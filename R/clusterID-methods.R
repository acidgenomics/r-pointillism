#' @name clusterID
#' @inherit bioverbs::clusterID
#' @note Updated 2019-07-31.
#'
#' @inheritParams basejump::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(seurat)
#' x <- clusterID(seurat)
#' head(x)
#' table(x)
NULL



#' @rdname clusterID
#' @name clusterID
#' @importFrom bioverbs clusterID
#' @usage clusterID(object, ...)
#' @export
NULL



clusterID.SingleCellExperiment <-  # nolint
    function(object) {
        object <- as(object, "SingleCellExperiment")
        x <- colData(object)[["ident"]]
        assert(is.factor(x))
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



clusterID.Seurat <-  # nolint
    clusterID.SingleCellExperiment



#' @rdname clusterID
#' @export
setMethod(
    f = "clusterID",
    signature = signature("Seurat"),
    definition = clusterID.Seurat
)
