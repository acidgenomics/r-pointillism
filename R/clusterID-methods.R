#' @name clusterID
#' @inherit bioverbs::clusterID
#' @inheritParams basejump::params
#' @examples
#' data(seurat)
#' x <- clusterID(seurat)
#' head(x)
#' table(x)
NULL



#' @rdname clusterID
#' @name clusterID
#' @importFrom bioverbs clusterID
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
