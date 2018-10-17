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
#' data(seurat_small)
#' x <- clusterID(seurat_small)
#' head(x)
#' table(x)
NULL



.clusterID.SCE <-  # nolint
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
    definition = .clusterID.SCE
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
