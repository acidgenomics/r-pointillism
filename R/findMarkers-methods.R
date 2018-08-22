#' Find Cluster-Specific Marker Genes
#'
#' @name findMarkers
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [diffExp()].
#'
#' @return `list` containing:
#' - `caller = "edgeR"`: `DGELRT`.
#' - `caller = "DESeq2"`: `DESeqResults`.
#'
#' @examples
#' object <- sce_small
#'
#' # edgeR
#' x <- findMarkers(object, caller = "edgeR")
#'
#' # DESeq2
#' x <- findMarkers(object, caller = "DESeq2")
NULL



#' @rdname findMarkers
#' @export
setMethod(
    "findMarkers",
    signature("SingleCellExperiment"),
    function(object, ...) {
        .assertHasZinbwave(object)
        ident <- clusterID(object)
        assert_is_factor(ident)
        clusters <- levels(ident)
        stopifnot(length(clusters) >= 2L)
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                message(paste("Cluster", cluster, "===="))
                numerator <- colnames(object)[which(ident == cluster)]
                denominator <- colnames(object)[which(ident != cluster)]
                diffExp(
                    object = object,
                    numerator = numerator,
                    denominator = denominator,
                    ...
                )
            }
        )
        names(list) <- clusters
        list
    }
)
