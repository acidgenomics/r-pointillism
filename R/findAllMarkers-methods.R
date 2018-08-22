# Consider adding multi-core support here using BiocParallel



#' Find All Cluster-Specific Markers
#'
#' @rdname findAllMarkers
#' @name findAllMarkers
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [diffExp()].
#'
#' @return `list`
#' - `caller = "edgeR"`: `DGELRT`.
#' - `caller = "DESeq2"`: `DESeqResults`.
#'
#' @examples
#' object <- sce_small
#'
#' # edgeR
#' x <- findAllMarkers(object, caller = "edgeR")
#'
#' # DESeq2
#' x <- findAllMarkers(object, caller = "DESeq2")
setMethod(
    "findAllMarkers",
    signature("SingleCellExperiment"),
    function(object, ...) {
        .assertHasZinbwave(object)
        ident <- ident(object)
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
