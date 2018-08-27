#' Find Cluster-Specific Marker Genes
#'
#' @note Cluster identity (`ident`) must be defined in `colData()` for this
#'   function to work.
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
#' x <- suppressMessages(
#'     findMarkers(sce_small, caller = "edgeR")
#' )
#' class(x)
#' lapply(x, class)
NULL



#' @rdname findMarkers
#' @export
setMethod(
    "findMarkers",
    signature("SingleCellExperiment"),
    function(object, ...) {
        # Object must contain pre-calculate ZINB weights.
        .assertHasZinbwave(object)
        # Get the cluster identities.
        ident <- clusterID(object)
        assert_is_factor(ident)
        clusters <- levels(ident)
        stopifnot(length(clusters) >= 2L)
        message(paste(length(clusters), "clusters detected"))
        # Loop across the clusters and calculate gene enrichment relative to
        # all of the other clusters combined.
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
