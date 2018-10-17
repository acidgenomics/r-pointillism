# TODO Consider adding `progress` option.



# @examples
# data(seurat_small)
# x <- findMarkers(seurat_small, caller = "edgeR")
# class(x)
# lapply(x, class)



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
NULL



.findMarkers.SCE <-  # nolint
    function(object, ...) {
        # Object must contain pre-calculate ZINB weights.
        .assertHasZinbwave(object)

        # Get the cluster identities.
        ident <- clusterID(object)
        assert_is_factor(ident)
        assert_has_names(ident)
        clusters <- levels(ident)
        stopifnot(length(clusters) >= 2L)
        message(paste(length(clusters), "clusters detected"))

        # Loop across the clusters and calculate gene enrichment relative to
        # all of the other clusters combined.
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                message(paste("Cluster", cluster, "===="))
                # Numerator: cells in the current cluster.
                numerator <- ident[which(ident == cluster)]
                stopifnot(all(numerator == cluster))
                numerator <- sort(names(numerator))
                assert_is_non_empty(numerator)
                # Denominator: cells in all other clusters.
                denominator <- sort(setdiff(colnames(object), numerator))
                assert_is_non_empty(denominator)
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



#' @rdname findMarkers
#' @export
setMethod(
    f = "findMarkers",
    signature = signature("SingleCellExperiment"),
    definition = .findMarkers.SCE
)



#' @rdname findMarkers
#' @export
setMethod(
    f = "findMarkers",
    signature = signature("seurat"),
    definition = getMethod("findMarkers", "SingleCellExperiment")
)
