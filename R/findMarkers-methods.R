## FIXME Consider moving this to AcidSingleCell.



#' @name findMarkers
#' @inherit AcidGenerics::findMarkers
#' @note Updated 2021-09-13.
#'
#' @inheritParams AcidRoxygen::params
#' @param clusters `character` or `NULL`.
#'   Cluster identifiers.
#'   Must correspond to values in [clusters()].
#'   Note that Seurat uses zero-indexed IDs by default (e.g. 0, 1, 2, ...).
#'   If left `NULL` (default), all clusters will be analyzed.
#'   Looping across clusters manually here can avoid memory issues on laptops
#'   and other machines with limited amounts of RAM.
#' @param ... Passthrough arguments to [diffExp()].
#'
#' @return `list` containing:
#' - `caller = "edgeR"`: `DGELRT`.
#' - `caller = "DESeq2"`: `DESeqResults`.
#'
#' @examples
#' data(Seurat, package = "AcidTest")
#'
#' ## Seurat ====
#' object <- Seurat
#' x <- findMarkers(object, caller = "edgeR")
#' class(x)
#' lapply(x, class)
NULL



## Updated 2020-09-03.
`findMarkers,SCE` <-  # nolint
    function(
        object,
        clusters = NULL,
        ...
    ) {
        assert(isCharacter(clusters, nullOK = TRUE))
        h1(sprintf("{.fun %s}", "findMarkers"))
        object <- as(object, "SingleCellExperiment")
        ## Get the cluster mappings. Following the Seurat nomenclature here of
        ## using "ident" to denote the cluster identifier mapping factor.
        ident <- clusters(object)
        assert(is.factor(ident), hasNames(ident))
        assert(
            length(levels(ident)) >= 2L,
            msg = "Object does not contain 2 or more clusters."
        )
        if (is.null(clusters)) {
            clusters <- levels(ident)
        }
        alertInfo(sprintf(
            "%d %s detected.",
            length(clusters),
            ngettext(
                n = length(clusters),
                msg1 = "cluster",
                msg2 = "clusters"
            )
        ))
        ## Loop across the clusters and calculate gene enrichment relative to
        ## all of the other clusters combined.
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                h2(sprintf("Cluster %s", as.character(cluster)))
                ## Numerator: cells in the current cluster.
                numerator <- ident[which(ident == cluster)]
                assert(all(numerator == cluster))
                numerator <- sort(names(numerator))
                assert(hasLength(numerator))
                ## Denominator: cells in all other clusters.
                denominator <- sort(setdiff(colnames(object), numerator))
                assert(hasLength(denominator))
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
    definition = `findMarkers,SCE`
)
