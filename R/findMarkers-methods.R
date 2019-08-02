#' @name findMarkers
#' @inherit bioverbs::findMarkers
#' @note Updated 2019-08-02.
#'
#' @inheritParams basejump::params
#' @param ... Passthrough arguments to [diffExp()].
#'
#' @return `list` containing:
#' - `caller = "edgeR"`: `DGELRT`.
#' - `caller = "DESeq2"`: `DESeqResults`.
#'
#' @examples
#' data(Seurat, package = "acidtest")
#'
#' ## Seurat ====
#' object <- Seurat
#' x <- findMarkers(object, caller = "edgeR")
#' class(x)
#' lapply(x, class)
NULL



#' @rdname findMarkers
#' @name findMarkers
#' @importFrom bioverbs findMarkers
#' @usage findMarkers(object, ...)
#' @export
NULL



## Updated 2019-08-02.
`findMarkers,SingleCellExperiment` <-  # nolint
    function(object, ...) {
        object <- as(object, "SingleCellExperiment")

        ## Get the cluster mappings. Following the Seurat nomenclature here of
        ## using "ident" to denote the cluster identifier mapping factor.
        ident <- clusters(object)
        assert(is.factor(ident), hasNames(ident))
        clusters <- levels(ident)
        assert(length(clusters) >= 2L)
        message(paste(length(clusters), "clusters detected"))

        ## Loop across the clusters and calculate gene enrichment relative to
        ## all of the other clusters combined.
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                message(paste("Cluster", cluster, "===="))
                ## Numerator: cells in the current cluster.
                numerator <- ident[which(ident == cluster)]
                assert(all(numerator == cluster))
                numerator <- sort(names(numerator))
                assert(isNonEmpty(numerator))
                ## Denominator: cells in all other clusters.
                denominator <- sort(setdiff(colnames(object), numerator))
                assert(isNonEmpty(denominator))
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
    definition = `findMarkers,SingleCellExperiment`
)



## Updated 2019-07-31.
`findMarkers,Seurat` <-  # nolint
    `findMarkers,SingleCellExperiment`



#' @rdname findMarkers
#' @export
setMethod(
    f = "findMarkers",
    signature = signature("Seurat"),
    definition = `findMarkers,Seurat`
)
