#' @name findMarkers
#' @inherit acidgenerics::findMarkers
#' @note Updated 2020-01-30.
#'
#' @inheritParams acidroxygen::params
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
#' @importFrom acidgenerics findMarkers
#' @usage findMarkers(object, ...)
#' @export
NULL



## Updated 2020-01-30.
`findMarkers,SingleCellExperiment` <-  # nolint
    function(object, ...) {
        cli_h1("{.fun findMarkers}")
        object <- as(object, "SingleCellExperiment")
        ## Get the cluster mappings. Following the Seurat nomenclature here of
        ## using "ident" to denote the cluster identifier mapping factor.
        ident <- clusters(object)
        assert(is.factor(ident), hasNames(ident))
        clusters <- levels(ident)
        assert(length(clusters) >= 2L)
        cli_alert_info(sprintf("%d clusters detected.", length(clusters)))
        ## Loop across the clusters and calculate gene enrichment relative to
        ## all of the other clusters combined.
        list <- lapply(
            X = clusters,
            FUN = function(cluster) {
                cli_h2(sprintf("Cluster %s", as.character(cluster)))
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
