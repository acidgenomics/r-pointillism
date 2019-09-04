#' @name plotKnownMarkers
#' @inherit bioverbs::plotKnownMarkers
#' @note Updated 2019-09-03.
#'
#' @inheritParams acidroxygen::params
#' @param markers Object.
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @examples
#' data(Seurat, package = "acidtest")
#' data(seuratKnownMarkers)
#'
#' ## Seurat ====
#' object <- Seurat
#' markers <- seuratKnownMarkers
#' plotKnownMarkers(
#'     object = object,
#'     markers = markers,
#'     reduction = "UMAP"
#' )
NULL



#' @rdname plotKnownMarkers
#' @name plotKnownMarkers
#' @importFrom bioverbs plotKnownMarkers
#' @usage plotKnownMarkers(object, markers, ...)
#' @export
NULL



## Updated 2019-09-03.
`plotKnownMarkers,SingleCellExperiment,KnownMarkers` <-  # nolint
    function(
        object,
        markers,
        reduction,
        headerLevel,
        ...
    ) {
        validObject(object)
        validObject(markers)
        assert(
            isSubset(unique(markers[["name"]]), rownames(object)),
            isScalar(reduction),
            isHeaderLevel(headerLevel)
        )
        markers <- as(markers, "DataFrame")
        cellTypes <- markers[["cellType"]]
        cellTypes <- unique(na.omit(as.character(cellTypes)))
        assert(isNonEmpty(cellTypes))
        list <- lapply(
            X = cellTypes,
            FUN = function(cellType) {
                genes <- markers[
                    markers[["cellType"]] == cellType,
                    "name",
                    drop = TRUE
                ]
                genes <- unique(na.omit(as.character(genes)))
                assert(isNonEmpty(genes))
                markdownHeader(
                    text = cellType,
                    level = headerLevel,
                    tabset = TRUE,
                    asis = TRUE
                )
                lapply(genes, function(gene) {
                    p <- plotMarker(
                        object = object,
                        genes = gene,
                        reduction = reduction,
                        ...
                    )
                    show(p)
                    invisible(p)
                })
            }
        )
        invisible(list)
    }

formals(`plotKnownMarkers,SingleCellExperiment,KnownMarkers`)[c(
    "headerLevel",
    "reduction",
    "BPPARAM"
)] <- list(
    headerLevel = headerLevel,
    reduction = reduction,
    BPPARAM = BPPARAM
)



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature(
        object = "SingleCellExperiment",
        markers = "KnownMarkers"
    ),
    definition = `plotKnownMarkers,SingleCellExperiment,KnownMarkers`
)



## Updated 2019-07-31.
`plotKnownMarkers,Seurat,KnownMarkers` <-  # nolint
    `plotKnownMarkers,SingleCellExperiment,KnownMarkers`



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature(
        object = "Seurat",
        markers = "KnownMarkers"
    ),
    definition = `plotKnownMarkers,Seurat,KnownMarkers`
)
