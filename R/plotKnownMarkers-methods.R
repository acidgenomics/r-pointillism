#' @name plotKnownMarkers
#' @inherit acidgenerics::plotKnownMarkers
#' @note Updated 2020-01-30.
#'
#' @inheritParams acidroxygen::params
#' @param markers Object.
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @examples
#' data(Seurat, package = "acidtest")
#' data(seurat_known_markers)
#'
#' ## Seurat ====
#' object <- Seurat
#' markers <- seurat_known_markers
#' plotKnownMarkers(
#'     object = object,
#'     markers = markers,
#'     reduction = "UMAP"
#' )
NULL



#' @rdname plotKnownMarkers
#' @name plotKnownMarkers
#' @importFrom acidgenerics plotKnownMarkers
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
        assert(hasLength(cellTypes))
        list <- lapply(
            X = cellTypes,
            FUN = function(cellType) {
                genes <- markers[
                    markers[["cellType"]] == cellType,
                    "name",
                    drop = TRUE
                ]
                genes <- unique(na.omit(as.character(genes)))
                assert(hasLength(genes))
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

args <- c("headerLevel", "reduction", "BPPARAM")
formals(`plotKnownMarkers,SingleCellExperiment,KnownMarkers`)[args] <-
    .formalsList[args]
rm(args)



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
