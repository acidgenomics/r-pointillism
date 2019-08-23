#' @name plotKnownMarkers
#' @inherit bioverbs::plotKnownMarkers
#' @note Updated 2019-08-23.
#'
#' @inheritParams acidroxygen::params
#' @inheritParams BiocParallel::bplapply
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



## Updated 2019-08-23.
`plotKnownMarkers,SingleCellExperiment,KnownMarkers` <-  # nolint
    function(
        object,
        markers,
        reduction,
        headerLevel,
        BPPARAM = BiocParallel::SerialParam(),  # nolint
        ...
    ) {
        validObject(object)
        validObject(markers)
        assert(
            isSubset(unique(markers[["name"]]), rownames(object)),
            isScalar(reduction),
            isHeaderLevel(headerLevel)
        )
        ## Safe to remove our nested ranges.
        markers <- as(markers, "DataFrame")
        markers[["ranges"]] <- NULL
        cellTypes <- markers %>%
            .[["cellType"]] %>%
            as.character() %>%
            na.omit() %>%
            unique()
        assert(isNonEmpty(cellTypes))
        list <- bplapply(
            X = cellTypes,
            FUN = function(cellType) {
                genes <- markers %>%
                    as_tibble() %>%
                    filter(cellType == !!cellType) %>%
                    pull("name") %>%
                    as.character() %>%
                    na.omit() %>%
                    unique()
                assert(isNonEmpty(genes))
                markdownHeader(
                    text = as.character(cellType),
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
            },
            BPPARAM = BPPARAM
        )
        invisible(list)
    }

formals(`plotKnownMarkers,SingleCellExperiment,KnownMarkers`)[c(
    "headerLevel",
    "reduction"
)] <- list(
    headerLevel = headerLevel,
    reduction = reduction
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
