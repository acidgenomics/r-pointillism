#' @name plotKnownMarkers
#' @inherit bioverbs::plotKnownMarkers
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @param markers `grouped_df`.
#'   Marker genes, grouped by `"cellType"`.
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @examples
#' data(seurat_small, known_markers_small)
#' plotKnownMarkers(
#'     object = seurat_small,
#'     markers = known_markers_small
#' )
NULL



#' @importFrom bioverbs plotKnownMarkers
#' @aliases NULL
#' @export
bioverbs::plotKnownMarkers



plotKnownMarkers.SingleCellExperiment <-  # nolint
    function(
        object,
        markers,
        reducedDim,
        headerLevel,
        progress = FALSE,
        ...
    ) {
        validObject(object)
        validObject(markers)
        assert(
            isSubset(unique(markers[["name"]]), rownames(object)),
            isScalar(reducedDim),
            isHeaderLevel(headerLevel),
            isFlag(progress)
        )
        if (isTRUE(progress)) {
            applyFun <- pblapply
        } else {
            applyFun <- lapply
        }

        # Safe to remove our nested ranges.
        markers <- as(markers, "DataFrame")
        markers[["ranges"]] <- NULL

        cellTypes <- markers %>%
            .[["cellType"]] %>%
            as.character() %>%
            na.omit() %>%
            unique()
        assert(isNonEmpty(cellTypes))

        list <- applyFun(cellTypes, function(cellType) {
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
                    reducedDim = reducedDim,
                    ...
                )
                show(p)
                invisible(p)
            })
        })

        invisible(list)
    }

formals(plotKnownMarkers.SingleCellExperiment)[c(
    "headerLevel",
    "reducedDim"
)] <- list(
    headerLevel = headerLevel,
    reducedDim = reducedDim
)



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature(
        object = "SingleCellExperiment",
        markers = "KnownMarkers"
    ),
    definition = plotKnownMarkers.SingleCellExperiment
)



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature(
        object = "seurat",
        markers = "KnownMarkers"
    ),
    definition = plotKnownMarkers.SingleCellExperiment
)
