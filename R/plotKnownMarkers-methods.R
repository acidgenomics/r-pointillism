#' Plot Known Markers
#'
#' @name plotKnownMarkers
#'
#' @inheritParams general
#' @param markers `grouped_df`. Marker genes, grouped by "`cellType`".
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' data(seurat_small, known_markers_small)
#' plotKnownMarkers(
#'     object = seurat_small,
#'     markers = known_markers_small
#' )
NULL



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
        assert_is_subset(
            x = unique(markers[["name"]]),
            y = rownames(object)
        )
        assert_is_scalar(reducedDim)
        assertIsHeaderLevel(headerLevel)
        assert_is_a_bool(progress)
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
        assert_is_non_empty(cellTypes)

        list <- applyFun(cellTypes, function(cellType) {
            genes <- markers %>%
                as_tibble() %>%
                filter(cellType == !!cellType) %>%
                pull("name") %>%
                as.character() %>%
                na.omit() %>%
                unique()
            assert_is_non_empty(genes)

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
        markers = "KnownSeuratMarkers"
    ),
    definition = plotKnownMarkers.SingleCellExperiment
)



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature(
        object = "seurat",
        markers = "KnownSeuratMarkers"
    ),
    definition = getMethod(
        f = "plotKnownMarkers",
        signature = signature(
            object = "SingleCellExperiment",
            markers = "KnownSeuratMarkers"
        )
    )
)
