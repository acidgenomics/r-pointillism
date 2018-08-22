#' Plot Known Markers Detected
#'
#' @name plotKnownMarkersDetected
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param markers `grouped_df`. Marker genes, grouped by "`cellType`".
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' plotKnownMarkersDetected(
#'     object = sce_small,
#'     markers = known_markers_small
#' )
NULL



#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature("SingleCellExperiment"),
    function(
        object,
        markers,
        reducedDim = c("TSNE", "UMAP"),
        headerLevel = 2L,
        ...
    ) {
        .assertIsKnownMarkersDetected(markers)
        reducedDim <- match.arg(reducedDim)
        assertIsAHeaderLevel(headerLevel)

        cellTypes <- markers %>%
            pull("cellType") %>%
            as.character() %>%
            na.omit() %>%
            unique()
        assert_is_non_empty(cellTypes)

        list <- pblapply(cellTypes, function(cellType) {
            genes <- markers %>%
                filter(cellType == !!cellType) %>%
                pull("geneID") %>%
                as.character() %>%
                na.omit() %>%
                unique()
            assert_is_non_empty(genes)

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
                    reducedDim = reducedDim,
                    ...
                )
                show(p)
                invisible(p)
            })
        })

        invisible(list)
    }
)



#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature("seurat"),
    getMethod("plotKnownMarkersDetected", "SingleCellExperiment")
)
