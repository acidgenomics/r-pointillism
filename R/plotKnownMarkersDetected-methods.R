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
#' markers <- head(known_markers_small, n = 1)
#' glimpse(markers)
#' plotKnownMarkersDetected(seurat_small, markers = head(markers, 1))
NULL



# Methods ======================================================================
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
        assert_has_rows(markers)
        stopifnot(is(markers, "grouped_df"))
        assert_has_rows(markers)
        assert_is_subset("cellType", colnames(markers))
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
