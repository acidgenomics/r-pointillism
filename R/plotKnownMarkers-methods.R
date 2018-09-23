# FIXME Fix the formals.



#' Plot Known Markers
#'
#' @name plotKnownMarkers
#' @family Plot Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param markers `grouped_df`. Marker genes, grouped by "`cellType`".
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' plotKnownMarkers(
#'     object = sce_small,
#'     markers = known_markers_small
#' )
NULL



# FIXME Rethink how this works.
.plotKnownMarkers.SCE <-  # nolint
    function(
        object,
        markers,
        reducedDim = c("TSNE", "UMAP"),
        headerLevel = 2L,
        ...
    ) {
        assert_is_all_of(markers, "CellTypeMarkers")
        reducedDim <- match.arg(reducedDim)
        assertIsHeaderLevel(headerLevel)

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



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature("SingleCellExperiment"),
    definition = .plotKnownMarkers.SCE
)



#' @rdname plotKnownMarkers
#' @export
setMethod(
    f = "plotKnownMarkers",
    signature = signature("seurat"),
    definition = getMethod(
        f = "plotKnownMarkers",
        signature = signature("SingleCellExperiment")
    )
)
