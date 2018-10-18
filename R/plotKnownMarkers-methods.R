# FIXME Fix the formals.
# FIXME Rethink how this works.



#' Plot Known Markers
#'
#' @name plotKnownMarkers
#' @family Plot Functions
#'
#' @inheritParams general
#' @param markers `grouped_df`. Marker genes, grouped by "`cellType`".
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



.plotKnownMarkers.SCE <-  # nolint
    function(
        object,
        markers,
        reducedDim = 1L,
        headerLevel = 2L,
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

        # Safe to remove our nested ranges.
        markers <- as(markers, "DataFrame")
        markers[["ranges"]] <- NULL

        cellTypes <- markers %>%
            .[["cellType"]] %>%
            as.character() %>%
            na.omit() %>%
            unique()
        assert_is_non_empty(cellTypes)

        list <- pblapply(cellTypes, function(cellType) {
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
                    reducedDim = reducedDim
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
    signature = signature(
        object = "SingleCellExperiment",
        markers = "KnownSeuratMarkers"
    ),
    definition = .plotKnownMarkers.SCE
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
