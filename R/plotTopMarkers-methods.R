#' Plot Top Markers
#'
#' The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @name plotTopMarkers
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams topMarkers
#' @param markers `grouped_df`. Marker genes, grouped by "`cluster`".
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' markers <- topMarkers(all_markers_small, n = 1)
#' glimpse(markers)
#' plotTopMarkers(seurat_small, markers = tail(markers, 1))
NULL



# Methods ======================================================================
#' @rdname plotTopMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature("SingleCellExperiment"),
    function(
        object,
        markers,
        n = 10L,
        direction = c("positive", "negative", "both"),
        coding = FALSE,
        reducedDim = c("TSNE", "UMAP"),
        headerLevel = 2L,
        ...
    ) {
        # Passthrough: n, direction, coding
        validObject(object)
        stopifnot(is(markers, "grouped_df"))
        stopifnot(.isSanitizedMarkers(markers))
        markers <- topMarkers(
            data = markers,
            n = n,
            direction = direction,
            coding = coding
        )
        reducedDim <- match.arg(reducedDim)
        assertIsAHeaderLevel(headerLevel)

        assert_is_subset("cluster", colnames(markers))
        clusters <- levels(markers[["cluster"]])

        list <- pblapply(clusters, function(cluster) {
            genes <- markers %>%
                filter(cluster == !!cluster) %>%
                pull("rowname")
            if (!length(genes)) {
                return(invisible())
            }
            if (length(genes) > 10L) {
                warning("Maximum of 10 genes per cluster is recommended")
            }

            markdownHeader(
                text = paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE
            )

            lapply(genes, function(gene) {
                markdownHeader(
                    text = gene,
                    level = headerLevel + 1L,
                    asis = TRUE
                )
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
