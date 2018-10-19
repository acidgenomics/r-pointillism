# FIXME Improve the formals.
# FIXME Need to rework this function to use class.
# FIXME Improve error handling for this:
# markers <- head(all_markers_small, n = 1L)



#' Plot Top Markers
#'
#' The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @name plotTopMarkers
#' @include globals.R
#'
#' @inheritParams general
#' @inheritParams topMarkers
#' @param markers `grouped_df`. Marker genes, grouped by "`cluster`".
#'
#' @return Show graphical output. Invisibly return `ggplot` `list`.
#'
#' @examples
#' data(seurat_small, all_markers_small)
#' plotTopMarkers(object = seurat_small, markers = all_markers_small)
NULL



.plotTopMarkers.SingleCellExperiment <-  # nolint
    function(
        object,
        markers,
        n = 1L,
        direction,
        reducedDim,
        headerLevel = 2L,
        ...
    ) {
        # Passthrough: n, direction, coding
        validObject(object)
        markers <- topMarkers(
            object = markers,
            n = n,
            direction = direction
        )
        assert_is_scalar(reducedDim)
        assertIsHeaderLevel(headerLevel)

        assert_is_subset("cluster", colnames(markers))
        clusters <- levels(markers[["cluster"]])

        # FIXME Add progress option.
        list <- pblapply(clusters, function(cluster) {
            genes <- markers %>%
                ungroup() %>%
                filter(cluster == !!cluster) %>%
                pull("name")
            if (!length(genes)) {
                message(paste0("No genes for cluster ", cluster, "."))
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
formals(.plotTopMarkers.SingleCellExperiment)[["direction"]] <- direction
formals(.plotTopMarkers.SingleCellExperiment)[["reducedDim"]] <- reducedDim



#' @rdname plotTopMarkers
#' @export
setMethod(
    f = "plotTopMarkers",
    signature = signature("SingleCellExperiment"),
    definition = .plotTopMarkers.SingleCellExperiment
)



#' @rdname plotTopMarkers
#' @export
setMethod(
    f = "plotTopMarkers",
    signature = signature("seurat"),
    definition = getMethod(
        f = "plotTopMarkers",
        signature = signature("SingleCellExperiment")
    )
)
