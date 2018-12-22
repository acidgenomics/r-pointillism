#' @name plotTopMarkers
#' @include globals.R
#' @inherit bioverbs::plotTopMarkers
#' @inheritParams basejump::params
#' @inheritParams topMarkers
#'
#' @details
#' The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @param markers `grouped_df`.
#'   Marker genes, grouped by "`cluster`".
#' @param ... Passthrough arguments to [plotMarker()].
#'
#' @examples
#' data(seurat_small, all_markers_small)
#' plotTopMarkers(object = seurat_small, markers = all_markers_small)
NULL



#' @importFrom bioverbs plotTopMarkers
#' @aliases NULL
#' @export
bioverbs::plotTopMarkers



plotTopMarkers.SingleCellExperiment <-  # nolint
    function(
        object,
        markers,
        n = 1L,
        direction,
        reducedDim,
        headerLevel = 2L,
        progress = FALSE,
        ...
    ) {
        # Passthrough: n, direction, coding
        validObject(object)
        markers <- topMarkers(
            object = markers,
            n = n,
            direction = direction
        )
        assert(
            isScalar(reducedDim),
            isHeaderLevel(headerLevel),
            isFlag(progress)
        )
        if (isTRUE(progress)) {
            applyFun <- pblapply
        } else {
            applyFun <- lapply
        }

        assert(isSubset("cluster", colnames(markers)))
        clusters <- levels(markers[["cluster"]])


        list <- applyFun(clusters, function(cluster) {
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
formals(plotTopMarkers.SingleCellExperiment)[["direction"]] <- direction
formals(plotTopMarkers.SingleCellExperiment)[["reducedDim"]] <- reducedDim



#' @rdname plotTopMarkers
#' @export
setMethod(
    f = "plotTopMarkers",
    signature = signature("SingleCellExperiment"),
    definition = plotTopMarkers.SingleCellExperiment
)



#' @rdname plotTopMarkers
#' @export
setMethod(
    f = "plotTopMarkers",
    signature = signature("seurat"),
    definition = plotTopMarkers.SingleCellExperiment
)
